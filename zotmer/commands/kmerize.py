"""
Usage:
    zot kmerize [options] <k> <output> <input>...

Kmerize FASTA or FASTQ inputs to produce a standard container object.

Arguments:
    <k>         the length of the k-mers. Recommended values: 10-30
    <output>    the name of the output file.
                recommended naming convention
                    - mykmers.k25 for a k-mer set of 25-mers
                    - mykmers.kf25 for a k-mer frequency set of 25-mers
                    - mykmers.e25 for an expanded k-mer set of 25-mers

Options:
    -m MEM      in-memory buffer size (in MB)
    -C BAITS    capture mode - use kmers from the given FASTA file.
    -D FRAC     subsample k-mers, using FRAC proportion of k-mers
    -S SEED     if -D is given, give a seed for determining the
                subspace (defaults to 0).
    -v          produce verbose progress messages
"""

import array
import os
import sys
import time

import docopt

from pykmer.basics import kmersList, render, sub
from pykmer.container.casket import casket
from pykmer.file import openFile, readFasta, readFastqBlock, tmpfile
from pykmer.misc import heap, radix_sort
from pykmer.timer import timer

from zotmer.library.files import countsWriter, kmerWriter, readKmersAndCounts, writeKmersAndCounts, writeKmersAndCounts2
import zotmer.library.kmers as zotk
from zotmer.library.reads import reads

def merge(xs, cs, ys, zs, ss):

    xZ = len(xs)
    yZ = len(ys)

    i = 0
    haveX = True
    j = 0
    haveY = True

    if i < xZ:
        x = xs[i]
        xc = cs[i]
        i += 1
    else:
        haveX = False

    if j < yZ:
        y = ys[j]
        yc = 0
        while j < yZ and ys[j] == y:
            yc += 1
            j += 1
    else:
        haveY = False

    while haveX and haveY:
        if x < y:
            zs.append(x)
            ss.append(xc)
            if i < xZ:
                x = xs[i]
                xc = cs[i]
                i += 1
            else:
                haveX = False
            continue

        if x > y:
            zs.append(y)
            ss.append(yc)
            if j < yZ:
                y = ys[j]
                yc = 0
                while j < yZ and ys[j] == y:
                    yc += 1
                    j += 1
            else:
                haveY = False
            continue

        # x == y
        zs.append(x)
        ss.append(xc + yc)

        if i < xZ:
            x = xs[i]
            xc = cs[i]
            i += 1
        else:
            haveX = False

        if j < yZ:
            y = ys[j]
            yc = 0
            while j < yZ and ys[j] == y:
                yc += 1
                j += 1
        else:
            haveY = False

    while haveX:
        zs.append(x)
        ss.append(xc)
        if i < xZ:
            x = xs[i]
            xc = cs[i]
            i += 1
        else:
            haveX = False

    while haveY:
        zs.append(y)
        ss.append(yc)
        if j < yZ:
            y = ys[j]
            yc = 0
            while j < yZ and ys[j] == y:
                yc += 1
                j += 1
        else:
            haveY = False

def merge2(xs, ys, h):
    moreXs = True
    try:
        x = xs.next()
    except StopIteration:
        moreXs = False
    moreYs = True
    try:
        y = ys.next()
    except StopIteration:
        moreYs = False

    while moreXs and moreYs:
        if x[0] < y[0]:
            yield x
            h[x[1]] = 1 + h.get(x[1], 0)
            try:
                x = xs.next()
            except StopIteration:
                moreXs = False
            continue
        if x[0] > y[0]:
            yield y
            h[y[1]] = 1 + h.get(y[1], 0)
            try:
                y = ys.next()
            except StopIteration:
                moreYs = False
            continue
        v = (x[0], x[1] + y[1])
        yield v
        h[v[1]] = 1 + h.get(v[1], 0)
        try:
            x = xs.next()
        except StopIteration:
            moreXs = False
        try:
            y = ys.next()
        except StopIteration:
            moreYs = False

    while moreXs:
        yield x
        h[x[1]] = 1 + h.get(x[1], 0)
        try:
            x = xs.next()
        except StopIteration:
            moreXs = False

    while moreYs:
        yield y
        h[y[1]] = 1 + h.get(y[1], 0)
        try:
            y = ys.next()
        except StopIteration:
            moreYs = False

class _kmerStream(object):
    def __init__(self, xs):
        self.more = True
        self.xs = xs
        self.x = None
        self.next()

    def next(self):
        try:
            self.x = self.xs.next()
        except StopIteration:
            self.more = False

    def __getitem__(self, i):
        return self.x[i]

    def __lt__(self, other):
        assert self.more
        assert other.more
        return self.x[0] < other.x[0]

def mergeN(xss, hist):
    ss = [_kmerStream(xs) for xs in xss]
    q = heap(ss)
    while len(q) > 0:
        v = q.front()
        x = v.x[0]
        c = 0
        n = 0
        while len(q) > 0 and v.x[0] == x:
            c += v.x[1]
            n += 1
            v.next()
            if v.more:
                q.modifyfront()
                v = q.front()
            else:
                q.pop()
                if len(q) > 0:
                    v = q.front()
        if c not in hist:
            hist[c] = 0
        hist[c] += 1
        yield (x, c)

class _kmerRadixBlockStream(object):
    def __init__(self, K, stream):
        self.B = 12
        self.S = 2*K - self.B
        self.more = True
        self.stream = stream
        self.radix = 0
        self.xs = []
        self.x = None
        try:
            self.x = self.stream.next()
        except StopIteration:
            self.more = False
        self.next()

    def done(self):
        return len(self.xs) == 0 and not self.more

    def next(self):
        self.xs = []
        while self.more and (self.x[0] >> self.S) == self.radix:
            self.xs.append(self.x)
            try:
                self.x = self.stream.next()
            except StopIteration:
                self.more = False
                self.x = None
                break
        self.radix += 1

    def __lt__(self, other):
        return self.radix < other.radix

def mergeNinto(K, xss, hist, z, nm = None):
    if nm is None:
        knm = 'kmers'
        cnm = 'counts'
    else:
        knm = nm + '-kmers'
        cnm = nm + '-counts'

    t = tmpfile()
    with z.add_stream(knm) as kf, open(t, 'w') as cf:
        with kmerWriter(kf) as kx, countsWriter(cf) as cx:
            ss = [_kmerRadixBlockStream(K, xs) for xs in xss]
            q = heap(ss)
            while len(q) > 0:
                v = q.front()
                r = v.radix
                ys = {}
                while len(q) > 0 and v.radix == r:
                    for (x, c) in v.xs:
                        ys[x] = c + ys.get(x, 0)
                    v.next()
                    if v.done():
                        q.pop()
                        if len(q) > 0:
                            v = q.front()
                    else:
                        q.modifyfront()
                        v = q.front()
                ys = ys.items()
                ys.sort()
                for (x, c) in ys:
                    hist[c] = 1 + hist.get(c, 0)
                    kx.append(x)
                    cx.append(c)
    z.add_file(cnm, t)
    os.remove(t)

class KmerAccumulator:
    def __init__(self):
        self.toc = {}
        self.z = 0

    def __len__(self):
        return self.z

    def clear(self):
        self.toc = {}
        self.z = 0

    def add(self, x):
        xh = x >> 32
        xl = x & 0xFFFFFFFF
        if xh not in self.toc:
            self.toc[xh] = array.array('I')
        self.toc[xh].append(xl)
        self.z += 1

    def addList(self, xs):
        for x in xs:
            xh = x >> 32
            xl = x & 0xFFFFFFFF
            if xh not in self.toc:
                self.toc[xh] = array.array('I')
            self.toc[xh].append(xl)
            self.z += 1

    def kmers(self):
        xhs = self.toc.keys()
        xhs.sort()
        for xh in xhs:
            x0 = xh << 32
            xls = self.toc[xh].tolist()
            xls.sort()
            for xl in xls:
                x = x0 | xl
                yield x

    def stat(self, brief=False):
        if brief == True:
            n = 0
            s = 0
            s2 = 0
            for (xh, xls) in self.toc.iteritems():
                l = len(xls)
                n += 1
                s += l
                s2 += l*l
            a = float(s) / float(n)
            v = float(s2) / float(n) - a*a
            print '%d\t%f\t%f' % (n, a, v)
            return
        h = {}
        for (xh, xls) in self.toc.iteritems():
            l = len(xls)
            h[l] = 1 + h.get(l, 0)

        h = h.items()
        h.sort()
        for (l,c) in h:
            print '%d\t%d' % (l, c)

class KmerAccumulator2:
    def __init__(self, K):
        self.K = K
        self.idxX = array.array('L', [])
        self.idxC = array.array('I', [])
        self.buf = []
        self.z0 = 0
        self.z = 0

    def __len__(self):
        return self.z

    def mem(self):
        return 8*len(self.idxX) + 4*len(self.idxC) + 8*len(self.buf)

    def clear(self):
        self.idxX = array.array('L', [])
        self.idxC = array.array('I', [])
        self.buf = []
        self.z0 = 0
        self.z = 0

    def add(self, x):
        self.z += 1

        self.buf.append(x)
        self.z0 += 1

        if self.z0 > len(self.idxX) and self.z0 > 128*1024*1024:
            self.flush()

    def addList(self, xs):
        z0 = len(self.buf)
        self.buf.extend(xs)
        z1 = len(self.buf)
        n = z1 - z0
        self.z += n
        self.z0 += n

        if self.z0 > len(self.idxX) and self.z0 > 128*1024*1024:
            self.flush()

    def flush(self):
        if len(self.buf) == 0:
            return

        #self.buf.sort()
        radix_sort(2*self.K, self.buf)
        ys = array.array('L', [])
        cs = array.array('I', [])
        merge(self.idxX, self.idxC, self.buf, ys, cs)
        self.idxX = ys
        self.idxC = cs
        self.buf = []
        self.z0 = 0

    def kmers(self):
        self.flush()
        for i in xrange(len(self.idxX)):
            yield (self.idxX[i], self.idxC[i])

    def kmersOnly(self):
        self.flush()
        return self.idxX

    def countsOnly(self):
        self.flush()
        return self.idxC

def hist(zs, h):
    for z in zs:
        h[z[1]] = 1 + h.get(z[1], 0)
        yield z

def histBlock(zss, h):
    for zs in zss:
        for c in zs[1]:
            h[c] = 1 + h.get(c, 0)
        yield zs

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    verbose = opts['-v']

    K = int(opts['<k>'])

    out = opts['<output>']

    Z = 1024*1024*32
    if opts['-m'] is not None:
        Z = 1024*1024*int(opts['-m'])

    buf = KmerAccumulator2(K)
    n = 0
    tmps = []
    acgt = [0, 0, 0, 0]
    m = 0

    d = None
    if opts['-D'] is not None:
        d = float(opts['-D'])

        S = 0
        if opts['-S'] is not None:
            S = int(opts['-S'])

        cacheYes = set([])
        cacheNo = set([])

    B = opts['-C']
    if B is not None:
        xs = set([])
        for (nm, seq) in readFasta(openFile(B)):
            xs |= set(kmersList(K, seq, True))
        B = xs

    tmpnm = tmpfile('.pmc')
    with casket(tmpnm, 'w') as z:
        nr = 0
        for itm in reads(opts['<input>'], K=K, pairs=False, reads=False, kmers=True, both=True, verbose=verbose):
            xs = itm.kmers[0]
            for x in xs:
                acgt[x&3] += 1
            if d is not None:
                for x in xs:
                    if x in cacheNo:
                        continue
                    if x not in cacheYes:
                        if not sub(S, d, x):
                            cacheNo.add(x)
                            continue
                        cacheYes.add(x)
                    buf.add(x)
                    m += 1
                    n += 1
                if len(cacheYes) > 1000000:
                    cacheYes = set([])
                if len(cacheNo) > 1000000:
                    cacheNo = set([])
            elif B is not None:
                found = False
                for x in xs:
                    if x in B:
                        found = True
                        break
                if found:
                    buf.addList(xs)
                    for x in xs:
                        m += 1
                        n += 1
            else:
                buf.addList(xs)
                for x in xs:
                    m += 1
                    n += 1

            nr += 1
            if (nr & 1023) == 0 and buf.mem() >= Z//2:
                fn = 'tmps-%d' % (len(tmps),)
                tmps.append(fn)
                writeKmersAndCounts2(z, buf.kmersOnly(), buf.countsOnly(), fn)
                buf.clear()
                n = 0

        if len(tmps) and len(buf):
            fn = 'tmps-%d' % (len(tmps),)
            tmps.append(fn)
            writeKmersAndCounts2(z, buf.kmersOnly(), buf.countsOnly(), fn)
            buf = []

    with zotk.kmers(out, 'w') as z:
        h = {}
        if len(tmps) == 0:
            for c in buf.countsOnly():
                h[c] = 1 + h.get(c, 0)
            writeKmersAndCounts2(z, buf.kmersOnly(), buf.countsOnly())
        elif len(tmps) == 1:
            with casket(tmpnm, 'r') as z0:
                writeKmersAndCounts(z, readKmersAndCounts(z0, tmps[0]))
        else:
            with casket(tmpnm, 'r') as z0:
                xss = [readKmersAndCounts(z0, t) for t in tmps]
                mergeNinto(K, xss, h, z)
        n = float(sum(acgt))
        acgt = [c/n for c in acgt]
        z.meta['K'] = K
        z.meta['kmers'] = 'kmers'
        z.meta['counts'] = 'counts'
        z.meta['hist'] = h
        z.meta['acgt'] = acgt
        z.meta['reads'] = nr
    os.remove(tmpnm)

if __name__ == '__main__':
    main(sys.argv[1:])
