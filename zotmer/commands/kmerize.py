"""
Usage:
    zot kmerize [-m MEM] [-C BAITS] [-D FRAC [-S SEED]] <k> <output> <input>...

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
"""

import array
import gzip
import os
import subprocess
import sys
import time

import docopt

from pykmer.basics import kmersList, render, sub
from pykmer.container import container
from pykmer.container.std import readKmersAndCounts, readKmersAndCountsBlock, writeKmersAndCounts, writeKmersAndCountsBlock
from pykmer.file import openFile, readFasta, readFastqBlock, tmpfile

from zotmer.library.merge import merge2

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
        for x in xs:
            self.z += 1

            self.buf.append(x)
            self.z0 += 1

        if self.z0 > len(self.idxX) and self.z0 > 128*1024*1024:
            self.flush()

    def flush(self):
        self.buf.sort()
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

    def kmersBlock(self):
        self.flush()
        B = 65536
        for i in xrange(0, len(self.idxX), B):
            yield (self.idxX[i:i+B], self.idxC[i:i+B])

    def stat(self, brief=False):
        pass

def pairs(xs):
    i = 0
    res = []
    while i + 1 < len(xs):
        res.append((xs[i], xs[i + 1]))
        i += 2
    if i < len(xs):
        res.append((xs[i], ))
    return res

def stripCompressionSuffix(nm):
    if nm.endswith('.gz'):
        return nm[:-3]
    return nm

def isFasta(nm):
    bnm = stripCompressionSuffix(nm)
    if bnm.endswith(".fa"):
        return True
    if bnm.endswith(".fasta"):
        return True
    if bnm.endswith(".fas"):
        return True
    if bnm.endswith(".fna"):
        return True

def mkParser(fn):
    if isFasta(fn):
        for (nm, seq) in readFasta(openFile(fn)):
            yield [(nm, seq)]
    else:
        for grps in readFastqBlock(openFile(fn)):
            yield [(grp[0], grp[1]) for grp in grps]

def mkPairs(xs):
    m = 0
    p = 0
    n = 0
    for x in xs:
        if x != p:
            if n > 0:
                m += 1
                yield (p, n)
            p = x
            n = 0
        n += 1
    if n > 0:
        m += 1
        yield (p, n)

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
    with container(tmpnm, 'w') as z:
        pass

    PN = 1024*1024

    nr = 0
    t0 = time.time()
    for fn in opts['<input>']:
        for rds in mkParser(fn):
            for (nm, seq) in rds:
                nr += 1
                if nr & (PN - 1) == 0:
                    t1 = time.time()
                    print >> sys.stderr, 'reads processed:', nr, (PN)/(t1 - t0), 'reads/second'
                    t0 = t1
                    buf.stat()
                xs = kmersList(K, seq, True)
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
                        acgt[x&3] += 1
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
                            acgt[x&3] += 1
                            m += 1
                            n += 1
                else:
                    buf.addList(xs)
                    for x in xs:
                        acgt[x&3] += 1
                        m += 1
                        n += 1
                if (nr % 1023) == 0 and buf.mem() >= Z//2:
                    fn = 'tmps-%d' % (len(tmps),)
                    print >> sys.stderr, "writing " + fn + "\t" + tmpnm
                    tmps.append(fn)
                    with container(tmpnm, 'a') as z:
                        writeKmersAndCountsBlock(K, buf.kmersBlock(), z, fn)
                    buf.clear()
                    n = 0

    t1 = time.time()
    print >> sys.stderr, 'reads processed:', nr, (nr % PN)/(t1 - t0), 'reads/second'

    if len(tmps) and len(buf):
        fn = 'tmps-%d' % (len(tmps),)
        #print >> sys.stderr, "writing " + fn + "\t" + tmpnm
        tmps.append(fn)
        with container(tmpnm, 'a') as z:
            writeKmersAndCountsBlock(K, buf.kmersBlock(), z, fn)
        buf = []

    while len(tmps) > 2:
        tmpnm2 = tmpfile('.pmc')
        tmps2 = []
        with container(tmpnm, 'r') as z0, container(tmpnm2, 'w') as z:
            ps = pairs(tmps)
            for p in ps:
                fn = 'tmps-%d' % (len(tmps2),)
                tmps2.append(fn)
                if len(p) == 1:
                    writeKmersAndCountsBlock(K, readKmersAndCountsBlock(z0, p[0]), z, fn)
                    continue
                h = {}
                merge2(z, K, readKmersAndCounts(z0, p[0]), readKmersAndCounts(z0, p[1]), h, fn)
        os.remove(tmpnm)
        tmpnm = tmpnm2
        tmps = tmps2

    with container(out, 'w') as z:
        h = {}
        if len(tmps) == 0:
            zs = histBlock(buf.kmersBlock(), h)
            writeKmersAndCountsBlock(K, zs, z)
        elif len(tmps) == 1:
            with container(tmpnm, 'r') as z0:
                writeKmersAndCountsBlock(K, histBlock(readKmersAndCountsBlock(z0, tmps[0]), h), z)
        else:
            assert len(tmps) == 2
            with container(tmpnm, 'r') as z0:
                merge2(z, K, readKmersAndCounts(z0, tmps[0]), readKmersAndCounts(z0, tmps[1]), h)
        n = float(sum(acgt))
        acgt = [c/n for c in acgt]
        z.meta['hist'] = h
        z.meta['acgt'] = acgt
        z.meta['reads'] = nr
    os.remove(tmpnm)

if __name__ == '__main__':
    main(sys.argv[1:])
