"""
Usage:
    zot merge <output> <input>...
"""

import os
import sys

import docopt

from pykmer.basics import render
from pykmer.file import tmpfile
from pykmer.container.casket import casket
from pykmer.misc import heap
from zotmer.library.kmers import kmers
from zotmer.library.files import countsWriter, kmerWriter, readKmersAndCounts, writeKmersAndCounts, writeKmersAndCounts2

def pairs(xs):
    i = 0
    while i + 1 < len(xs):
        yield (xs[i], xs[i + 1])
        i += 2
    if i < len(xs):
        yield (xs[i], )

def merge(xs, ys):
    moreXs = True
    moreYs = True

    try:
        x = xs.next()
        xk = x[0]
    except StopIteration:
        moreXs = False
    try:
        y = ys.next()
        yk = y[0]
    except StopIteration:
        moreYs = False

    while moreXs and moreYs:
        if xk < yk:
            yield x
            try:
                x = xs.next()
                xk = x[0]
            except StopIteration:
                moreXs = False
            continue

        if xk > yk:
            yield y
            try:
                y = ys.next()
                yk = y[0]
            except StopIteration:
                moreYs = False
            continue

        assert xk == yk
        yield (xk, x[1] + y[1])

        try:
            x = xs.next()
            xk = x[0]
        except StopIteration:
            moreXs = False
        try:
            y = ys.next()
            yk = y[0]
        except StopIteration:
            moreYs = False

    while moreXs:
        yield x
        try:
            x = xs.next()
        except StopIteration:
            moreXs = False

    while moreYs:
        yield y
        try:
            y = ys.next()
        except StopIteration:
            moreYs = False

def hist(zs, h, acgt):
    for z in zs:
        h[z[1]] = 1 + h.get(z[1], 0)
        acgt[z[0]&3] += 1
        yield z

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

def mergeNinto(K, xss, hist, acgt, z, nm = None):
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
                    acgt[x&3] += c
                    kx.append(x)
                    cx.append(c)
    z.add_file(cnm, t)
    os.remove(t)

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = None

    out = opts['<output>']

    px = list(pairs(opts['<input>']))
    if len(px) == 1:
        with kmers(out, 'w') as z:
            h = {}
            acgt = [0, 0, 0, 0]
            ix = px[0]
            if len(ix) == 1:
                with kmers(ix[0], 'r') as z0:
                    K = z0.meta['K']
                    xs = readKmersAndCounts(z0)
                    zs = hist(xs, h, acgt)
                    writeKmersAndCounts(z, xs)
            else:
                with kmers(ix[0], 'r') as z0, kmers(ix[1], 'r') as z1:
                    K = z0.meta['K']
                    K1 = z1.meta['K']
                    if K1 != K:
                        print >> sys.stderr, "mismatched K"
                        sys.exit(1)
                    xs = readKmersAndCounts(z0)
                    ys = readKmersAndCounts(z1)
                    zs = hist(merge(xs, ys), h, acgt)
                    writeKmersAndCounts(z, zs)
            n = float(sum(acgt))
            acgt = [c/n for c in acgt]
            z.meta['hist'] = h
            z.meta['acgt'] = acgt
        return

    tmps = []
    tmpnm = tmpfile('.pmc')
    with casket(tmpnm, 'w') as z:
        for ix in px:
            if len(ix) == 1:
                nm = 'tmp-' + str(len(tmps))
                tmps.append(nm)
                with kmers(ix[0], 'r') as z0:
                    if K is None:
                        K = z0.meta['K']
                    else:
                        K0 = z0.meta['K']
                        if K0 != K:
                            print >> sys.stderr, "mismatched K"
                            sys.exit(1)
                    xs = readKmersAndCounts(z0)
                    writeKmersAndCounts(z, xs, nm)
            else:
                nm = 'tmp-' + str(len(tmps))
                tmps.append(nm)
                with kmers(ix[0], 'r') as z0, kmers(ix[1], 'r') as z1:
                    if K is None:
                        K = z0.meta['K']
                    else:
                        K0 = z0.meta['K']
                        if K0 != K:
                            print >> sys.stderr, "mismatched K"
                            sys.exit(1)
                        K1 = z1.meta['K']
                        if K1 != K:
                            print >> sys.stderr, "mismatched K"
                            sys.exit(1)
                    xs = readKmersAndCounts(z0)
                    ys = readKmersAndCounts(z1)
                    writeKmersAndCounts(z, merge(xs, ys), nm)

    assert K is not None

    with kmers(out, 'w') as z:
        h = {}
        acgt = [0, 0, 0, 0]
        with casket(tmpnm, 'r') as z0:
            xss = [readKmersAndCounts(z0, t) for t in tmps]
            mergeNinto(K, xss, h, acgt, z)
        n = float(sum(acgt))
        acgt = [c/n for c in acgt]
        z.meta['hist'] = h
        z.meta['acgt'] = acgt

    os.remove(tmpnm)

if __name__ == '__main__':
    main(sys.argv[1:])
