"""
Usage:
    zot vars [-r ref] <input>...

Options:
    -r ref      reference k-mers
"""

from zotmer.library.basics import can, fasta, lcp, render
from zotmer.library.file import readFasta
from zotmer.library.misc import heap
from zotmer.library.stats import logBinGe, logBinLe
from zotmer.library.exceptions import MismatchedK
from zotmer.library.kmers import kmers
from zotmer.library.files import readKmersAndCounts


import docopt
import math
import sys

def getK(ins):
    k = None
    for fn in ins:
        with kmers(fn, 'r') as z:
            k0 = z.meta['K']
            if k is None:
                k = k0
            elif k != k0:
                raise MismatchedK(k, k0)
    return k

def group(K, J, n, xs):
    py = 0
    grp = []
    S = 2*(K-J)
    for (x,c) in xs:
        y = x >> S
        if y != py:
            if len(grp) > 0:
                yield (py, n, grp)
            py = y
            grp = []
        grp.append((x,c))
    if len(grp) > 0:
        yield (py, n, grp)

class Group:
    def __init__(self, K, J, n, xs):
        self.grps = group(K, J, n, xs)
        self.curr = None
        self.more = True
        self.next()

    def valid(self):
        return self.more

    def next(self):
        try:
            self.curr = self.grps.next()
        except StopIteration:
            self.more = False

    def this(self):
        if not self.valid():
            print self
        assert self.valid()
        return self.curr

    def __lt__(self, other):
        assert self.valid()
        assert other.valid()
        return self.this()[0] < other.this()[0]

    def __str__(self):
        return repr(self.curr) + '\t' + repr(self.more)

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = getK(opts['<input>'])
    J = K - 1
    M = (1 << (2 * (K - J))) - 1

    if opts['-r'] is not None:
        with kmers(opts['-r'], 'r') as z:
            xs = list(group(K, J, 0, readKmersAndCounts(z)))

        for fn in opts['<input>']:
            with kmers(fn, 'r') as z:
                samXs = readKmersAndCounts(z)
                i = 0
                for (yCtx, _, yGrp) in group(K, J, 0, samXs):
                    while i < len(xs) and xs[i][0] < yCtx:
                        i += 1
                    assert i < len(xs)
                    assert xs[i][0] == yCtx
                    gt = float(sum([c for (x,c) in xs[i][2]]))
                    gx = [0 for j in xrange(M+1)]
                    for (x,c) in xs[i][2]:
                        gx[x&M] = c
                    st = sum([c for (x,c) in yGrp])
                    sx = [0 for j in xrange(M+1)]
                    for (x,c) in yGrp:
                        sx[x&M] = c
                    ss = []
                    b = 0
                    for j in xrange(M+1):
                        p = float(gx[j])/gt
                        v = 0.0
                        if 0.0 < p and p < 1.0:
                            v = logBinGe(p, st, sx[j])
                            if v < -10:
                                b |= 1 << j
                        ss.append('%3.2g' % (v,))
                    if b > 0:
                        print '%s\t%s\t%s' % (render(J, yCtx), fasta(b), '\t'.join(ss))
                    i += 1
        return

    # Parse files in parallel to get global distribution

    N = len(opts['<input>'])
    h = heap.heap()
    i = 0
    for fn in opts['<input>']:
        (_, xs) = kfset.read(fn)
        i += 1
        h.push(Group(K, J, i, xs))

    while len(h) > 0:
        xfs = []
        g = h.pop()
        gy = g.this()[0]
        xfs.append(g.this())
        g.next()
        if g.valid():
            h.push(g)
        for x in h.xs:
            assert x.valid()
        while len(h) > 0 and h.front().this()[0] == gy:
            g = h.pop()
            xfs.append(g.this())
            g.next()
            if g.valid():
                h.push(g)
            for i in xrange(len(h.xs)):
                assert h.xs[i].valid()

        ds = []
        gc = [0 for i in xrange(M+1)]
        for (_, n, xc) in xfs:
            t = sum([c for (x,c) in xc])
            d = [0 for i in xrange(M+1)]
            for (x,c) in xc:
                j = x & M
                gc[j] += c
                d[j] = c
            ds.append((n, d))

        res = ['*' for i in xrange(N)]
        seen = set([])
        gt = float(sum(gc))
        for (n, d) in ds:
            t = sum(d)
            b = [0 for i in xrange((M+1)/4)]
            for i in xrange(M+1):
                p = float(gc[i])/gt
                if 0.0 < p and p < 1.0:
                    #vL = logBinLe(p, t, d[i])
                    #vG = logBinGe(p, t, d[i])
                    #v = min(vL, vG)
                    v = logBinGe(p, t, d[i])
                    if v > -10:
                        w = i >> 2
                        j = i & 3
                        b[w] |= 1 << j
            res[n-1] = ''.join([fasta(b0) for b0 in b])
            seen.add(res[n-1])
        if len(seen) > 1:
            print '%s\t%s' % (render(J, gy), '\t'.join(res))

if __name__ == '__main__':
    main(sys.argv[1:])
