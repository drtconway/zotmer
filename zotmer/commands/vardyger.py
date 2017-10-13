
"""
zot vardyger - search for tandem duplications

Usage:
    zot vardyger [options] <sequences> <input>...

Options:
    -k K            value of K to use - must be even. [default: 24]
    -s              strict path validation
    -v              produce verbose output
"""

import array
import math
import sys
import zipfile

import docopt
import tqdm
import yaml

from pykmer.basics import ham, kmersList, kmersWithPosList, rc, render
from pykmer.file import openFile, readFasta, readFastq
from zotmer.library.reads import reads
from zotmer.library.debruijn import interpolate, interpolateBetween, pathBetween

def pairs(xs):
    N = len(xs)
    for i in xrange(N):
        for j in xrange(i+1, N):
            yield (xs[j], xs[i])
            yield (xs[i], xs[j])

def findDup(refFst, refLst, xFst, xLst):
    for (c, ys) in xFst.iteritems():
        if c not in refFst:
            continue
        ds = refFst[c]
        for (b, d) in pairs(ys):
            if d not in ds:
                continue
            if b not in refLst:
                continue
            for a in refLst[b]:
                if a == c:
                    continue
                if a not in refFst:
                    continue
                yield (a, b, c, d)

def offset(xs, ys, d):
    for x in xs:
        for y in ys:
            if x + d == y:
                yield (x,y)

def positions(posIdx, J, a, b, c, d):
    pas = posIdx.get(a, [])
    pbs = posIdx.get(b, [])
    pabs = list(offset(pas, pbs, J))
    pcs = posIdx.get(c, [])
    pds = posIdx.get(d, [])
    pcds = list(offset(pcs, pds, J))

    res = None
    for (pa, pb) in pabs:
        for (pc, pd) in pcds:
            if pc <= pb:
                continue
            if res is None:
                res = []
            res.append((pa, pb, pc, pd))
    return res

def renderPath(K, xs):
    if len(xs) == 0:
        return ''
    res = [render(K, xs[0])]
    for x in xs[1:]:
        res.append('ACGT'[x&3])
    return ''.join(res)

def renderPath1(K, xs, X):
    J = K // 2
    #for i in range(len(xs)):
    #    print '%d\t%s\t%d' % (i, render(K, xs[i]), X.get(xs[i], 0))
    print '%d\t%s' % (len(xs) - 1, renderPath(K, xs)[J:-J])

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = int(opts['-k'])
    if (K & 1) != 0:
        print >> sys.stderr, "K must be even."
        return

    verbose = opts['-v']

    J = K // 2
    S = 2*(K - J)
    Mj = (1 << (2*J)) - 1

    names = []
    seqs = []
    bait = {}
    wtFst = []
    wtLst = []
    posIdx = []
    with openFile(opts['<sequences>']) as f:
        for (nm, seq) in readFasta(f):
            n = len(names)
            names.append(nm)
            seqs.append(seq)
            wf = {}
            wl = {}
            for x in kmersList(K, seq, False):
                if x not in bait:
                    bait[x] = set([])
                bait[x].add(n)

                y0 = x >> S
                y1 = x & Mj
                #print '- %s\t%s\t%s' % (render(K, x), render(J, y0), render(J, y1))

                if y0 not in wf:
                    wf[y0] = set([])
                wf[y0].add(y1)

                if y1 not in wl:
                    wl[y1] = set([])
                wl[y1].add(y0)

            wtFst.append(wf)
            wtLst.append(wl)
            
            px = {}
            for (x,p) in kmersWithPosList(J, seq, False):
                if x not in px:
                    px[x] = []
                px[x].append(p)
            posIdx.append(px)

    N = len(names)

    X = [{} for n in range(N)]
    for itm in reads(opts['<input>'], K=K, reads=False, kmers=True, both=True, verbose=verbose):
        xs = itm[0][0]
        hits = set([])
        for x in xs:
            if x in bait:
                hits |= bait[x]
        for n in hits:
            for x in xs:
                if x not in X[n]:
                    X[n][x] = 0
                X[n][x] += 1

    for n in range(N):
        xs = X[n]

        fst = {}
        lst = {}
        for (x,c) in xs.iteritems():
            if c < 5:
                continue
            y0 = x >> S
            y1 = x & Mj

            if y0 not in fst:
                fst[y0] = []
            fst[y0].append(y1)

            if y1 not in lst:
                lst[y1] = []
            lst[y1].append(y0)

        for (a, b, c, d) in findDup(wtFst[n], wtLst[n], fst, lst):
            pps = positions(posIdx[n], J, a, b, c, d)
            if pps is None:
                continue
            for pp in pps:
                ab = a << S | b
                cb = c << S | b
                cd = c << S | d

                dd = pp[2] - pp[0]

                if opts['-s']:
                    fstPath = pathBetween(K, xs, ab, cb, dd+1)
                    sndPath = pathBetween(K, xs, cb, cd, dd+1)

                    if fstPath is None:
                        continue
                    if sndPath is None:
                        continue

                    if fstPath[J:-J] != sndPath[J:-J]:
                        continue

                pa = pp[0]
                pc = pp[2]

                cab = xs.get(ab, 0)
                ccb = xs.get(cb, 0)
                ccd = xs.get(cd, 0)

                m = (cab + ccd) / 2.0
                # Assume the true std dev is 10% of the mean
                w = ccb / m

                hgvs = '%s:c.%d_%ddup' % (names[n], pa, pc - 1)
                print '%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%g\t%g' % (hgvs, render(K, ab), xs.get(ab, 0), render(K, cb), xs.get(cb, 0), render(K, cd), xs.get(cd, 0), dd, m, w)

if __name__ == '__main__':
    main(sys.argv[1:])
