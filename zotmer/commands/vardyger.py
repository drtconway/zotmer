
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

from pykmer.basics import ham, kmersList, kmersWithPosList, kmersWithPosLists, rc, render
from pykmer.file import openFile, readFasta, readFastq
from zotmer.library.debruijn import interpolate, interpolateBetween, pathBetween
from zotmer.library.hgvs import Duplication
from zotmer.library.reads import reads

def pairs(xs):
    N = len(xs)
    for i in xrange(N):
        for j in xrange(i+1, N):
            yield (xs[j], xs[i])
            yield (xs[i], xs[j])

def findDupDeNovo(xFst, xLst):
    for (c, ys) in xFst.iteritems():
        for (b, d) in pairs(ys):
            for a in xLst[b]:
                if a == c:
                    continue
                yield (a, b, c, d)

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

def nearest(m, xs, y):
    dMin = m
    posMin = []
    for (x,p) in xs.iteritems():
        d = ham(x, y)
        if d > dMin:
            continue
        if d < dMin:
            dMin = d
            posMin = []
        posMin += p
    return posMin

def offset(xs, ys, d):
    for x in xs:
        for y in ys:
            if x + d == y:
                yield (x,y)

def positions(posIdx, J, a, b, c, d):
    D = 2
    pas = nearest(D, posIdx, a)
    pbs = nearest(D, posIdx, b)
    pabs = list(offset(pas, pbs, J))
    pcs = nearest(D, posIdx, c)
    pds = nearest(D, posIdx, d)
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

def locate(K, idx, seq):
    hits = {}
    (xs, ys) = kmersWithPosLists(K, seq)
    for (zs, s) in [(xs, True), (ys, False)]:
        for (x,p) in zs:
            if x not in idx:
                continue
            p -= 1
            for q in idx[x]:
                qp = (q - p, s)
                if qp not in hits:
                    hits[qp] = 0
                hits[qp] += 1

    cand = None
    for ((p,s),c) in hits.items():
        if cand is None:
            cand = (c, [(p, s)])
        elif c > cand[0]:
            cand = (c, [(p, s)])
        elif c == cand[0]:
            cand[1].append((p, s))

    if cand is None:
        return []

    res = []
    (c, pss) = cand
    for (p0, s) in pss:
        if s:
            zs = [(x, p0 + p - 1) for (x,p) in xs]
        else:
            zs = [(y, p0 + p - 1) for (y,p) in ys]
        res += zs
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
    seqs = {}
    bait = {}
    wtFst = []
    wtLst = []
    posIdx = []
    rds = []
    with openFile(opts['<sequences>']) as f:
        for (nm, seq) in readFasta(f):
            n = len(names)
            names.append(nm)
            seqs[nm] = seq
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

            rds.append([])

    N = len(names)

    L = None
    X = [{} for n in range(N)]
    for itm in reads(opts['<input>'], K=K, reads=True, kmers=True, both=True, verbose=verbose):
        rd = itm[0][0]
        L = len(rd)

        xs = itm[1][0][0]
        hits = set([])
        for x in xs:
            if x in bait:
                hits |= bait[x]
        for n in hits:
            for x in xs:
                if x not in X[n]:
                    X[n][x] = 0
                X[n][x] += 1
            rds[n].append(rd)

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

        for (a, b, c, d) in findDupDeNovo(fst, lst):
            #print [render(J, w) for w in [a, b, c, d]]
            #continue
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
                pb = pp[1]
                pc = pp[2]
                pd = pp[3]

                cab = xs.get(ab, 0)
                ccb = xs.get(cb, 0)
                ccd = xs.get(cd, 0)

                m = (cab + ccd) / 2.0
                # Assume the true std dev is 10% of the mean
                w = ccb / m

                hgvs = '%s:c.%d_%ddup' % (names[n], pb, pd - 1)

                v = Duplication(names[n], pb, pd-1)
                v.setSequenceFactory(seqs)

                ctx = v.context(2*L)

                idx = {}
                for (x,p) in kmersWithPosList(K, ctx[1], False):
                    if x not in idx:
                        idx[x] = []
                    idx[x].append(p)
                res = {}
                for fq in rds[n]:
                    for (x,p) in locate(K, idx, fq[1]):
                        if p not in res:
                            res[p] = {}
                        if x not in res[p]:
                            res[p][x] = 0
                        res[p][x] += 1

                for (p,ys) in sorted(res.items()):
                    for (y,c) in sorted(ys.items()):
                        print '%d\t%s\t%d' % (p, render(K, y), c)
                print '%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%g\t%g' % (hgvs, render(K, ab), xs.get(ab, 0), render(K, cb), xs.get(cb, 0), render(K, cd), xs.get(cd, 0), dd, m, w)

if __name__ == '__main__':
    main(sys.argv[1:])
