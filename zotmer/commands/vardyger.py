
"""
zot vardyger - search for tandem duplications

Usage:
    zot vardyger [options] <sequences> <input>...

Options:
    -A              show anchors in reads
    -a              report all detected duplications, not just in-frame ones.
    -k K            value of K to use - must be even. [default: 24]
    -m MIN          minimum k-mer coverage for reporting anchors [default: 10]
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
from zotmer.library.debruijn import interpolate
from zotmer.library.hgvs import Duplication
from zotmer.library.reads import reads
from zotmer.library.align import revComp

def pairs(xs, ys):
    N = len(xs)
    M = len(ys)
    for i in xrange(N):
        for j in xrange(M):
            if xs[i] != ys[j]:
                yield (xs[i], ys[j])

def findDupDeNovo(xFst, xLst):
    for (c, ys) in xFst.iteritems():
        for (b, d) in pairs(ys, ys):
            for a in xLst[b]:
                if a == c:
                    continue
                yield (a, b, c, d)

def findDup(refFst, refLst, xFst, xLst):
    cs = set(xFst.keys()) & set(refFst.keys())
    Ws = set(refLst.keys())
    for c in cs:
        ys = xFst[c]
        Ys = set(ys)
        ds = sorted(Ys & set(refFst[c]))
        bs = sorted(Ys & Ws)
        for (b, d) in pairs(bs, ds):
            for a in refLst[b]:
                if a == c:
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
    D = 1
    pas = nearest(D, posIdx, a)
    pbs = nearest(D, posIdx, b)
    pabs = list(offset(pas, pbs, J))
    pcs = nearest(D, posIdx, c)
    pds = nearest(D, posIdx, d)
    pcds = list(offset(pcs, pds, J))

    res = None
    for (pa, pb) in pabs:
        for (pc, pd) in pcds:
            #print (pa, pb, pc, pd)
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

def remapReads(K, L, rds, v):
    ctx = v.context(2*L)
    idx = {}
    for (x,p) in kmersWithPosList(K, ctx[1], False):
        if x not in idx:
            idx[x] = []
        idx[x].append(p)

    res = {}
    for fq in rds:
        for (x,p) in locate(K, idx, fq[1]):
            if p not in res:
                res[p] = {}
            if x not in res[p]:
                res[p][x] = 0
            res[p][x] += 1

    for (p,ys) in sorted(res.items()):
        for (y,c) in sorted(ys.items()):
            print '%d\t%s\t%d' % (p, render(K, y), c)

def showAnchoredReads(K, anks, rds):
    A = set(anks.keys())

    res = {}
    for rd in rds:
        (xs, ys) = kmersWithPosLists(K, rd[1])

        fwd = set([])
        for (x,p) in xs:
            if x in A:
                fwd.add((x,p))

        rev = set([])
        for (y,p) in ys:
            if y in A:
                rev.add((y,p))

        if len(fwd) > 0:
            p0 = min([p for (x,p) in fwd])
            fwd = frozenset([(x,p-p0) for (x,p) in fwd])
            seq = rd[1]
            if fwd not in res:
                res[fwd] = {}
            k = (p0,seq)
            if k not in res[fwd]:
                res[fwd][k] = 0
            res[fwd][k] += 1

        if len(rev) > 0:
            p0 = min([p for (y,p) in rev])
            rev = frozenset([(y,p-p0) for (y,p) in rev])
            seq = revComp(rd[1])
            if rev not in res:
                res[rev] = {}
            k = (p0,seq)
            if k not in res[rev]:
                res[rev][k] = 0
            res[rev][k] += 1

    for s in sorted(res.keys()):
        q = max([p for (p,seq) in res[s].keys()])
        lab = ','.join(sorted([anks[x] for (x,p) in s]))
        hdr = [lab,'\t\t', (q-1)*' ']
        pp = None
        for (p,x) in sorted([(p1,x1) for (x1,p1) in s]):
            if pp is None:
                hdr += [render(K, x)]
            else:
                d = p - pp - K
                hdr += [d*' ', render(K, x)]
            pp = p
        print ''.join(hdr)
        vv = []
        for (k,c) in sorted(res[s].items()):
            (p0, seq) = k
            vv.append((p0, c, seq))
        for (p0, c, seq) in sorted(vv):
            print '%d\t%d\t%s%s' % (p0, c, (q-p0)*' ', seq)
        print

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

    minCov = int(opts['-m'])

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

            for (a, b, c, d) in findDup(wtFst[n], wtLst[n], wtFst[n], wtLst[n]):
                pps = positions(posIdx[n], J, a, b, c, d)
                if pps is None:
                    continue
                for pp in pps:
                    ab = a << S | b
                    cb = c << S | b
                    cd = c << S | d
                    dd = pp[2] - pp[0]
                    print >> sys.stderr, 'warning: phantom dumplication: %s-%s-%s (%d)' % (render(K, ab), render(K, cb), render(K, cd), dd)

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

    hdrShown = False
    vn = 0
    for n in range(N):
        xs = X[n]

        fst = {}
        lst = {}
        for (x,c) in xs.iteritems():
            #if c < 5:
            #    continue
            y0 = x >> S
            y1 = x & Mj

            if y0 not in fst:
                fst[y0] = []
            fst[y0].append(y1)

            if y1 not in lst:
                lst[y1] = []
            lst[y1].append(y0)

        #for (a, b, c, d) in findDupDeNovo(fst, lst):
        for (a, b, c, d) in findDup(wtFst[n], wtLst[n], fst, lst):
            #continue
            pps = positions(posIdx[n], J, a, b, c, d)
            if pps is None:
                continue
            for pp in pps:
                ab = a << S | b
                cb = c << S | b
                cd = c << S | d
                #print [(render(J, w), p) for (w,p) in zip([a, b, c, d], pps)]

                dd = pp[2] - pp[0]

                if not opts['-a'] and dd % 3 != 0:
                    continue

                if opts['-s']:
                    fstPath = interpolate(K, xs, ab, cb, dd+1)
                    sndPath = interpolate(K, xs, cb, cd, dd+1)

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

                if cab < minCov:
                    continue
                if ccb < minCov:
                    continue
                if ccd < minCov:
                    continue

                m = (cab + ccd) / 2.0
                # Assume the true std dev is 10% of the mean
                w = ccb / m

                hgvs = '%s:c.%d_%ddup' % (names[n], pb, pd - 1)
                v = Duplication(names[n], pb, pd-1, seqs)
                if opts['-A']:
                    showAnchoredReads(K, {ab:'AB', cb:'CB', cd:'CD'}, rds[n])

                vn += 1

                hdrs = ['n']
                fmts = ['%d']
                outs = [vn]

                hdrs += ['left', 'leftCov']
                fmts += ['%s','%d']
                outs += [render(K, ab), cab]

                hdrs += ['mid', 'midCov']
                fmts += ['%s','%d']
                outs += [render(K, cb), ccb]

                hdrs += ['right', 'rightCov']
                fmts += ['%s','%d']
                outs += [render(K, cd), ccd]

                hdrs += ['len']
                fmts += ['%d']
                outs += [dd]

                hdrs += ['vaf']
                fmts += ['%g']
                outs += [w]

                hdrs += ['hgvs']
                fmts += ['%s']
                outs += [hgvs]

                if not hdrShown:
                    hdrShown = True
                    print '\t'.join(hdrs)
                print '\t'.join(fmts) % tuple(outs)

if __name__ == '__main__':
    main(sys.argv[1:])
