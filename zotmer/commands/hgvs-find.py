"""
zot hgvs-find - search for known variants in read data.

Usage:
    zot hgvs-find -R [options] <variant>
    zot hgvs-find -X [options] <index> [<variant>...]
    zot hgvs-find [options] <index> <input>...

Options:
    -R              report the ambiguity associated with a variant
    -X              index HGVS variants
    -a              include output for all variants, not just positive results
    -f FILENAME     read variants from a file
    -k K            value of k to use [default: 25]
    -g PATH         directory of FASTQ reference sequences
    -t BEDFILE      test an index against a set of genomic regions
    -v              produce verbose output
"""

import array
import math
import sys

import docopt
import yaml

from pykmer.basics import ham, kmer, kmers, kmersList, kmersWithPos, murmer, rc, render
from pykmer.file import openFile, readFasta, readFastq
from pykmer.misc import unionfind
from pykmer.sparse import sparse
from zotmer.library.hgvs import parseHGVS, refSeq2Hg19
from zotmer.library.kmers import kmers
from zotmer.library.files import readKmersAndCounts
from zotmer.library.rope import rope
from zotmer.library.trim import trim

def merge(rs):
    r = None
    for r0 in rs:
        if r is None:
            r = r0
            continue
        if r[1] < r0[0]:
            yield r
            r = r0
        else:
            r = (min(r[0], r0[0]), max(r[1], r0[1]))
    if r is not None:
        yield r

def lcp(xs):
    z = xs[0]
    for x in xs[1:]:
        i = 0
        while i < len(z) and i < len(x) and x[i] == z[i]:
            i += 1
        z = z[:i]
    return z

def mkRnks(xs):
    zs = {}
    for x in xs:
        if x not in zs:
            zs[x] = 0
        zs[x] += 1
    zs = zs.items()
    zs.sort()

    rs = {}
    n = 1
    for (z, c) in zs:
        r = n + (c - 1.0) / 2.0
        rs[z] = r
        n += c
    return rs

def mannWhitney(xs, ys):
    rs = mkRnks(xs + ys)

    nx = len(xs)
    ny = len(ys)

    rx = sum([rs[x] for x in xs])
    #ry = sum([rs[y] for y in ys])

    ux = rx - nx * (ny + 1) / 2.0
    #uy = ry - ny * (nx + 1) / 2.0

    mu = nx * ny / 2.0
    sd = math.sqrt(nx*ny*(nx + ny + 1)/12.0)

    zx = (ux - mu) / sd
    #zy = (uy - mu) / sd

    px = 0.5*math.erfc(zx/math.sqrt(2))
    #py = 0.5*math.erfc(zy/math.sqrt(2))

    return (ux, zx, px)

def kullbackLeibler(xs, ys):
    assert len(xs) == len(ys)

    xs = [max(1, x) for x in xs]
    tx = float(sum(xs))
    px = [x/tx for x in xs]

    ys = [max(1, y) for y in ys]
    ty = float(sum(ys))
    py = [y/ty for y in ys]

    d = 0.0
    for (x,y) in zip(px, py):
        if x == 0.0:
            continue
        assert y > 0
        d += x * math.log(x/y)
    return d

def hist(xs):
    r = {}
    for x in xs:
        if x not in r:
            r[x] = 0
        r[x] += 1
    return r

def kolmogorovSmirnov(hx, hy):
    zs = list(set(hx.keys() + hy.keys()))
    zs.sort()

    cx = 0
    tx = float(sum(hx.values()))
    cy = 0
    ty = float(sum(hy.values()))

    dl = 0.0
    dr = 0.0
    for z in zs:
        if z in hx:
            cx += hx[z]
        if z in hy:
            cy += hy[z]
        px = cx / tx
        py = cy / ty
        dl = max(dl, px - py)
        dr = max(dr, py - px)
        #print '%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g' % (z, cx, cy, px, py, px - py, py - px, dl, dr)
    return (dl, dr)

def pairs(xs):
    assert len(xs) & 1 == 0
    i = 0
    while i + 1 < len(xs):
        yield (xs[i], xs[i + 1])
        i += 2

def both(xs, ys):
    while True:
        try:
            x = xs.next()
            y = ys.next()
            yield(x, y)
        except StopIteration:
            return

def neigh(K, x, d):
    if d == 0:
        return

    for j in xrange(K):
        for y in [1,2,3]:
            z = x ^ (y << (2*j))
            if d == 1:
                yield z
            else:
                for w in neigh(K, z, d - 1):
                    if ham(x, w) == d:
                        yield w


def ball(K, xs, d):
    ys = set(xs)
    i = 0
    while i < d:
        i += 1
        for x in xs:
            ys |= set(neigh(K, x, i))
    return ys

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = int(opts['-k'])

    if opts['-R']:
        x = parseHGVS(opts['<variant>'][0])

        d = "."
        if opts['-g']:
            d = opts['-g']

        acc = x['accession']
        if acc not in refSeq2Hg19:
            print >> sys.stderr, "accession %s not supported" % (acc,)
            return
        h = refSeq2Hg19[acc]

        with openFile(d + "/" + h + ".fa.gz") as f:
            for (nm,seq) in readFasta(f):
                p = x['position'] - 1
                wt = seq[p-K:p + K]
                mut = wt[:K] + x['variant'] + wt[K+1:]
                wtXs = set(kmersList(K, wt, True))
                mutXs = set(kmersList(K, mut, True))
                wtYs = list(wtXs - mutXs)
                wtYs.sort()
                mutYs = list(mutXs - wtXs)
                mutYs.sort()

        wt = dict([(y,[]) for y in wtYs])
        mut = dict([(y,[]) for y in mutYs])

        for c in refSeq2Hg19.values():
            print >> sys.stderr, 'scanning %s' % (c,)
            with openFile(d + "/" + c + ".fa.gz") as f:
                for (nm,seq) in readFasta(f):
                    for (y,p) in kmersWithPos(K, seq, True):
                        if y in wt:
                            wt[y].append((c,p))
                        if y in mut:
                            mut[y].append((c,p))

        for (y,ps) in wt.iteritems():
            for p in ps:
                print 'wt\t%s\t%s\t%d' % (render(K, y), p[0], p[1])

        for (y,ps) in mut.iteritems():
            for p in ps:
                print 'mut\t%s\t%s\t%d' % (render(K, y), p[0], p[1])

        return

    if opts['-X']:
        variants = opts['<variant>']
        if opts['-f']:
            with openFile(opts['-f']) as f:
                variants += f.read().split()
        vx = {}
        for v in variants:
            x = parseHGVS(v)
            if x is None:
                print >> sys.stderr, "unable to parse %s" % (v,)
                continue
            if x['type'] != 'substitution':
                print >> sys.stderr, "only substitutions are supported at this time: %s" % (v,)
                continue
            x['hgvs'] = v
            acc = x['accession']
            if acc not in vx:
                vx[acc] = []
            vx[acc].append(x)

        d = "."
        if opts['-g']:
            d = opts['-g']

        rs = []
        for (acc, vs) in vx.iteritems():
            if acc not in refSeq2Hg19:
                print >> sys.stderr, "accession %s not supported" % (acc,)
                continue
            h = refSeq2Hg19[acc]
            with openFile(d + "/" + h + ".fa.gz") as f:
                for (nm,seq) in readFasta(f):
                    for v in vs:
                        p = v['position'] - 1
                        wt = seq[p-K:p + K]
                        mut = wt[:K] + v['variant'] + wt[K+1:]
                        wtXs = set(kmersList(K, wt, True))
                        mutXs = set(kmersList(K, mut, True))
                        wtYs = list(wtXs - mutXs)
                        wtYs.sort()
                        mutYs = list(mutXs - wtXs)
                        mutYs.sort()

                        r = {}
                        r['hgvs'] = v['hgvs']
                        r['pos'] = [render(K, y) for y in mutYs]
                        r['neg'] = [render(K, y) for y in wtYs]
                        rs.append(r)
        with open(opts['<index>'], 'w') as f:
            yaml.safe_dump(rs, f)

        return

    with open(opts['<index>']) as f:
        itms = yaml.load(f)

    posIdx = {}
    posRes = {}
    negIdx = {}
    negRes = {}
    hs = set([])
    univ = set([])
    for itm in itms:
        h = itm['hgvs']
        hs.add(h)
        posRes[h] = {}
        negRes[h] = {}
        for x in itm['pos']:
            y = kmer(x)
            univ.add(y)
            if y not in posIdx:
                posIdx[y] = []
            posIdx[y].append(h)
            posRes[h][y] = 0
        for x in itm['neg']:
            y = kmer(x)
            univ.add(y)
            if y not in negIdx:
                negIdx[y] = []
            negIdx[y].append(h)
            negRes[h][y] = 0

    if opts['-t']:
        d = "."
        if opts['-g']:
            d = opts['-g']

        roi = {}
        with openFile(opts['-t']) as f:
            for l in f:
                if l[0] == '#':
                    continue
                t = l.split('\t')
                if t[0] == 'track':
                    continue
                c = t[0]
                s = int(t[1])
                e = int(t[2])
                if c not in roi:
                    roi[c] = set([])
                roi[c].add((s, e))

        kx = {}
        for x in univ:
            kx[x] = 0

        roiCs = roi.keys()
        roiCs.sort()
        for c in roiCs:
            print >> sys.stderr, 'scanning %s' % (c,)
            rs = list(roi[c])
            rs.sort()
            with openFile(d + "/" + c + ".fa.gz") as f:
                for (nm,seq) in readFasta(f):
                    w = rope.atom(seq)
                    for (s,e) in merge(rs):
                        v = rope.substr(w, s, e)
                        for x in kmersList(K, v[:], True):
                            if x in kx:
                                kx[x] += 1

        for h in hs:
            xs = posRes[h].keys()
            posHist = {}
            for x in xs:
                c = kx.get(x, 0)
                if c not in posHist:
                    posHist[c] = 0
                posHist[c] += 1
            ys = negRes[h].keys()
            negHist = {}
            for y in ys:
                c = kx.get(y, 0)
                if c not in negHist:
                    negHist[c] = 0
                negHist[c] += 1
            posU = posHist.get(1, 0)
            posT = float(sum(posHist.values()))
            posT = max(1.0, posT)
            posAvg = 0.0
            for (c,f) in posHist.items():
                posAvg += c*f
            posAvg /= float(sum(posHist.values()))
            negU = negHist.get(1, 0)
            negT = float(sum(negHist.values()))
            negT = max(1.0, negT)
            negAvg = 0.0
            for (c,f) in negHist.items():
                negAvg += c*f
            negAvg /= float(sum(negHist.values()))
            if opts['-v']:
                hh = set(posHist.keys() + negHist.keys())
                hh = list(hh)
                hh.sort()
                for k in hh:
                    pc = posHist.get(k, 0)
                    nc = negHist.get(k, 0)
                    print '%s\t%g\t%g\t%g\t%g\t%d\t%d\t%d' % (h, posU/posT, posAvg, negU/negT, negAvg, k, pc, nc)
            else:
                print '%s\t%g\t%g\t%g\t%g' % (h, posU/posT, posAvg, negU/negT, negAvg)

        return

    kx = {}
    for (fn1, fn2) in pairs(opts['<input>']):
        with openFile(fn1) as f1, openFile(fn2) as f2:
            for fq1, fq2 in both(readFastq(f1), readFastq(f2)):
                xs = kmersList(K, fq1[1]) + kmersList(K, fq2[1])
                found = False
                for x in xs:
                    if x in univ:
                        found = True
                        break
                if not found:
                    continue
                for x in xs:
                    if x not in kx:
                        kx[x] = 0
                    kx[x] += 1
    trim(K, kx)

    zh = {}
    for (x, c) in kx.iteritems():
        if c not in zh:
            zh[c] = 0
        zh[c] += 1
        if x in posIdx:
            for h in posIdx[x]:
                posRes[h][x] = c
        if x in negIdx:
            for h in negIdx[x]:
                negRes[h][x] = c
    nz = float(len(kx))

    if True:
        for h in hs:
            xs = posRes[h].values()
            nx = float(len(xs))
            xh = hist(xs)
            ys = negRes[h].values()
            ny = float(len(ys))
            yh = hist(ys)

            alpha = 0.05

            (_posL, posR) = kolmogorovSmirnov(xh, zh)
            posW = math.sqrt((nx + nz)/(nx*nz))
            posC = math.sqrt(-0.5*math.log(alpha/2.0)) * posW
            posA = posR / posC

            (_negL, negR) = kolmogorovSmirnov(yh, zh)
            negW = math.sqrt((ny + nz)/(ny*nz))
            negC = math.sqrt(-0.5*math.log(alpha/2.0)) * negW
            negA = negR / negC

            if opts['-v'] and (opts['-a'] or posA > 1.0 or negA > 1.0):
                vx = {}
                for (c,f) in zh.iteritems():
                    if c not in vx:
                        vx[c] = [0, 0, 0]
                    vx[c][0] = f
                for (c,f) in xh.iteritems():
                    if c not in vx:
                        vx[c] = [0, 0, 0]
                    vx[c][1] = f
                for (c,f) in yh.iteritems():
                    if c not in vx:
                        vx[c] = [0, 0, 0]
                    vx[c][2] = f
                vx = vx.items()
                vx.sort()
                tx = [0, 0, 0]
                for i in xrange(len(vx)):
                    (c, fs) = vx[i]
                    tx[0] += fs[0]
                    tx[1] += fs[1]
                    tx[2] += fs[2]
                    print '%s\t%d\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g' % (h, c, 'Glo', fs[0], tx[0]/nz, posR, negR, posC, negC, posA, negA)
                    print '%s\t%d\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g' % (h, c, 'Pos', fs[1], tx[1]/nx, posR, negR, posC, negC, posA, negA)
                    print '%s\t%d\t%s\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g' % (h, c, 'Neg', fs[2], tx[2]/ny, posR, negR, posC, negC, posA, negA)
            elif opts['-a']:
                print '%s\t%g\t%g\t%g\t%g\t%g\t%g' % (h, posR, negR, posC, negC, posA, negA)
            elif posA > 1.0:
                print '%s\t%g\t%g\t%g\t%g\t%g\t%g' % (h, posR, negR, posC, negC, posA, negA)
            sys.stdout.flush()

if __name__ == '__main__':
    main(sys.argv[1:])
