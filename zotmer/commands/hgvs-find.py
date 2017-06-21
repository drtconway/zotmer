"""
zot hgvs-find - search for known variants in read data.

Usage:
    zot hgvs-find -X [options] <index> [<variant>...]
    zot hgvs-find [options] <index> <input>...

Options:
    -X              index HGVS variants
    -a              include output for all variants, not just positive results
    -f FILENAME     read variants from a file
    -k K            value of k to use [default: 25]
    -g PATH         directory of FASTQ reference sequences
    -w WIDTH        context width [default: 1000]
    -v              produce verbose output
"""

import array
import math
import sys

import docopt
import yaml

from pykmer.basics import ham, kmer, kmersList, murmer, rc, render
from pykmer.file import openFile, readFasta, readFastq
from pykmer.misc import unionfind
from pykmer.sparse import sparse
from zotmer.library.hgvs import parseHGVS, refSeq2Hg19
from zotmer.library.kmers import kmers
from zotmer.library.files import readKmersAndCounts

def trim(K, kx):
    xs = kx.keys()
    xs.sort()
    xs = sparse(2*K, array.array('L', xs))

    J = K // 2
    K4 = 1 << (2*K)
    i = 0
    uf = unionfind()
    while i < xs.count():
        x = xs.select(i)
        w = x >> (2*J)
        y0 = w << (2*J)
        y1 = (w+1) << (2*J)
        r0 = xs.rank(y0)
        if y1 < K4:
            r1 = xs.rank(y1)
        else:
            r1 = xs.count()
        for j in xrange(r0, r1):
            y0 = xs.select(j)
            for k in xrange(j+1, r1):
                y1 = xs.select(k)
                d = ham(y0, y1)
                if d < 3:
                    uf.union(y0, y1)
                    uf.union(rc(K, y0), rc(K, y1))
        i = r1

    idx = {}
    for i in xrange(xs.count()):
        x = xs.select(i)
        a = uf.find(x)
        if a not in idx:
            idx[a] = []
        idx[a].append(x)

    for xs in idx.values():
        ds = set([])
        for i in xrange(len(xs)):
            x = xs[i]
            xf = float(kx[x])
            for j in xrange(i+1, len(xs)):
                y = xs[j]
                if y in ds:
                    continue
                d = ham(x, y)
                if d >= 3:
                    continue
                yf = float(kx[y])
                v = 1.0/(xf + yf)
                if yf*v < 0.05:
                    #print '%d\t%s\t%d\t%s\t%d' % (d, render(K, x), kx[x], render(K, y), kx[y])
                    ds.add(y)
                elif xf*v < 0.05:
                    #print '%d\t%s\t%d\t%s\t%d' % (d, render(K, y), kx[y], render(K, x), kx[x])
                    ds.add(x)
                    break
        for d in ds:
            del kx[d]

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
    W = int(opts['-w'])

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
