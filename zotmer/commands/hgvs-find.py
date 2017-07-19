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
    -s              merge strands rather than counting them separately
    -t BEDFILE      test an index against a set of genomic regions
    -v              produce verbose output
"""

import array
import math
import sys

import docopt
import yaml

from pykmer.basics import ham, kmer, kmers, kmersList, kmersWithPos, kmersWithPosList, murmer, rc, render
from pykmer.file import openFile, readFasta, readFastq
from pykmer.misc import unionfind
from pykmer.sparse import sparse
from pykmer.stats import logGammaP, logGammaQ
from zotmer.library.hgvs import parseHGVS, refSeq2Hg19
from zotmer.library.files import readKmersAndCounts
from zotmer.library.rope import rope
from zotmer.library.trim import trim

def sqr(x):
    return x*x

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

def kullbackLeibler(Xs, Ys):
    m = 1 + int(0.5 + math.log(max(Xs.keys() + Ys.keys())))

    xs = [0 for i in range(m)]
    for (f, c) in Xs.items():
        j = int(0.5 + math.log(f))
        xs[j] += c
    tx = float(sum(xs))
    px = [x/tx for x in xs]

    ys = [0 for i in range(m)]
    for (f, c) in Ys.items():
        j = int(0.5 + math.log(f))
        ys[j] += c
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
    return (dl, dr)

def chiSquared(Xs, Ys):
    m = 1 + int(0.5 + math.log(max(Xs.keys() + Ys.keys())))

    xs = [0 for i in range(m)]
    for (f, c) in Xs.items():
        j = int(0.5 + math.log(f))
        xs[j] += c
    tx = float(sum(xs))
    px = [x/tx for x in xs]

    ys = [0 for i in range(m)]
    for (f, c) in Ys.items():
        j = int(0.5 + math.log(f))
        ys[j] += c
    ty = float(sum(ys))
    py = [y/ty for y in ys]

    c2 = 0.0
    for (x,y) in zip(px, py):
        assert y > 0
        c2 += sqr(x - y)/y
    return (c2, m)

def chiSquaredPval(x, n, lowerTail=True, logP=False):
    lp = 0
    if lowerTail:
        lp = logGammaP(n/2.0, x/2.0)
    else:
        lp = logGammaQ(n/2.0, x/2.0)
    if logP:
        return lp
    else:
        return math.exp(lp)

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

def applyVariant(sf, v):
    acc = v['accession']
    seq = sf[acc]
    s = rope.atom(seq)
    if v['type'] == 'substitution':
        p = v['position']
        l = rope.substr(s, 0, p)
        m = rope.atom(v['variant'])
        r = rope.substr(s, p + 1, len(s))
        return (s, rope.join([l, m, r]))
    elif v['type'] == 'insertion':
        p = v['after-position'] + 1
        q = v['before-position']
        assert p == q
        l = rope.substr(s, 0, p)
        m = rope.atom(v['sequence'])
        r = rope.substr(s, p, len(s))
        return (s, rope.join([l, m, r]))
    elif v['type'] == 'deletion':
        p = v['first-position']
        q = v['last-position'] + 1
        l = rope.substr(s, 0, p)
        r = rope.substr(s, q, len(s))
        return (s, rope.concat(l, r))
    return None

def context0(k, v, L, sf):
    if v['type'] == 'substitution':
        p = v['position']
        wt = (p - L + 1, p + L)
        mut = wt
    elif v['type'] == 'insertion':
        p = v['after-position'] + 1
        q = v['before-position']
        wt = (p - L + 1, q + L - 1)
        mut = (p - L + 1, q + L - 1 + len(v['sequence']))
    elif v['type'] == 'deletion':
        p = v['first-position']
        q = v['last-position'] + 1
        wt = (p - L + 1, q + L)
        mut = (p - L + 1, p + L)
    else:
        return None

    (r, s) = applyVariant(sf, v)

    wtXs = kmersWithPosList(k, r[wt[0]:wt[1]], False)
    mutXs = kmersWithPosList(k, s[mut[0]:mut[1]], False)

    return (wtXs, mutXs)

def context(k, v, sf):
    if v['type'] == 'substitution':
        p = v['position']
        wt = (p - k + 1, p + k)
        mut = wt
    elif v['type'] == 'insertion':
        p = v['after-position'] + 1
        q = v['before-position']
        wt = (p - k + 1, q + k - 1)
        mut = (p - k + 1, q + k - 1 + len(v['sequence']))
    elif v['type'] == 'deletion':
        p = v['first-position']
        q = v['last-position'] + 1
        wt = (p - k + 1, q + k)
        mut = (p - k + 1, p + k)
    else:
        return None

    (r, s) = applyVariant(sf, v)

    wtXs = set(kmers(k, r[wt[0]:wt[1]], True))
    mutXs = set(kmers(k, s[mut[0]:mut[1]], True))
    com = wtXs & mutXs
    wtXs -= com
    mutXs -= com

    wtXs = list(wtXs)
    wtXs.sort()

    mutXs = list(mutXs)
    mutXs.sort()

    print '%s\t%d\t%d\t%d' % (v['hgvs'], len(wtXs), len(mutXs), len(com))
    return (wtXs, mutXs)

class SequenceFactory(object):
    def __init__(self, home):
        self.home = home
        self.prevAcc = None
        self.prevSeq = None

    def __getitem__(self, acc):
        if acc != self.prevAcc:
            if acc not in refSeq2Hg19:
                print >> sys.stderr, "accession %s not available." % (acc)
            assert acc in refSeq2Hg19
            h = refSeq2Hg19[acc]

            with openFile(self.home + "/" + h + ".fa.gz") as f:
                for (nm,seq) in readFasta(f):
                    self.prevAcc = acc
                    self.prevSeq = seq
                    break
        return self.prevSeq

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = int(opts['-k'])

    d = "."
    if opts['-g']:
        d = opts['-g']
    sf = SequenceFactory(d)

    if opts['-R']:
        v = parseHGVS(opts['<variant>'][0])
        (wt, mut) = context0(K, v, 100, sf)
        idx = {}
        for (x, p) in wt:
            idx[x] = []
        for (x, p) in mut:
            idx[x] = []

        for c in refSeq2Hg19.values():
            print >> sys.stderr, 'scanning %s' % (c,)
            with openFile(d + "/" + c + ".fa.gz") as f:
                for (nm,seq) in readFasta(f):
                    for (y,p) in kmersWithPos(K, seq, True):
                        if y in idx:
                            idx[y].append((c,p))

        for (x,p) in wt:
            for (c,p0) in idx[x]:
                print 'wt\t%s\t%d\t%s\t%d' % (render(K, x), p, c, p0)

        for (x,p) in mut:
            for (c,p0) in idx[x]:
                print 'mut\t%s\t%d\t%s\t%d' % (render(K, x), p, c, p0)

        return

    if opts['-X']:
        variants = opts['<variant>']
        if opts['-f']:
            with openFile(opts['-f']) as f:
                variants += f.read().split()

        wanted = set(['substitution', 'insertion', 'deletion'])
        vx = {}
        for v in variants:
            x = parseHGVS(v)
            if x is None:
                print >> sys.stderr, "unable to parse %s" % (v,)
                continue
            if x['type'] not in wanted:
                print >> sys.stderr, "variant type not suppprted: %s" % (v,)
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
            for v in vs:
                (wtYs, mutYs) = context(K, v, sf)
                if len(wtYs) == 0 or len(mutYs) == 0:
                    print >> sys.stderr, "variant has no distinguishing k-mers: %s" % (v['hgvs'],)
                    continue
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

    combineStrands = False
    if opts['-s']:
        combineStrands = True

    kx = {}
    rn = 0
    M = (1 << 18) - 1
    for (fn1, fn2) in pairs(opts['<input>']):
        with openFile(fn1) as f1, openFile(fn2) as f2:
            for fq1, fq2 in both(readFastq(f1), readFastq(f2)):
                rn += 1
                if (rn & M) == 0 and opts['-v']:
                    print >> sys.stderr, 'read pairs processed: %d' % (rn,)
                xs = kmersList(K, fq1[1], combineStrands) + kmersList(K, fq2[1], combineStrands)
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

    #trim(K, kx)

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
            xs = [max(1, x) for x in posRes[h].values()]
            nx = float(len(xs))
            xh = hist(xs)
            ys = [max(1, y) for y in negRes[h].values()]
            ny = float(len(ys))
            yh = hist(ys)

            (posD, posDF) = chiSquared(xh, zh)
            posP = chiSquaredPval(posD, posDF, lowerTail=False)
            (negD, negDF) = chiSquared(yh, zh)
            negP = chiSquaredPval(negD, negDF, lowerTail=False)

            if opts['-v']:
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
                    print '%s\t%d\t%s\t%d\t%g\t%g\t%g' % (h, c, 'Glo', fs[0], tx[0]/nz, posD, negD)
                    print '%s\t%d\t%s\t%d\t%g\t%g\t%g' % (h, c, 'Pos', fs[1], tx[1]/nx, posD, negD)
                    print '%s\t%d\t%s\t%d\t%g\t%g\t%g' % (h, c, 'Neg', fs[2], tx[2]/ny, posD, negD)
            elif opts['-a']:
                print '%s\t%g\t%g\t%g\t%g' % (h, posD, negD, posP, negP)
            sys.stdout.flush()

    if False:
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
