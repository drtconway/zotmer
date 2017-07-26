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
    -e              estimate cutoffs for KL statistic
    -f FILENAME     read variants from a file
    -k K            value of k to use [default: 25]
    -l              print long form results
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
from pykmer.stats import logGammaP, logGammaQ, logLowerGamma
from zotmer.library.hgvs import applyVariant, parseHGVS, refSeq2Hg19
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

def hist(xs):
    r = {}
    for x in xs:
        if x not in r:
            r[x] = 0
        r[x] += 1
    return r

def quantiles(Xs, n):
    N = sum(Xs.values())
    cp = {}
    t = 0
    for x in sorted(Xs.keys()):
        c = Xs[x]
        t += Xs[x]
        cp[x] = float(t - 1)/float(N)
    r = [0 for i in range(n)]
    for (x,p) in cp.items():
        q = int(p*n)
        r[q] = max(r[q], x)
    x0 = 0
    for q in range(n):
        r[i] = max(r[i], x0)
        x0 = r[i]
    return r

def binLog(Xs, n, w = None):
    if w is not None:
        (lx, ux) = w
    else:
        lx = min(Xs.keys())
        ux = max(Xs.keys())
        assert lx > 0

    ly = math.log(1 + lx)
    uy = math.log(1 + ux + 1)
    ry = uy - ly
    r = {}
    for (x, c) in Xs.items():
        y = int(n * (math.log(1 + x) - ly) / ry)
        if y not in r:
            r[y] = 0
        r[y] += c
    return r

def binSqrt(Xs, n, w = None):
    if w is not None:
        (lx, ux) = w
    else:
        lx = min(Xs.keys())
        ux = max(Xs.keys())
        assert lx >= 0

    ly = math.sqrt(lx)
    uy = math.sqrt(ux + 1)
    ry = uy - ly
    r = {}
    for (x, c) in Xs.items():
        y = int(n * (math.sqrt(x) - ly) / ry)
        if y not in r:
            r[y] = 0
        r[y] += c
    return r

def binHist(Xs, tx = None):
    if tx is None:
        tx = quantiles(Xs, 10)
    r = []
    j = 0
    n = 0
    for x in sorted(Xs.keys()):
        while x > tx[j]:
            r.append(n)
            n = 0
            j += 1
        n += Xs[x]
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

def kullbackLeibler(Xs, Ys):
    bs = Ys.keys()
    bs.sort()

    zx = float(sum(Xs.values()))
    px = {}
    for b in bs:
        c = Xs.get(b, 0)
        px[b] = c/zx

    zy = float(sum(Ys.values()))
    py = {}
    for b in bs:
        c = Ys.get(b, 0)
        py[b] = c/zy

    d = 0.0
    for b in bs:
        x = px[b]
        y = py[b]
        if x == 0.0:
            continue
        assert y > 0
        d += x * math.log(x/y)
    return d

def estimateGamma(xs):
    N = float(len(xs))
    s = math.log(sum(xs)/N) - sum([math.log(x) for x in xs])/N
    k = (3 - s + math.sqrt(sqr(s - 3) + 24*s)) / (12*s)
    t = sum(xs) / (k*N)
    return (k, t)

def estimateGammaCutoff(xs, p):
    (k, t) = estimateGamma(xs)
    print >> sys.stderr, 'estimated gamma parameters: %g\t%g' % (k, t)

    gk = math.lgamma(k)

    p = math.log(p)

    xl = 1e-100
    xh = max(xs) + 1
    for i in xrange(20):
        x = (xh + xl) / 2.0
        q = logLowerGamma(k, x/t) - gk
        #print '(%g, %g)\t%g\t%g\t%g' % (xl, xh, x, q, p)
        if q > p:
            xh = x
        elif q < p:
            xl = x
        else:
            break
    return x

def chiSquared(Xs, Ys):
    bs = Ys.keys()
    bs.sort()

    zx = float(sum(Xs.values()))
    px = {}
    for b in bs:
        c = Xs.get(b, 0)
        px[b] = c/zx

    zy = float(sum(Ys.values()))
    py = {}
    for b in bs:
        c = Ys.get(b, 0)
        py[b] = c/zy

    c2 = 0.0
    for b in bs:
        x = px[b]
        y = py[b]
        #print '%d\t%g\t%g\t%g' % (b, x, y, sqr(x - y)/y)
        c2 += sqr(x - y)/y
    return (c2, len(bs))

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

def context0(k, v, L, sf):
    if v['type'] == 'substitution':
        p = v['position'] - 1
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

    acc = v['accession']
    seq = sf[acc]
    s = rope.atom(seq)

    r = applyVariant(s, v)

    #print (' '*(L-1)) + '.'
    #print r[wt[0]:wt[1]]
    #print s[mut[0]:mut[1]]

    wtXs = kmersWithPosList(k, r[wt[0]:wt[1]], False)
    mutXs = kmersWithPosList(k, s[mut[0]:mut[1]], False)

    return (wtXs, mutXs)

def context(k, v, sf):
    if v['type'] == 'substitution':
        p = v['position'] - 1
        wt = (p - k + 1, p + k)
        mut = wt
    elif v['type'] == 'insertion':
        p = v['after-position']
        q = v['before-position'] - 1
        assert p == q
        wt = (p - k + 1, q + k - 1)
        mut = (p - k + 1, q + k - 1 + len(v['sequence']))
    elif v['type'] == 'deletion':
        p = v['first-position'] - 1
        q = v['last-position']
        wt = (p - k + 1, q + k)
        mut = (p - k + 1, p + k)
    else:
        return None

    acc = v['accession']
    seq = sf[acc]
    s = rope.atom(seq)

    r = applyVariant(s, v)

    wtXs = set(kmers(k, s[wt[0]:wt[1]], True))
    mutXs = set(kmers(k, r[mut[0]:mut[1]], True))
    com = wtXs & mutXs
    wtXs -= com
    mutXs -= com

    wtXs = list(wtXs)
    wtXs.sort()

    mutXs = list(mutXs)
    mutXs.sort()

    print '%s\t%s\t%d\t%d\t%d' % (v['hgvs'], v['type'], len(wtXs), len(mutXs), len(com))
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
                print >> sys.stderr, "variant type not supported: %s" % (v,)
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
                r['wt'] = [render(K, y) for y in wtYs]
                r['mut'] = [render(K, y) for y in mutYs]
                rs.append(r)

        with open(opts['<index>'], 'w') as f:
            yaml.safe_dump(rs, f)

        return

    with open(opts['<index>']) as f:
        itms = yaml.load(f)

    wtIdx = {}
    wtRes = {}
    mutIdx = {}
    mutRes = {}
    hs = set([])
    univ = set([])
    for itm in itms:
        h = itm['hgvs']
        hs.add(h)
        wtRes[h] = {}
        for x in itm['wt']:
            y = kmer(x)
            univ.add(y)
            if y not in wtIdx:
                wtIdx[y] = []
            wtIdx[y].append(h)
            wtRes[h][y] = 0
        mutRes[h] = {}
        for x in itm['mut']:
            y = kmer(x)
            univ.add(y)
            if y not in mutIdx:
                mutIdx[y] = []
            mutIdx[y].append(h)
            mutRes[h][y] = 0

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
            xs = wtRes[h].keys()
            wtHist = {}
            for x in xs:
                c = kx.get(x, 0)
                if c not in wtHist:
                    wtHist[c] = 0
                wtHist[c] += 1
            wtU = wtHist.get(1, 0)
            wtT = float(sum(wtHist.values()))
            wtT = max(1.0, wtT)
            wtAvg = 0.0
            ys = mutRes[h].keys()
            mutHist = {}
            for y in ys:
                c = kx.get(y, 0)
                if c not in mutHist:
                    mutHist[c] = 0
                mutHist[c] += 1
            for (c,f) in wtHist.items():
                wtAvg += c*f
            wtAvg /= float(sum(wtHist.values()))
            mutU = mutHist.get(1, 0)
            mutT = float(sum(mutHist.values()))
            mutT = max(1.0, mutT)
            mutAvg = 0.0
            for (c,f) in mutHist.items():
                mutAvg += c*f
            mutAvg /= float(sum(mutHist.values()))
            if opts['-v']:
                hh = set(mutHist.keys() + wtHist.keys())
                hh = list(hh)
                hh.sort()
                for k in hh:
                    pc = mutHist.get(k, 0)
                    nc = wtHist.get(k, 0)
                    print '%s\t%g\t%g\t%g\t%g\t%d\t%d\t%d' % (h, mutU/mutT, mutAvg, wtU/wtT, wtAvg, k, pc, nc)
            else:
                print '%s\t%g\t%g\t%g\t%g' % (h, mutU/mutT, mutAvg, wtU/wtT, wtAvg)

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

    B = 10

    zh0 = {}
    for (x, c) in kx.iteritems():
        if c not in zh0:
            zh0[c] = 0
        zh0[c] += 1
        if x in wtIdx:
            for h in wtIdx[x]:
                wtRes[h][x] = c
        if x in mutIdx:
            for h in mutIdx[x]:
                mutRes[h][x] = c
    nz = float(len(kx))
    zz = (0, max(zh0.keys()))
    zh = binLog(zh0, B, zz)
    zp = [zh[b]/float(nz) for b in sorted(zh.keys())]

    hdrShown = False

    kls = []

    for h in hs:
        hdrs = ['hgvs']
        fmts = ['%s']
        outs = [h]

        xs = [x for x in wtRes[h].values()]
        nx = float(len(xs))
        xh0 = hist(xs)
        xh = binLog(xh0, B, zz)
        (wtD, wtDF) = chiSquared(xh, zh)
        wtP = chiSquaredPval(wtD, wtDF, lowerTail=False)

        ys = [y for y in mutRes[h].values()]
        ny = float(len(ys))
        yh0 = hist(ys)
        yh = binLog(yh0, B, zz)
        (mutD, mutDF) = chiSquared(yh, zh)
        mutP = chiSquaredPval(mutD, mutDF, lowerTail=False)

        if True:
            hdrs += ['wtChi2', 'mutChi2', 'chi2DF', 'wtPval', 'mutPval']
            fmts += ['%g', '%g', '%d', '%g', '%g']
            outs += [wtD, mutD, mutDF, wtP, mutP]


        if True:
            wtKL = kullbackLeibler(xh, zh)
            mutKL = kullbackLeibler(yh, zh)

            kls += [wtKL, mutKL]

            hdrs += ['wtKL', 'mutKL']
            fmts += ['%g', '%g']
            outs += [wtKL, mutKL]

        if True:
            wtKS = kolmogorovSmirnov(xh, zh)[1]
            mutKS = kolmogorovSmirnov(yh, zh)[1]

            hdrs += ['wtKS', 'mutKS']
            fmts += ['%g', '%g']
            outs += [wtKS, mutKS]

        if not opts['-e'] and opts['-l']:
            vx = {}
            #for (c,f) in zh0.iteritems():
            #    if c not in vx:
            #        vx[c] = [0, 0, 0]
            #    vx[c][0] = f
            for (c,f) in xh0.iteritems():
                if c not in vx:
                    vx[c] = [0, 0, 0]
                vx[c][1] = f
            for (c,f) in yh0.iteritems():
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
                hdrs1 = hdrs + ['coverage', 'wtCnt', 'wtCum', 'mutCnt', 'mutCum']
                fmts1 = fmts + ['%d', '%d', '%g', '%d', '%g']
                outs1 = outs + [c, fs[1], tx[1]/nx, fs[2], tx[2]/ny]
                if not hdrShown:
                    hdrShown = True
                    print '\t'.join(hdrs1)
                print '\t'.join(fmts1) % tuple(outs1)
        elif not opts['-e'] and opts['-a']:
            if not hdrShown:
                hdrShown = True
                print '\t'.join(hdrs)
            print '\t'.join(fmts) % tuple(outs)
        sys.stdout.flush()

    if opts['-e']:
        kls.sort()
        m = len(kls) // 2
        p = 0.01
        negCut = estimateGammaCutoff(kls[:m], 1-p)
        posCut = estimateGammaCutoff(kls[m:], p)
        if opts['-v']:
            print 'estimated cutoff to reject a negative call at p < %g is %g' % (p, negCut)
            print 'estimated cutoff to reject a positive call at p < %g is %g' % (p, posCut)
        else:
            print p, negCut, posCut

if __name__ == '__main__':
    main(sys.argv[1:])
