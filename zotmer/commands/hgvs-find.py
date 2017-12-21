"""
zot hgvs-find - search for known variants in read data.

Usage:
    zot hgvs-find -X [options] <index> [<variant>...]
    zot hgvs-find [options] <index> <input>...

Options:
    -X              index HGVS variants
    -A ALPHA        alpha level for Kolmogorov-Smirnov test [default: 0.001]
    -B BINFUNC      binning function to use [default: none]
    -S SCALE        scaling factor for binning function [default: 1.0]
    -f FILENAME     read variants from a file
    -F FORMAT       a format string for printing results [default: ks]
    -k K            value of k to use [default: 25]
    -g PATH         directory of FASTQ reference sequences
    -s              single stranded k-mer extraction
    -v              produce verbose output

Format Strings

Format strings control the fields that get included in the output.
The value is a comma separated list. The order of fields is fixed,
so the order in the format string is irrelevant.  The supported
fields are:

    ks              Kolmogorov-Smirnov statistics
    hist            produce k-mer frequency histograms
"""

import array
import math
import random
import sys
import zipfile

import docopt
import tqdm
import yaml

from pykmer.basics import ham, kmers, kmersList, kmersWithPosList, murmer, rc, render
from pykmer.file import openFile, readFasta, readFastq
from pykmer.stats import counts2cdf, pdf2cdf, logBinEq, logBinLe, ksDistance2, logGammaP, logGammaQ, logLowerGamma
from zotmer.library.debruijn import interpolate
from zotmer.library.hgvs import makeHGVS, refSeq2Hg19
from zotmer.library.reads import reads

AluY = 'GGCCAGGAGTGGTGGCTCAGCCTATAATCCCAGCACTTTGGGAGGCTGAGGCAGGCAGATCACGAGGTCAGGAGATCGAGACCACCCTGGCTAACATGGTGAAACCCTGTCTCCACTAAAAATAGAAAAAAATTAGCATGGTGTGGTGGCATGCACCTGTAGTCCCAGCTGTTAGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCAGAGAGGCAGAGATTGCAGTGAGCCAAGATCGTGCCACTGCACTCCAGCCAGGGCGACAGAGTGAGACTCCGTCTAAA'

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

class summarizer(object):
    def __init__(self):
        self.n = 0
        self.s = 0
        self.s2 = 0

    def add(self, x):
        self.n += 1
        self.s += x
        self.s2 += x * x

    def mean(self):
        return float(self.s)/float(self.n)

    def var(self):
        m = self.mean()
        return float(self.s2)/float(self.n) - m*m

    def sd(self):
        return math.sqrt(self.var())

def computeBias(K, zs, verbose = False):
    S = summarizer()
    for (x, xc) in zs.iteritems():
        y = rc(K, x)
        if y < x:
            continue
        yc = zs.get(y, 0)

        if xc > yc:
            a = xc
            b = yc
        else:
            a = yc
            b = xc
        apb = a + b
        if apb > 0:
            v = float(a) / float(apb)
        else:
            v = 0.5
        if verbose:
            print '%s\t%s\t%d\t%d\t%g' % (render(K, x), render(K, y), xc, yc, v)
        S.add(v)
    return (S.mean(), S.var())

def nullModelPdf(K):
    return [math.exp(logBinEq(0.25, K, j)) for j in range(K+1)]

def binLe(p, N, k):
    return math.exp(logBinLe(p, N, k))

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

def histQuantiles(Xs, n):
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

def binLin(Xs, S):
    r = {}
    for (x, c) in Xs.items():
        y = int(x*S)
        if y not in r:
            r[y] = 0
        r[y] += c
    return r

def binLog(Xs, S):
    r = {}
    for (x, c) in Xs.items():
        y = int(math.log(1+x)*S)
        if y not in r:
            r[y] = 0
        r[y] += c
    return r

def binSqrt(Xs, S):
    r = {}
    for (x, c) in Xs.items():
        y = int(math.sqrt(x)*S)
        if y not in r:
            r[y] = 0
        r[y] += c
    return r

def binHist(Xs, tx, S):
    if tx == 'log':
        return binLog(Xs, S)
    if tx == 'sqrt':
        return binSqrt(Xs, S)
    return binLin(Xs, S)

def histMean(hx):
    s = 0
    n = 0
    for (c,f) in hx.items():
        s += c*f
        n += f
    n = max(1,n)
    return float(s) / float(n)

def histMedian(hx):
    t = float(max(1, sum([c*f for (c,f) in hx.items()])))
    s = 0
    for (c,f) in hx.items():
        s += c*f
        if s / t >= 0.5:
            return c
    return 0

def makeRanks(xns):
    n = float(len(xns))
    r = [None for v in xns]
    for i in range(len(xns)):
        (v,j) = xns[i]
        r[j] = (i/n,v)
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

def wassersteinMetric(hx, hy, M = None):
    zs = list(set(hx.keys() + hy.keys()))
    zs.sort()

    if M is None:
        M = float(zs[-1])

    cx = 0
    tx = float(sum(hx.values()))
    cy = 0
    ty = float(sum(hy.values()))

    a = 0.0
    pz = 0
    for z in zs:
        if z in hx:
            cx += hx[z]
        if z in hy:
            cy += hy[z]
        px = cx / tx
        py = cy / ty
        d = max(0, py - px)
        a += d * (z - pz)
        pz = z
    return a / M

def differentialArea(hx, hy):
    zs = list(set(hx.keys() + hy.keys()))
    zs.sort()

    cx = 0
    tx = float(sum(hx.values()))
    cy = 0
    ty = float(sum(hy.values()))

    a = 0.0
    ppy = 0.0
    for z in zs:
        if z in hx:
            cx += hx[z]
        if z in hy:
            cy += hy[z]
        px = cx / tx
        py = cy / ty
        a += px * (py - ppy)
        ppy = py
    return a

def kullbackLeibler(Xs, Ys):
    bs = Ys.keys()
    bs.sort()

    zx = float(sum(Xs.values()))
    zx = max(1.0, zx)
    px = {}
    for b in bs:
        c = Xs.get(b, 0)
        px[b] = c/zx

    zy = float(sum(Ys.values()))
    zy = max(1.0, zy)
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

def gammaPDF(m, x):
    (k, t) = m
    num = math.log(x) * (k - 1) - x/t
    den = math.lgamma(k) + k * math.log(t)
    return math.exp(num - den)

def gammaCDF(m, x):
    (k, t) = m
    gk = math.lgamma(k)
    return math.exp(logLowerGamma(k, x/t) - gk)

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

def posKmer(K, x, p):
    return x | (p << (2*K))

def unpackPosKmer(K, v):
    M = (1 << (2*K)) - 1
    x = v&M
    p = v >> (2*K)
    return (x, p)

def renderPath(K, xs):
    if len(xs) == 0:
        return ''
    r = [render(K, xs[0])]
    for x in xs[1:]:
        r.append("ACGT"[x&3])
    return ''.join(r)

def scoreAllele(K, pthRes, lhs, seq, rhs):
    pth = pthRes[2]

    r = [0 for i in range(K+1)]

    if seq is None:
        r[0] += len(pth)
        return r

    xs = kmersList(K, lhs[-(K+pthRes[3]-1):] + seq + rhs[:(K+pthRes[4]-1)])

    for d in [ham(x,y) for (x,y) in zip(pth, xs)]:
        r[d] += 1

    return r

def makeIndexedVariant(v, K):
    r = {}
    r['hgvs'] = str(v)
    rng = v.range()
    big = v.loadAccession()
    lhsSt = rng[0] - 2*K
    lhsEn = rng[0]
    lhs = big[lhsSt:lhsEn]
    rhsSt = rng[1]
    rhsEn = rng[1] + 2*K
    rhs = big[rhsSt:rhsEn]
    r['lhsFlank'] = lhs
    r['rhsFlank'] = rhs
    r['wtSeq'] = big[rng[0]:rng[1]]
    if v.anonymous() or v.cryptic():
        r['mutSeq'] = None
    else:
        r['mutSeq'] = v.sequence()

    return r

def checkVariant(v, K):
    (wtSeq, mutSeq) = v.context(K)

    # A regular variant
    wtXs = set(kmers(K, wtSeq, False))
    mutXs = set(kmers(K, mutSeq, False))
    com = wtXs & mutXs
    wtXs -= com
    mutXs -= com
    if len(wtXs) == 0 or len(mutXs) == 0:
        print >> sys.stderr, "variant has no distinguishing k-mers: %s" % (str(v),)
        print >> sys.stderr, wtSeq
        print >> sys.stderr, mutSeq
        return False
    return True

def addProbesToIdx(n, xs, idx):
    for x in xs:
        if x not in idx:
            idx[x] = set([])
        idx[x].add(n)

def makeMetaIdx(idx):
    top = {}
    bot = {}
    for (x, ys) in idx.iteritems():
        h = 17
        for y in sorted(ys):
            h = murmer(y, h)
        top[x] = h
        bot[h] = ys
    return (top, bot)

def findHits(idx, xs, rx):
    for x in xs:
        if x in idx:
            v = idx[x]
            if len(rx) == 0 or rx[-1] != v:
                rx.append(v)

def debruijn(K, x, y):
    M = (1 << 2*(K-1)) - 1
    return (x&M) == (y >> 2)

def findAnchors(K, seq, mx, isLhs, D):
    xps = kmersWithPosList(K, seq, False)
    xps = [(x, p-1) for (x,p) in xps]

    # Find the highest coverage k-mer that intersects.
    xc = 0
    for (x,p) in xps:
        if x in mx and mx[x] >= xc:
            xc = mx[x]

    if xc == 0:
        return set([])

    # Seeds should be in the same order of magnitude
    # as the highest-coverage seed.
    t = int(math.exp(math.log(xc) - 1.5))

    xs = {}
    for (x,c) in mx.iteritems():
        if c < t:
            continue
        xs[x] = c

    ys = xs.keys()

    zs = {}
    for (x,p) in xps:
        zs[p] = set([])

        if x not in xs:
            continue

        for y in ys:
            d = ham(x, y)
            if d > D:
                continue
            zs[p].add(y)

    e = set([])
    res = set([])
    if isLhs:
        for (x,p) in xps:
            ss = zs.get(p, e)
            tt = zs.get(p-1, e)
            for s in ss:
                res.add((s,p))
                for t in tt:
                    if debruijn(K, t, s):
                        res.discard((t,p-1))
    else:
        for (x,p) in xps[::-1]:
            ss = zs.get(p, e)
            tt = zs.get(p+1, e)
            for s in ss:
                res.add((s,p))
                for t in tt:
                    if debruijn(K, s, t):
                        res.discard((t,p+1))

    if isLhs:
        l = len(seq) - K
        res = [(x, l - p) for (x,p) in res]
    else:
        res = list(res)
    res.sort()

    return res

def findAllele(K, mx, lx, lxp, rx, rxp, z):
    ll = z + lxp + rxp + 1 + K
    pth = interpolate(K, mx, lx, rx, ll)
    if pth is None:
        return None
    resSeq = renderPath(K, pth)
    resSeq = resSeq[lxp+K:-(rxp+K)]
    pthMin = min([mx[x] for x in pth[(lxp+1):-(rxp+1)]])
    return (resSeq, pthMin, pth[(lxp+1):-(rxp+1)], lxp, rxp)

class SequenceFactory(object):
    def __init__(self, home):
        self.home = home
        self.prevAcc = None
        self.prevSeq = None

    def __getitem__(self, acc):
        if acc != self.prevAcc:
            if acc in refSeq2Hg19:
                h = refSeq2Hg19[acc]
            else:
                h = acc

            with openFile(self.home + "/" + h + ".fa.gz") as f:
                for (nm,seq) in readFasta(f):
                    self.prevAcc = acc
                    self.prevSeq = seq
                    break
        return self.prevSeq

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = int(opts['-k'])

    verbose = opts['-v']

    d = "."
    if opts['-g']:
        d = opts['-g']
    sf = SequenceFactory(d)

    if opts['-X']:
        variants = opts['<variant>']
        if opts['-f']:
            with openFile(opts['-f']) as f:
                variants += f.read().split()

        vx = {}
        for v in variants:
            x = makeHGVS(v)
            if x is None:
                print >> sys.stderr, "unable to parse %s" % (v,)
                continue
            x.setSequenceFactory(sf)
            acc = x.accession()
            if acc not in vx:
                vx[acc] = []
            vx[acc].append(x)

        rs = []
        for (acc, vs) in vx.iteritems():
            for v in vs:
                checkVariant(v, K)
                r = makeIndexedVariant(v, K)
                if r is not None:
                    rs.append(r)

        with open(opts['<index>'], 'w') as f:
            yaml.safe_dump(rs, f, default_flow_style=False)

        return

    fmt = set(opts['-F'].split(','))

    if verbose:
        print >> sys.stderr, "loading index."

    with open(opts['<index>']) as f:
        hgvsVars = yaml.load(f)

    V = len(hgvsVars)

    probeIdx = {}

    for n in range(V):
        itm = hgvsVars[n]
        h = itm['hgvs']
        v = makeHGVS(h)
        itm['var'] = v
        probes = set(kmersList(K, itm['lhsFlank'], True) + kmersList(K, itm['rhsFlank'], True))
        addProbesToIdx(n, probes, probeIdx)

    (probeHashes, hashSets) = makeMetaIdx(probeIdx)

    if verbose:
        print >> sys.stderr, "done."

    combineStrands = True
    if opts['-s']:
        combineStrands = False

    kmerHits = [{} for n in range(V)]

    rn = 0
    #for itm in reads(opts['<input>'], K=K, paired=False, reads=True, kmers=True, both=True, verbose=verbose):
    for itm in reads(opts['<input>'], K=K, paired=True, reads=True, kmers=True, both=True, verbose=verbose):
        rn += 1
        probeHits = []
        for xs in itm.kmers:
            findHits(probeHashes, xs, probeHits)
            rn = -rn

        if len(probeHits) == 0:
            continue

        setHits = set([])
        for h in set(probeHits):
            setHits |= hashSets[h]

        for n in setHits:
            for xs in itm.kmers:
                for x in xs:
                    kmerHits[n][x] = 1 + kmerHits[n].get(x, 0)

    D = 4

    Q = 10

    nullCdf = pdf2cdf(nullModelPdf(K))

    ps = [1-binLe(0.25, K, d) for d in range(K+1)]

    hdrShown = False
    for n in range(V):
        itm = hgvsVars[n]
        v = itm['var']
        h = itm['hgvs']

        mx = kmerHits[n]

        lhsFlank = itm['lhsFlank']
        rhsFlank = itm['rhsFlank']

        lhs = findAnchors(K, lhsFlank, mx, True, D)
        rhs = findAnchors(K, rhsFlank, mx, False, D)

        if len(lhs)*len(rhs) > 10:
            print >> sys.stderr, 'warning: highly ambiguous anchors for', str(v)
            print >> sys.stderr, '%d\tlhs anchors, and %d rhs anchors' % (len(lhs), len(rhs))

        wtSeq = itm['wtSeq']
        wtZ = len(wtSeq)

        mutSeq = itm['mutSeq']
        mutZ = v.size()

        alleles = {}
        alleles['wt'] = []
        alleles['mut'] = []
        for (lx, lxp) in lhs:
            for (rx, rxp) in rhs:
                if wtZ == mutZ:
                    pthRes =  findAllele(K, mx, lx, lxp, rx, rxp, wtZ)
                    if pthRes is None:
                        continue
                    pthSeq = pthRes[0]
                    if pthSeq == wtSeq:
                        alleles['wt'].append(pthRes)
                        continue
                    if mutSeq is None:
                        alleles['mut'].append(pthRes)
                    elif pthSeq == mutSeq:
                        alleles['mut'].append(pthRes)
                else:
                    pthRes =  findAllele(K, mx, lx, lxp, rx, rxp, wtZ)
                    if pthRes is not None and pthRes[0] == wtSeq:
                        alleles['wt'].append(pthRes)

                    pthRes =  findAllele(K, mx, lx, lxp, rx, rxp, mutZ)
                    if pthRes is None:
                        continue
                    pthSeq = pthRes[0]
                    if mutSeq is None:
                        alleles['mut'].append(pthRes)
                    elif pthSeq == mutSeq:
                        alleles['mut'].append(pthRes)

        wtMin = 0
        wtHs = [0 for j in range(K+1)]
        for pthRes in alleles['wt']:
            wtA =  scoreAllele(K, pthRes, lhsFlank, wtSeq, rhsFlank)
            if pthRes[1] >= wtMin:
                wtMin = pthRes[1]
                wtHs = wtA
        wtCdf = counts2cdf(wtHs)
        wtD = ksDistance2(wtCdf, nullCdf)[0]
        wtQ = 1
        for d in range(K+1):
            wtQ *= math.pow(ps[d], wtHs[d])
        wtQ = 1 - wtQ

        mutMin = 0
        mutHs = [0 for j in range(K+1)]
        for pthRes in alleles['mut']:
            mutA =  scoreAllele(K, pthRes, lhsFlank, mutSeq, rhsFlank)
            if pthRes[1] >= mutMin:
                mutMin = pthRes[1]
                mutHs = mutA
        mutCdf = counts2cdf(mutHs)
        mutD = ksDistance2(mutCdf, nullCdf)[0]
        mutQ = 1
        for d in range(K+1):
            mutQ *= math.pow(ps[d], mutHs[d])
        mutQ = 1 - mutQ

        hdrs = ['n']
        fmts = ['%d']
        outs = [n]

        wtAllele = ((wtMin > Q) and (wtQ < 0.05))
        mutAllele = ((mutMin > Q) and (mutQ < 0.05))
        resV =  1*wtAllele + 2*mutAllele
        res = ['null', 'wt', 'mut', 'wt/mut'][resV]

        hdrs += ['res']
        fmts += ['%s']
        outs += [res]

        hdrs += ['wtMin', 'mutMin']
        fmts += ['%d', '%d']
        outs += [wtMin, mutMin]

        hdrs += ['wtD', 'mutD']
        fmts += ['%g', '%g']
        outs += [wtD, mutD]

        hdrs += ['wtQ', 'mutQ']
        fmts += ['%g', '%g']
        outs += [wtQ, mutQ]

        hdrs += ['hgvs']
        fmts += ['%s']
        outs += [h]

        if not hdrShown:
            hdrShown = True
            print '\t'.join(hdrs)
        print '\t'.join(fmts) % tuple(outs)

if __name__ == '__main__':
    main(sys.argv[1:])
