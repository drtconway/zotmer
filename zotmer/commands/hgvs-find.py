"""
zot hgvs-find - search for known variants in read data.

Usage:
    zot hgvs-find -X [options] <index> [<variant>...]
    zot hgvs-find [options] <index> <input>...

Options:
    -T BEDFILE      test regions at indexing time
    -X              index HGVS variants
    -A ALPHA        alpha level for Kolmogorov-Smirnov test [default: 0.001]
    -B BINFUNC      binning function to use [default: none]
    -C MINCOV       minimum coverage to call an allele [default: 10]
    -c FILENAME     capture reads into named zipfile
    -D MAXHAM       maximum hamming distance to allow in flanks for calling an allele [default: 2]
    -f FILENAME     read variants from a file
    -F FORMAT       a format string for printing extra result statistics [default: rds,vaf]
    -k K            value of k to use [default: 25]
    -g PATH         directory of FASTQ reference sequences
    -o FILENMANE    write the output to the named file rather than stdout
    -s              single stranded k-mer extraction
    -v              produce verbose output
    -V VAF          minimum vaf for calling an allele [default: 0.05]
    -w WIDTH        [indexing] width of context for capture [default: 75]
    -W WIDTH        [indexing] width of context for verification [default: 1024]

Format Strings

Format strings control the fields that get included in the output.
The value is a comma separated list. The order of fields is fixed,
so the order in the format string is irrelevant.  The supported
fields are:

    binom           Score alleles against a binomial model for Hamming distances
    hist            produce k-mer frequency histograms
    ks              Kolmogorov-Smirnov statistics
    rds             count the number of reads captured for each variant
    vaf             compute an approximate VAF
"""

import array
import math
import random
import sys
import zipfile

import docopt
import yaml

from pykmer.basics import ham, kmers, kmersList, kmersWithPosList, murmer, rc, render
from pykmer.file import openFile, readFasta, readFastq
from pykmer.misc import unionfind
from pykmer.stats import counts2cdf, pdf2cdf, logBinEq, logBinLe, ksDistance2, logGammaP, logGammaQ, logLowerGamma
from zotmer.library.capture import capture
from zotmer.library.debruijn import fixed_path_kmers, paths
from zotmer.library.hgvs import hg19ToRefSeq, makeHGVS, refSeq2Hg19
from zotmer.library.reads import reads

def sqr(x):
    return x*x

def outputFile(nm):
    if nm is None:
        return sys.stdout
    return openFile(nm, 'w')

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

def hamming(p, q):
    i = 0
    z = len(p)
    d = 0
    while i < z:
        if p[i] != q[i] and p[i] != 'N' and q[i] != 'N':
            d += 1
        i += 1
    return d

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

def readBED(f):
    res = {}
    first = True
    for l in f:
        t = l.split()
        if first:
            first = False
            if t[0] == 'track' or t[0] =='browser':
                continue
        c = t[0]
        if c in hg19ToRefSeq:
            c = hg19ToRefSeq[c]
        s = int(t[1])
        e = int(t[2])
        n = None
        if len(t) > 3:
            n = t[3]
        v = (s, e, n)
        if c not in res:
            res[c] = []
        res[c].append(v)
    return res

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

def context(r, W, seq):
    Z = len(seq)
    lhsSt = r[0] - W
    lhsSt = max(0, lhsSt)
    lhsEn = r[0]
    lhs = seq[lhsSt:lhsEn].upper()
    rhsSt = r[1]
    rhsEn = r[1] + W
    rhsEn = min(rhsEn, Z)
    rhs = seq[rhsSt:rhsEn].upper()
    return (lhs,rhs)

def makeIndexedVariant(v, K, W, V):
    r = {}
    r['hgvs'] = str(v)
    rng = v.range()
    seq = v.loadAccession()

    (lhs, rhs) = context(rng, W, seq)
    r['lhsFlank'] = lhs
    r['rhsFlank'] = rhs

    r['wtSeq'] = seq[rng[0]:rng[1]].upper()
    if v.anonymous() or v.cryptic():
        r['mutSeq'] = None
    else:
        r['mutSeq'] = v.sequence().upper()

    #(lhsVal, rhsVal) = context(rng, V, seq)
    #r['lhsVal'] = lhsVal
    #r['rhsVal'] = rhsVal

    return r

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

def consAllele(lhs, mid, rhs, lp, rp):
    if mid is None:
        return None
    return lhs[-lp:] + mid + rhs[:rp]

def consensus_kmer(K, xs, X):
    R = [[0,0,0,0] for i in range(K)]
    for x in xs:
        c = X[x]
        for i in range(K):
            R[i][x&3] += c
            x >>= 2
    c = None
    r = 0
    for i in range(K):
        cMax = -1
        bMax = 0
        for b in range(4):
            if R[i][b] > cMax:
                bMax = b
                cMax = R[i][b]
        r |= bMax << (2*i)
        if c is None:
            c = cMax
        else:
            c = min(c, cMax)
    return (r, c)

def findSimpleAllele(K, lhs, mid, rhs, mx):
    allele = consAllele(lhs, mid, rhs, K-1, K-1)
    xs = kmersList(K, allele)
    m = min([mx.get(x, 0) for x in xs])
    res = {}
    res['allele'] = mid
    res['covMin'] = m
    res['path'] = xs
    res['lhsPos'] = 0
    res['rhsPos'] = 0
    return res

def findFairlySimpleAllele(K, D, lhs, mid, rhs, mx):
    allele = consAllele(lhs, mid, rhs, K-1, K-1)
    xs = kmersList(K, allele)

    Y = fixed_path_kmers(K, D, mx, xs)
    if Y is None:
        return

    ancPath = []
    cov = []
    for i in range(len(xs)):
        ys = Y[i]
        (y,c) = consensus_kmer(K, ys, mx)
        ancPath.append(y)
        cov.append(c)

    seq = renderPath(K, ancPath)

    res = {}
    res['ancPath'] = ancPath
    res['allele'] = seq
    res['covMin'] = min(cov)
    res['path'] = ancPath
    res['lhsPos'] = 0
    res['rhsPos'] = 0
    yield res

def findAllele(K, mx, lx, lxp, rx, rxp, z):
    ll = z + lxp + rxp + 1 + K
    for pth in paths(K, mx, lx, rx, ll):
        seq = renderPath(K, pth)
        seq = seq[lxp+K:-(rxp+K)]
        res = {}
        res['ancPath'] = pth
        res['allele'] = seq
        res['covMin'] = min([mx[x] for x in pth[(lxp+1):-(rxp+1)]])
        res['path'] = pth[1:-1]
        res['lhsPos'] = lxp
        res['rhsPos'] = rxp
        yield res

def findPath(K, mx, lx, lxp, rx, rxp):
    for pth in paths(K, mx, lx, rx, slice(1, 250)):
        seq = renderPath(K, pth)
        seq = seq[lxp+K:-(rxp+K)]
        res = {}
        res['ancPath'] = pth
        res['allele'] = seq
        res['covMin'] = min([mx[x] for x in pth[(lxp+1):-(rxp+1)]])
        res['path'] = pth[1:-1]
        res['lhsPos'] = lxp
        res['rhsPos'] = rxp
        yield res

class Scorer(object):
    def __init__(self, K):
        self.K = K
        self.nullCdf = pdf2cdf(nullModelPdf(K))
        self.binModel =  [1-binLe(0.25, K, d) for d in range(K+1)]

    def score(self, res, lhs, seq, rhs):
        Km1 = self.K - 1
        Kp1 = self.K + 1

        xs = res['path']

        r = [0 for i in range(Kp1)]

        if seq is None:
            n = len(res['allele'])
            allele = consAllele(lhs, n*'N', rhs, Km1+res['lhsPos'], Km1+res['rhsPos'])
            r[0] += len(xs)
        else:
            allele = consAllele(lhs, seq, rhs, Km1+res['lhsPos'], Km1+res['rhsPos'])
            ys = kmersList(self.K, allele)
            for d in [ham(x,y) for (x,y) in zip(xs, ys)]:
                r[d] += 1

        res['hammingProfile'] = r
        res['hammingCdf'] = counts2cdf(r)
        res['ksDist'] = ksDistance2(res['hammingCdf'], self.nullCdf)[0]
        res['hamming'] = hamming(allele, renderPath(self.K, xs))

        p = 1.0
        for d in xrange(Kp1):
            p *= math.pow(self.binModel[d], r[d])
        res['binom'] = 1.0 - p

def isBetter(lhs, rhs):
    if lhs['covMin'] < rhs['covMin']:
        return False
    if lhs['covMin'] > rhs['covMin']:
        return True
    if lhs['hamming'] > rhs['hamming']:
        return False
    return True

def cat(xss):
    for xs in xss:
        for x in xs:
            yield x

class AlleleFinder(object):
    def __init__(self, K, D, v, mx, lhsFlank, rhsFlank, wtSeq, mutSeq, wtZ, mutZ):
        self.K = K
        self.M = ((1 << (2*K)) - 1)
        self.S = 2*(K-1)
        self.D = D
        self.v = v
        self.mx = mx
        self.lhsFlank = lhsFlank
        self.rhsFlank = rhsFlank
        self.wtSeq = wtSeq
        self.mutSeq = mutSeq
        self.wtZ = wtZ
        self.mutZ = mutZ

    def simpleAlleles(self):
        ax = findSimpleAllele(self.K, self.lhsFlank, self.wtSeq, self.rhsFlank, self.mx)
        if ax is not None:
            yield ('wt', ax)
        if self.mutSeq is not None:
            ax = findSimpleAllele(self.K, self.lhsFlank, self.mutSeq, self.rhsFlank, self.mx)
            if ax is not None:
                yield ('mut', ax)

    def definiteAlleles(self):
        for ax in findFairlySimpleAllele(self.K, self.D, self.lhsFlank, self.wtSeq, self.rhsFlank, self.mx):
            yield ('wt', ax)
        if self.mutSeq is not None:
            for ax in findFairlySimpleAllele(self.K, self.D, self.lhsFlank, self.mutSeq, self.rhsFlank, self.mx):
                yield ('mut', ax)

    def bridgingAlleles(self):
        lhs = findAnchors(self.K, self.lhsFlank, self.mx, True, self.D)
        rhs = findAnchors(self.K, self.rhsFlank, self.mx, False, self.D)

        if len(lhs)*len(rhs) > 10:
            print >> sys.stderr, 'warning: highly ambiguous anchors for', str(self.v)
            print >> sys.stderr, '%d\tlhs anchors, and %d rhs anchors' % (len(lhs), len(rhs))

        maxBridge = 4*self.K

        for (lx, lxp) in lhs:
            for (rx, rxp) in rhs:
                if self.wtZ == self.mutZ:
                    if self.wtZ > maxBridge:
                        print >> sys.stderr, 'warning: cannot bridge alleles because they are too long (%d) for %s' % (self.wtZ, self.v)
                        continue
                    for pthRes in findAllele(self.K, self.mx, lx, lxp, rx, rxp, self.wtZ):
                        pthSeq = pthRes['allele']
                        if pthSeq == self.wtSeq:
                            yield ('wt', pthRes)
                            continue
                        if self.mutSeq is None or pthSeq == self.mutSeq:
                            yield ('mut', pthRes)
                else:
                    if self.wtZ > maxBridge:
                        print >> sys.stderr, 'warning: cannot bridge wt allele because it is too long (%d) for %s' % (self.wtZ, self.v)
                        continue
                    else:
                        for pthRes in findAllele(self.K, self.mx, lx, lxp, rx, rxp, self.wtZ):
                            pthSeq = pthRes['allele']
                            if pthSeq == self.wtSeq:
                                yield ('wt', pthRes)

                    if self.mutZ > maxBridge:
                        print >> sys.stderr, 'warning: cannot bridge mut allele because it is too long (%d) for %s' % (self.mutZ, self.v)
                        continue
                    else:
                        for pthRes in findAllele(self.K, self.mx, lx, lxp, rx, rxp, self.mutZ):
                            pthSeq = pthRes['allele']
                            if self.mutSeq is None or pthSeq == self.mutSeq:
                                yield ('mut', pthRes)

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

    verbose = opts['-v']

    K = int(opts['-k'])

    D = int(opts['-D'])

    Q = int(opts['-C'])

    V = float(opts['-V'])

    d = "."
    if opts['-g']:
        d = opts['-g']
    sf = SequenceFactory(d)

    if opts['-X']:
        Wcap = int(opts['-w'])
        Wval = int(opts['-W'])

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
                r = makeIndexedVariant(v, K, Wcap, Wval)
                if r is not None:
                    rs.append(r)

        with open(opts['<index>'], 'w') as f:
            yaml.safe_dump(rs, f, default_flow_style=False)

        return

    capt = False
    zipname = None
    if opts['-c']:
        capt = True
        zipname = opts['-c']

    fmt = set([])
    if opts['-F']:
        fmt = set(opts['-F'].split(','))

    if verbose:
        print >> sys.stderr, "loading index."

    with open(opts['<index>']) as f:
        hgvsVars = yaml.load(f)

    NV = len(hgvsVars)

    combineStrands = True
    if opts['-s']:
        combineStrands = False

    cap = capture(K, reads=capt, kmers=True, verbose=verbose)

    for n in range(NV):
        itm = hgvsVars[n]
        h = itm['hgvs']
        v = makeHGVS(h)
        itm['var'] = v
        lhs = itm['lhsFlank']
        rhs = itm['rhsFlank']
        wt = itm['wtSeq']
        mut = itm['mutSeq']
        bait = [lhs, wt, rhs]
        if mut is not None:
            bait += ['N']
            bait += [lhs, mut, rhs]
        bait = ''.join(bait)
        n0 = cap.addBait(h, bait)
        assert n0 == n

    if verbose:
        print >> sys.stderr, "done."

    rn = 0
    for itm in reads(opts['<input>'], K=K, paired=True, reads=True, kmers=False, both=True, verbose=verbose):
        rn += 1
        cap.addReadPairAndKmers(itm.reads[0], itm.reads[1])

    if capt:
        cap.saveReads(zipname)

    scorer = Scorer(K)

    globHist = {}

    for n in range(NV):
        mx = cap.capKmers[n]
        for c in mx.itervalues():
            if c < Q:
                continue
            if c not in globHist:
                globHist[c] = 0
            globHist[c] += 1

    with outputFile(opts['-o']) as out:
        hdrShown = False
        for n in range(NV):
            itm = hgvsVars[n]
            v = itm['var']
            h = itm['hgvs']

            mx = cap.capKmers[n]

            nr = cap.capReadCounts[n]

            if 'kmers' in fmt:
                for (x,c) in mx.iteritems():
                    print '%d\t%s\t%d' % (n, render(K, x), c)

            lhsFlank = itm['lhsFlank']
            rhsFlank = itm['rhsFlank']

            alleles = {}
            alleles['wt'] = []
            alleles['mut'] = []

            wtSeq = itm['wtSeq']
            wtZ = len(wtSeq)

            mutSeq = itm['mutSeq']
            mutZ = v.size()

            cs = [c for (x,c) in mx.iteritems() if c >= Q]
            cs.sort()
            nk = len(cs)
            if nk == 0:
                cs = [0]

            q10 = cs[1*len(cs)//10]
            q50 = cs[5*len(cs)//10]
            q90 = cs[9*len(cs)//10]

            af = AlleleFinder(K, D, v, mx, lhsFlank, rhsFlank, wtSeq, mutSeq, wtZ, mutZ)
            finders = []
            if not v.anonymous():
                finders.append(af.definiteAlleles())
            else:
                finders.append(af.bridgingAlleles())

            j = 0
            for (t, a) in cat(finders):
                assert t == 'wt' or t == 'mut'
                alleles[t].append(a)
                j += 1

            wtRes = {}
            wtRes['covMin'] = 0
            wtRes['binom'] = 1.0
            wtRes['ksDist'] = 0.0
            wtRes['hamming'] = 0
            wtRes['path'] = []
            for pthRes in alleles['wt']:
                scorer.score(pthRes, lhsFlank, wtSeq, rhsFlank)
                if isBetter(pthRes, wtRes):
                    wtRes = pthRes

            mutRes = {}
            mutRes['covMin'] = 0
            mutRes['binom'] = 1.0
            mutRes['ksDist'] = 0.0
            mutRes['hamming'] = 0
            mutRes['path'] = []
            for pthRes in alleles['mut']:
                scorer.score(pthRes, lhsFlank, mutSeq, rhsFlank)
                if isBetter(pthRes, mutRes):
                    mutRes = pthRes

            if True:
                wtXs = [mx.get(x, 0) for x in wtRes['path']]
                if len(wtXs) == 0:
                    wtXs = [0]
                wtXs.sort()
                wtCount = sum(wtXs)
                wtLen = len(wtXs)
                wtMean = float(wtCount)/float(wtLen)
                wtMedian = wtXs[wtLen//2]

                mutXs = [mx.get(x, 0) for x in mutRes['path']]
                if len(mutXs) == 0:
                    mutXs = [0]
                mutXs.sort()
                mutCount = sum(mutXs)
                mutLen = len(mutXs)
                mutMean = float(mutCount)/float(mutLen)
                mutMedian = mutXs[mutLen//2]

                totX  = max([1.0, float(wtMedian+mutMedian), float(q90)])
                wtVaf = wtMedian/totX
                mutVaf = mutMedian/totX

            hdrs = ['n']
            fmts = ['%d']
            outs = [n]

            wtAllele = ((wtRes['covMin'] > Q) and (wtRes['hamming'] < 4)) and (wtVaf > V)
            mutAllele = ((mutRes['covMin'] > Q) and (mutRes['hamming'] < 4)) and (mutVaf > V)
            resV = 1*wtAllele + 2*mutAllele
            res = ['null', 'wt', 'mut', 'wt/mut'][resV]

            hdrs += ['res']
            fmts += ['%s']
            outs += [res]

            if 'rds' in fmt:
                hdrs += ['numReads']
                fmts += ['%d']
                outs += [nr]

            hdrs += ['numKmers', 'covQ10', 'covQ50', 'covQ90']
            fmts += ['%d', '%d', '%d', '%d']
            outs += [nk, q10, q50, q90]

            hdrs += ['wtMin', 'mutMin']
            fmts += ['%d', '%d']
            outs += [wtRes['covMin'], mutRes['covMin']]

            hdrs += ['wtHam', 'mutHam']
            fmts += ['%d', '%d']
            outs += [wtRes['hamming'], mutRes['hamming']]

            if 'ks' in fmt:
                hdrs += ['wtD', 'mutD']
                fmts += ['%g', '%g']
                outs += [wtRes['ksDist'], mutRes['ksDist']]

            if 'binom' in fmt:
                hdrs += ['wtQ', 'mutQ']
                fmts += ['%g', '%g']
                outs += [wtRes['binom'], mutRes['binom']]

            if 'vaf' in fmt:
                hdrs += ['wtVaf', 'mutVaf']
                fmts += ['%g', '%g']
                outs += [wtVaf, mutVaf]

            hdrs += ['hgvs']
            fmts += ['%s']
            outs += [h]

            if not hdrShown:
                hdrShown = True
                print >> out, '\t'.join(hdrs)
            print >> out, '\t'.join(fmts) % tuple(outs)
            out.flush()

if __name__ == '__main__':
    main(sys.argv[1:])
