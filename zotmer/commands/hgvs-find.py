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
    -L LENGTH       nominal read length (for indexing) [default: 100]
    -p FILE         pulldown reads for variant hits into the named zip-file.
    -s              single stranded k-mer extraction
    -v              produce verbose output

Format Strings

Format strings control the fields that get included in the output.
The value is a comma separated list. The order of fields is fixed,
so the order in the format string is irrelevant.  The supported
fields are:

    ks              Kolmogorov-Smirnov statistics
"""

import array
import math
import random
import sys
import zipfile

import docopt
import tqdm
import yaml

from pykmer.basics import ham, kmer, kmers, kmersList, kmersWithPos, kmersWithPosList, kmersWithPosLists, lcp, murmer, rc, render
from pykmer.file import openFile, readFasta, readFastq
from pykmer.misc import unionfind
from pykmer.sparse import sparse
from pykmer.stats import counts2cdf, pdf2cdf, logBinEq, ksDistance2, logGammaP, logGammaQ, logLowerGamma
from zotmer.library.align import glocalAlignment, revComp
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

        #xh = murmer(x, 17)
        #yh = murmer(y, 17)
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

class PositionalKmerIndex(object):
    def __init__(self, K):
        self.K = K
        self.N = 0
        self.idx = {}

    def addSeq(self, seq0, seq1 = None):
        n = self.N
        self.N += 1
        for (x,p) in kmersWithPosList(self.K, seq0, False):
            p -= 1
            if x not in self.idx:
                self.idx[x] = []
            self.idx[x].append((n, p))
        if seq1 is not None:
            for (x,p) in kmersWithPosList(self.K, seq1, False):
                p -= 1
                if x not in self.idx:
                    self.idx[x] = []
                self.idx[x].append((n, p))
        return n

    def find(self, seq):
        hits = {}
        (xs, ys) = kmersWithPosLists(self.K, seq)
        for (zs, s) in [(xs, True), (ys, False)]:
            for (x,p) in zs:
                if x not in self.idx:
                    continue
                p -= 1
                for (n,q) in self.idx[x]:
                    qp = (q - p, s)
                    if n not in hits:
                        hits[n] = {}
                    if qp not in hits[n]:
                        hits[n][qp] = 0
                    hits[n][qp] += 1

        cand = {}
        for (n, ps) in hits.items():
            for ((p,s),c) in ps.items():
                if n not in cand:
                    cand[n] = (c, [(p, s)])
                elif c > cand[n][0]:
                    cand[n] = (c, [(p, s)])
                elif c == cand[n][0]:
                    cand[n][1].append((p, s))

        res = []
        for n in cand.keys():
            (c, pss) = cand[n]
            for (p0, s) in pss:
                if s:
                    zs = [posKmer(self.K, x, p0 + p - 1) for (x,p) in xs]
                else:
                    zs = [posKmer(self.K, y, p0 + p - 1) for (y,p) in ys]
                res.append((n, zs))
        return res

def nearestPosKmers(K, seq, xps):
    idx = {}
    res = {}
    for (x,p) in kmersWithPosList(K, seq, False):
        p -= 1
        assert p not in idx
        idx[p] = x
        res[p] = [K+1, 0]

    for xp in xps.iterkeys():
        (x,p) = unpackPosKmer(K, xp)
        if p not in idx:
            continue
        y = idx[p]
        d = ham(x, y)
        if d > res[p][0]:
            continue
        if d < res[p][0]:
            res[p] = [d, xps[xp]]
        else:
            res[p][1] = max(res[p][1], xps[xp])
    return sorted([(p, d, c) for (p, (d,c)) in res.items()])

def posKmersToKmers(K, xps):
    xs = {}
    for xp in xps.iterkeys():
        (x,p) = unpackPosKmer(K, xp)
        c = xps[xp]
        xs[x] = c + xs.get(x, 0)
    return xs

def renderPath(K, xs):
    if len(xs) == 0:
        return ''
    r = [render(K, xs[0])]
    for x in xs[1:]:
        r.append("ACGT"[x&3])
    return ''.join(r)

def subtractFlanks(K, flankKmers, xps):
    xs = {}
    for xp in xps.iterkeys():
        (x,p) = unpackPosKmer(K, xp)
        if x in flankKmers:
            continue
        c = xps[xp]
        xs[x] = c + xs.get(x, 0)
    return xs

def makeIndexedVariant(v, K, L):
    (wtSeq, mutSeq) = v.context(K)

    r = {}
    r['hgvs'] = str(v)
    r['wtSeq'] = wtSeq
    r['mutSeq'] = mutSeq

    if v.anonymous():
        # An anonymous sequence, so extract the flanking k-mers
        rng = v.range()
        big = v.loadAccession()
        lhsSt = rng[0] - 2*K
        lhsEn = rng[0]
        lhs = big[lhsSt:lhsEn]
        rhsSt = rng[1]
        rhsEn = rng[1] + 2*K
        rhs = big[rhsSt:rhsEn]
        r['mutLhs'] = lhs
        r['mutRhs'] = rhs
        r['mutSeq'] = None
        return r

    if v.cryptic() is not None:
        # A cryptic sequence insertion, so extract the flanking k-mers
        rng = v.range()
        big = v.loadAccession()
        lhsSt = rng[0] - 1 - 2*K
        lhsEn = rng[0] - 1
        lhs = big[lhsSt:lhsEn]
        rhsSt = rng[1] - 1
        rhsEn = rng[1] - 1 + 2*K
        rhs = big[rhsSt:rhsEn]
        r['mutLhs'] = lhs
        r['mutRhs'] = rhs
        return r

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
        return None
    return r

def addProbesToIdx(n, xs, idx):
    for x in xs:
        if x not in idx:
            idx[x] = set([])
        idx[x].add(n)

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

    L = int(opts['-L'])

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
                r = makeIndexedVariant(v, K, L)
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

    wtIdx = {}
    mutIdx = {}
    mutNegIdx = {}

    wtPosIdx = PositionalKmerIndex(K)
    mutPosIdx = PositionalKmerIndex(K)

    for n in range(V):
        itm = hgvsVars[n]
        h = itm['hgvs']
        v = makeHGVS(h)
        itm['var'] = v
        wtProbes = set(kmersWithPosList(K, itm['wtSeq'], False))
        if v.anonymous() or v.cryptic():
            mutProbes = set(kmers(K, itm['mutLhs'], False)) | set(kmers(K, itm['mutRhs'], False))
            mutNegProbes = wtProbes - mutProbes
        else:
            mutProbes = set(kmersWithPosList(K, itm['mutSeq'], False))
            com = wtProbes - mutProbes
            wtProbes -= com
            mutProbes -= com
            mutNegProbes = set([])
        itm['wtProbes'] = wtProbes
        itm['mutProbes'] = mutProbes
        itm['mutNegProbes'] = mutNegProbes

        #for x in mutProbes:
        #    print n, render(K, x)
        addProbesToIdx(n, wtProbes, wtIdx)
        addProbesToIdx(n, mutProbes, mutIdx)
        addProbesToIdx(n, mutNegProbes, mutNegIdx)

        wtPosIdx.addSeq(itm['wtSeq'])
        if itm['mutSeq'] is not None:
            mutPosIdx.addSeq(itm['mutSeq'])
        else:
            mutPosIdx.addSeq(itm['mutLhs'], itm['mutRhs'])

    if verbose:
        print >> sys.stderr, "done."

    combineStrands = True
    if opts['-s']:
        combineStrands = False

    pulldown = None
    if opts['-p']:
        pulldown = {}

    wtRes = [{} for n in range(V)]
    mutRes = [{} for n in range(V)]

    rn = 0
    #for itm in reads(opts['<input>'], K=K, paired=False, reads=True, kmers=True, both=True, verbose=verbose):
    for itm in reads(opts['<input>'], K=K, paired=False, reads=True, kmers=False, both=True, verbose=verbose):
        rn += 1
        #fq = itm[0][0]         # Read, End-0
        #xs = itm[1][0][0]      # Kmers, End-0, Strand-0
        fq = itm[0]      # Read, End-0

        for (n, pxs) in wtPosIdx.find(fq[1]):
            for px in pxs:
                wtRes[n][px] = 1 + wtRes[n].get(px, 0)

        for (n, pxs) in mutPosIdx.find(fq[1]):
            v = hgvsVars[n]['var']
            for px in pxs:
                mutRes[n][px] = 1 + mutRes[n].get(px, 0)
            if pulldown is not None:
                if n not in pulldown:
                    pulldown[n] = []
                pulldown[n] += fq

    Q = 10
    nullCdf = pdf2cdf(nullModelPdf(K))

    hdrShown = False
    for n in range(V):
        v = hgvsVars[n]['var']
        h = hgvsVars[n]['hgvs']

        wtSeq = hgvsVars[n]['wtSeq']
        wtHits = nearestPosKmers(K, wtSeq, wtRes[n])
        wtMin = 0
        if len(wtHits) > 0:
            wtMin = min([c for (p, d, c) in wtHits])
        wtHs = [0 for j in range(K+1)]
        for (p, d, c) in wtHits:
            if d <= K:
                wtHs[d] += 1
        wtCdf = counts2cdf(wtHs)
        wtD = ksDistance2(wtCdf, nullCdf)[0]

        if v.cryptic() is not None:
            # Decide if it looks like there as an ALU
            # between the endpoints
            lhs = hgvsVars[n]['mutLhs']
            rhs = hgvsVars[n]['mutRhs']
            alu = kmersList(K, AluY, True)
            flanks = set(kmersList(K, lhs, True)) | set(kmersList(K, rhs, True)) | set(kmersList(K, wtSeq, True))
            mx = subtractFlanks(K, flanks, mutRes[n])
            for (x,c) in mx.iteritems():
                j = 0
                zs = []
                for y in alu:
                    j0 = lcp(K, x, y)
                    if j0 < j:
                        continue
                    if j0 > j:
                        j = j0
                        zs = []
                    zs.append(y)
                z = random.choice(zs)
                print '%d\t%s\t%d\t%d\t%s' % (n, render(K, x), c, j, render(K, z))
            mutHs = [0 for j in range(K+1)]
            mutMin = None
            for x in kmersList(K, AluY, False):
                # Hamming
                #d = K+1
                #z = None
                #for y in mx.iterkeys():
                #    d0 = ham(x, y)
                #    if d0 < d:
                #        d = d0
                #        z = y
                # Longest Prefix
                d = 0
                z = None
                for y in mx.iterkeys():
                    d0 = lcp(K, x, y)
                    if d0 > d:
                        d = d0
                        z = y
                if z is not None:
                    mutHs[d] += mx[z]
                    if mutMin is None:
                        mutMin = mx[z]
                    else:
                        mutMin = min(mutMin, mx[z])
            if mutMin is None:
                mutMin = 0
            mutCdf = counts2cdf(mutHs)
            mutD = ksDistance2(mutCdf, nullCdf)[0]
        else:
            if v.anonymous():
                # Look for a path linking the ends
                mx = posKmersToKmers(K, mutRes[n])

                lhs = hgvsVars[n]['mutLhs']
                # Scan from the right hand end of the lhs
                # to find the highest coverage extant k-mer
                lx = None
                lxc = 0
                lxp = 0
                p0 = len(lhs) - K
                for x in kmersList(K, lhs):
                    # resolve equals by preferring the right-most
                    if x in mx and mx[x] >= lxc:
                        lx = x
                        lxc = mx[x]
                        lxp = p0
                    p0 -= 1

                rhs = hgvsVars[n]['mutRhs']
                # Scan from the left hand end of the rhs
                # to find the highest coverage extant k-mer
                rx = None
                rxc = 0
                rxp = 0
                p0 = 0
                for x in kmersList(K, rhs):
                    # resolve equals by preferring the left-most
                    if x in mx and mx[x] > rxc:
                        rx = x
                        rxc = mx[x]
                        rxp = p0
                    p0 += 1

                mutPth = None
                if lx is not None and rx is not None:
                    print lxp, render(K, lx), rxp, render(K, rx), v.size()
                    ll = v.size() + lxp + rxp + 1 + K
                    mutPth = interpolate(K, mx, lx, rx, ll)

                if mutPth is None:
                    mutSeq = ''
                    mutPth = []
                    mutMin = 0
                else:
                    mutSeq = renderPath(K, mutPth)
                    mutSeq = mutSeq[lxp+K:-(rxp+K)]
                    assert len(mutSeq) == v.size()
                    mutMin = min([mx[x] for x in mutPth[1:-1]])
                mutHs = [0 for j in range(K+1)]
                for x in mutPth:
                    mutHs[0] += 1
            else:
                mutSeq = hgvsVars[n]['mutSeq']
                mutHits = nearestPosKmers(K, mutSeq, mutRes[n])
                mutMin = 0
                if len(mutHits) > 0:
                    mutMin = min([c for (p, d, c) in mutHits])
                mutHs = [0 for j in range(K+1)]
                for (p, d, c) in mutHits:
                    if d <= K:
                        mutHs[d] += 1

            mutCdf = counts2cdf(mutHs)
            mutD = ksDistance2(mutCdf, nullCdf)[0]

        hdrs = ['n']
        fmts = ['%d']
        outs = [n]

        wtAllele = ((wtMin > Q) and (wtD > 0.8))
        mutAllele = ((mutMin > Q) and (mutD > 0.8))
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

        if 'pos' in fmt:
            wtLen = len(wtSeq) - K + 1
            mutLen = len(mutSeq) - K + 1
            seqLen = max(wtLen, mutLen)
            i = 0
            j = 0
            for p in range(seqLen):
                wd = 0
                wc = 0
                if i < len(wtHits) and wtHits[i][0] == p:
                    wd = wtHits[i][1]
                    if wd == K+1:
                        wd = 0
                    wc = wtHits[i][2]
                    i += 1
                md = 0
                mc = 0
                if j < len(mutHits) and mutHits[j][0] == p:
                    md = mutHits[j][1]
                    if md == K+1:
                        md = 0
                    mc = mutHits[j][2]
                    j += 1

                hdrs1 = hdrs + ['pos', 'wtHam', 'wtCnt', 'mutHam', 'mutCnt']
                fmts1 = fmts + ['%d', '%d', '%d', '%d', '%d']
                outs1 = outs + [p, wd, wc, md, mc]

                hdrs1 += ['hgvs']
                fmts1 += ['%s']
                outs1 += [h]

                if not hdrShown:
                    hdrShown = True
                    print '\t'.join(hdrs1)
                print '\t'.join(fmts1) % tuple(outs1)
            continue

        if 'ham' in fmt:
            for j in range(K+1):
                hdrs1 = hdrs + ['ham', 'wtCnt', 'wtCum', 'mutCnt', 'mutCum', 'nullCum']
                fmts1 = fmts + ['%d', '%d', '%g', '%d', '%g', '%g']
                outs1 = outs + [j, wtHs[j], wtCdf[j], mutHs[j], mutCdf[j], nullCdf[j]]

                hdrs1 += ['hgvs']
                fmts1 += ['%s']
                outs1 += [h]

                if not hdrShown:
                    hdrShown = True
                    print '\t'.join(hdrs1)
                print '\t'.join(fmts1) % tuple(outs1)
            continue

        hdrs += ['hgvs']
        fmts += ['%s']
        outs += [h]

        if not hdrShown:
            hdrShown = True
            print '\t'.join(hdrs)
        print '\t'.join(fmts) % tuple(outs)

    if pulldown is not None:
        with zipfile.ZipFile(opts['-p'], 'w', zipfile.ZIP_DEFLATED) as zf:
            for n in sorted(pulldown.keys()):
                h = hgvsVars[n]['hgvs']
                rds = pulldown[n]
                zf.writestr(h + '.fastq', '\n'.join(rds))

if __name__ == '__main__':
    main(sys.argv[1:])
