"""
zot hgvs-find - search for known variants in read data.

Usage:
    zot hgvs-find -R [options] <variant>
    zot hgvs-find -X [options] <index> [<variant>...]
    zot hgvs-find [options] <index> <input>...

Options:
    -R              report the ambiguity associated with a variant
    -X              index HGVS variants
    -A ALPHA        alpha level for Kolmogorov-Smirnov test [default: 0.001]
    -B BINFUNC      binning function to use [default: none]
    -S SCALE        scaling factor for binning function [default: 1.0]
    -f FILENAME     read variants from a file
    -F FORMAT       a format string for printing results [default: ks]
    -k K            value of k to use [default: 25]
    -g PATH         directory of FASTQ reference sequences
    -p FILE         pulldown reads for variant hits into the named zip-file.
    -s              single stranded k-mer extraction
    -v              produce verbose output

Format Strings

Format strings control the fields that get included in the output.
The value is a comma separated list. The order of fields is fixed,
so the order in the format string is irrelevant.  The supported
fields are:

    chi2            Chi2 statistics
    kl              Kullback-Leibler statistics
    ks              Kolmogorov-Smirnov statistics
    hist            Print long format results with the histograms
    glo             Include global frequency counts in long format results

"""

import array
import math
import sys
import zipfile

import docopt
import progressbar
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

def neigh(K, x):
    r = []
    for j in xrange(K):
        for y in [1,2,3]:
            z = x ^ (y << (2*j))
            r.append(z)
    return r

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

    if False:
        print '>'+v['hgvs']
        print s[wt[0]-100:wt[1]+100]
    if False:
        print '>'+v['hgvs']
        print r[mut[0]-100:mut[1]+100]

    wtXs = set(kmers(k, s[wt[0]:wt[1]], False))
    mutXs = set(kmers(k, r[mut[0]:mut[1]], False))
    com = wtXs & mutXs
    wtXs -= com
    mutXs -= com

    wtXs = list(wtXs)
    wtXs.sort()

    mutXs = list(mutXs)
    mutXs.sort()

    #print '%s\t%s\t%d\t%d\t%d' % (v['hgvs'], v['type'], len(wtXs), len(mutXs), len(com))
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

    fmt = set(opts['-F'].split(','))

    if opts['-v']:
        print >> sys.stderr, "loading index."

    with open(opts['<index>']) as f:
        itms = yaml.load(f)

    wtProbes = []
    wtRes = []
    mutProbes = []
    mutRes = []
    hs = []
    wtIdx = {}
    mutIdx = {}
    for itm in itms:
        n = len(hs)
        h = itm['hgvs']
        hs.append(h)
        wtProbes.append(set([]))
        wtRes.append({})
        for x in itm['wt']:
            y = kmer(x)
            if y not in wtIdx:
                wtIdx[y] = set([])
            wtIdx[y].add(n)
            wtProbes[-1].add(y)
        mutProbes.append(set([]))
        mutRes.append({})
        for x in itm['mut']:
            y = kmer(x)
            if y not in mutIdx:
                mutIdx[y] = set([])
            mutIdx[y].add(n)
            mutProbes[-1].add(y)

    if opts['-v']:
        print >> sys.stderr, "done."

    combineStrands = True
    if opts['-s']:
        combineStrands = False

    pulldown = None
    if opts['-p']:
        pulldown = {}

    rn = 0
    M = (1 << 12) - 1

    bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
    bar.start()
    for fn in opts['<input>']:
        with openFile(fn) as f:
            for fq in readFastq(f):
                rn += 1
                if (rn & M) == 0:
                    bar.update(rn/1000.0)
                if (rn & M) == 0 and opts['-v']:
                    print >> sys.stderr, 'read pairs processed: %d' % (rn,)
                xs = kmersList(K, fq[1], True)

                wtHits = None
                for x in xs:
                    if x in wtIdx:
                        if wtHits is None:
                            wtHits = set([])
                        wtHits |= wtIdx[x]
                if wtHits is not None:
                    for n in wtHits:
                        wtXs = wtRes[n]
                        for x in xs:
                            wtXs[x] = 1 + wtXs.get(x, 0)
                mutHits = None
                for x in xs:
                    if x in mutIdx:
                        if mutHits is None:
                            mutHits = set([])
                        mutHits |= mutIdx[x]
                if mutHits is not None:
                    for n in mutHits:
                        mutXs = mutRes[n]
                        for x in xs:
                            mutXs[x] = 1 + mutXs.get(x, 0)
                        if pulldown is not None:
                            if n not in pulldown:
                                pulldown[n] = []
                            pulldown[n] += fq
    bar.finish()

    B = opts['-B']
    S = float(opts['-S'])

    wtMeans = []
    wtMedians = []
    mutMeans = []
    mutMedians = []
    for n in range(len(hs)):
        xs = [x for x in wtRes[n].values()]
        xh0 = hist(xs)
        xh = binHist(xh0, B, S)
        wtMeans.append((histMean(xh), n))
        wtMedians.append((histMedian(xh), n))

        ys = [y for y in mutRes[n].values()]
        yh0 = hist(ys)
        yh = binHist(yh0, B, S)
        mutMeans.append((histMean(yh), n))
        mutMedians.append((histMedian(yh), n))

    wtMeans.sort()
    wtMeanRanks = makeRanks(wtMeans)
    wtMedians.sort()
    wtMedianRanks = makeRanks(wtMedians)

    mutMeans.sort()
    mutMeanRanks = makeRanks(mutMeans)
    mutMedians.sort()
    mutMedianRanks = makeRanks(mutMedians)

    Q = 10

    zf = None
    if pulldown is not None:
        zf = zipfile.ZipFile(opts['-p'], 'w', zipfile.ZIP_DEFLATED)

    hdrShown = False
    for n in range(len(hs)):

        wtH = [0 for k in xrange(K+1)]
        for x in wtProbes[n]:
            dMax = K+1
            yMax = 0
            for y in wtRes[n].keys():
                d = ham(x, y)
                if d < dMax:
                    dMax = d
                    yMax = y
            if dMax <= K:
                print '%d\twt\t%s\t%d\t%d' % (n, render(K, x), dMax, wtRes[n][yMax])
                #wtH[dMax] += wtRes[n][yMax]
                wtH[dMax] += 1

        mutH = [0 for k in xrange(K+1)]
        for x in mutProbes[n]:
            dMax = K+1
            yMax = 0
            for y in mutRes[n].keys():
                d = ham(x, y)
                if d < dMax:
                    dMax = d
                    yMax = y
            if dMax <= K:
                print '%d\tmut\t%s\t%d\t%d' % (n, render(K, x), dMax, mutRes[n][yMax])
                #mutH[dMax] += mutRes[n][yMax]
                mutH[dMax] += 1

        continue
        wtHMean = float(sum([d*wtH[d] for d in xrange(K+1)]))/float(max(1,sum(wtH)))
        mutHMean = float(sum([d*mutH[d] for d in xrange(K+1)]))/float(max(1,sum(mutH)))

        #print '\t'.join(map(str, [n] + wtH))
        #print '\t'.join(map(str, [n] + mutH))

        h = hs[n]
        hdrs = ['n']
        fmts = ['%d']
        outs = [n]

        resV =  1*(wtMedianRanks[n][1] > Q) + 2*(mutMedianRanks[n][1] > Q)
        res = ['null', 'wt', 'mut', 'wt/mut'][resV]

        hdrs += ['res', 'wtMedian', 'wtMedianRank', 'mutMedian', 'mutMedianRank']
        fmts += ['%s', '%g', '%g', '%g', '%g']
        outs += [res, wtMedianRanks[n][1], wtMedianRanks[n][0], mutMedianRanks[n][1], mutMedianRanks[n][0]]

        resV =  1*(wtMeanRanks[n][1] > Q) + 2*(mutMeanRanks[n][1] > Q)
        res = ['null', 'wt', 'mut', 'wt/mut'][resV]

        hdrs += ['resMean', 'wtMean', 'wtMeanRank', 'mutMean', 'mutMeanRank']
        fmts += ['%s', '%g', '%g', '%g', '%g']
        outs += [res, wtMeanRanks[n][1], wtMeanRanks[n][0], mutMeanRanks[n][1], mutMeanRanks[n][0]]

        hdrs += ['wtHMean', 'mutHMean']
        fmts += ['%g', '%g']
        outs += [wtHMean, mutHMean]

        if 'hist' in fmt:
            xs = [x for x in wtRes[n].values()]
            xs.sort()
            nx = float(len(xs))
            nx = max(1.0, nx)
            xh0 = hist(xs)
            xh = binHist(xh0, B, S)

            ys = [y for y in mutRes[n].values()]
            ys.sort()
            ny = float(len(ys))
            ny = max(1.0, ny)
            yh0 = hist(ys)
            yh = binHist(yh0, B, S)

            wtT = 0
            wtP = 0.0
            mutT = 0
            mutP = 0.0
            cs = list(set(xh.keys() + yh.keys()))
            cs.sort()
            for c in cs:
                wtC = xh.get(c, 0)
                wtT += wtC
                wtP = wtT / nx
                mutC = yh.get(c, 0)
                mutT += mutC
                mutP = mutT / ny

                hdrs1 = hdrs + ['coverage', 'wtCnt', 'wtCum', 'mutCnt', 'mutCum']
                fmts1 = fmts + ['%d', '%d', '%g', '%d', '%g']
                outs1 = outs + [c, wtC, wtP, mutC, mutP]

                hdrs1 += ['hgvs']
                fmts1 += ['%s']
                outs1 += [h]

                if not hdrShown:
                    hdrShown = True
                    print '\t'.join(hdrs1)
                print '\t'.join(fmts1) % tuple(outs1)
        else:
            hdrs += ['hgvs']
            fmts += ['%s']
            outs += [h]

            if not hdrShown:
                hdrShown = True
                print '\t'.join(hdrs)
            print '\t'.join(fmts) % tuple(outs)

        if pulldown is not None and mutMedianRanks[n][1] > Q:
            rds = pulldown[n]
            zf.writestr(h + '.fastq', '\n'.join(rds))

    if zf is not None:
        zf.close()

if __name__ == '__main__':
    main(sys.argv[1:])
