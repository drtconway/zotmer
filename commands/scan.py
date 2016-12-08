from base import Cmd

from pykmer.basics import fasta, kmersWithPos, ham, lcp, rc, render
from pykmer.file import readFasta
import pykmer.kset as kset
import pykmer.kfset as kfset
from pykmer.stats import counts2cdf, ksDistance2, log1mexp
from pykmer.uf import uf
from pykmer.sparse import sparse
from pykmer.exceptions import MismatchedK

import array
import math
import sys

def findLcp(idx, x):
    l = 0
    h = len(idx)
    while h > l:
        m = (h + l) // 2
        y = idx[m][0]
        if y == x:
            return m
        if y < x:
            l = m + 1
        else:
            h = m - 1
    return l

def resolveLcp(k, idx, lcps, x):
    i = findLcp(idx, x)
    z = min(i+1, len(idx))
    for h in range(z):
        j = z - h - 1
        l = lcp(k, x, idx[j][0])
        if l > lcps[j][0]:
            lcps[j] = (l, x)
        else:
            break
    for j in range(i + 1, len(idx)):
        l = lcp(k, x, idx[j][0])
        if l > lcps[j][0]:
            lcps[j] = (l, x)
        else:
            break

def succ(K, xs, x):
    m = (1 << (2*K)) - 1
    y0 = (x << 2) & m
    y1 = y0 + 4
    r0 = xs.rank(y0)
    r1 = xs.rank(y1)
    return [xs.select(i) for i in xrange(r0, r1)]

def los(K, x, y):
    m = (1 << (2*K)) - 1
    for j in range(1, K):
        m >>= 2
        if (x & m) == (y >> (2*j)):
            return K - j
    return 0

def indexSuffixes(K, xs):
    idx = [{} for j in range(K+1)]
    for i in range(len(xs)):
        x = xs[i]
        m = (1 << (2*K)) - 1
        for j in range(1, K):
            x >>= 2
            l = K - j
            if x not in idx[l]:
                idx[l][x] = []
            idx[l][x].append(i)
    return idx

def maxOverlap(K, x, i, idx):
    m = (1 << (2*K)) - 1
    for j in range(1, K):
        m >>= 2
        l = K - j
        y = x & m
        if y in idx[l]:
            for k in idx[l][y]:
                if k > i:
                    return (l, k)
    return (0, i)

def overlap(K, xs, i):
    m = (1 << (2*K)) - 1
    for j in range(1, K):
        if i + j >= len(xs):
            return (0, i)
        m >>= 2
        if (xs[i] & m) == (xs[i + j] >> (2*j)):
            return (K - j, i + j)
    return (0, i)

def gcat(k, x):
    cs = [0, 0, 0, 0]
    for i in range(k):
        cs[x&3] += 1
        x >>= 2
    return (float(cs[0] + cs[3])/k, float(cs[1] + cs[2])/k)

def logChi2LowerPval(n, x):
    m = n * (1 - 2.0/(9*n))**3
    if x < m:
        a = logLowerGamma(n/2.0, x/2.0)
        b = logGamma(n/2.0)
        c = a - b
        if c > 0:
            c = 0
        #print n, x, a, b, c
        return c
    else:
        return math.log1p(-math.exp(logChi2CDF(n, x)))

def null(g, s, j):
    if j == 0:
        return 0.0
    return math.exp(s*math.log1p(-math.pow(g,j)))

def logNull(g, s, j):
    return s*math.log1p(-math.pow(g,j))

class Scan(Cmd):
    def run(self, opts):
        gfn = opts['<genes>']
        K = kfset.read(opts['<input>'][0])[0]['K']
        idx = []
        lens = []
        nms = []
        n = 0
        qacgt = [0, 0, 0, 0]
        with open(gfn) as f:
            for (nm, seq) in readFasta(f):
                if len(seq) < K:
                    print >> sys.stderr, "warning: `%s' skipped" % (nm,)
                    continue
                nms.append(nm)
                lens.append(len(seq))
                for (x,p) in kmersWithPos(K, seq, True):
                    idx.append((x, n, p))
                    qacgt[x&3] += 1
                n += 1
        idx.sort()
        qat = float(qacgt[1]+qacgt[2])/float(qacgt[0]+qacgt[1]+qacgt[2]+qacgt[3])
        J = 3
        s = 2*(K - J)
        for fn in opts['<input>']:
            m = 0
            lcps = [(0, 0) for i in range(len(idx))]
            (meta, xs) = kfset.read(fn)
            sacgt = [0, 0, 0, 0]
            px = 0
            pxn = 0
            X = array.array('L', [])
            for (x,_) in xs:
                xp = (x >> s)
                if px != xp:
                    print >> sys.stderr, '%s\t%d' % (render(J, px), pxn)
                    px = xp
                    pxn = 0
                sacgt[x&3] += 1
                pxn += 1
                resolveLcp(K, idx, lcps, x)
                X.append(x)
                m += 1
            X = sparse(1 << (2*K), X)

            sat = float(sacgt[1]+sacgt[2])/float(sacgt[0]+sacgt[1]+sacgt[2]+sacgt[3])
            g = qat * sat - 0.5 * (qat + sat) + 0.5
            print >> sys.stderr, "g =", g
            nm = [null(g, m, j) for j in range(0, K+1)]

            cnt = [[0 for j in xrange(K+1)] for i in xrange(len(nms))]
            alleles = [[] for i in xrange(len(nms))]
            for i in range(len(idx)):
                (x, n, p) = idx[i]
                (l, y) = lcps[i]
                cnt[n][l] += 1
                alleles[n].append((p, x, y, l))

            for i in xrange(len(nms)):
                qc = math.log(0.05/(lens[i] - K + 1)/2)
                a = [[0, 0, 0] for j in xrange(lens[i] - K + 1)]
                for (p, x, y, l) in alleles[i]:
                    if p < 0:
                        q = -(p + 1)
                        x = rc(K, x)
                        y = rc(K, y)
                    else:
                        q = p - 1
                    if l > a[q][2]:
                        a[q] = [x, y, l]
                m = (1 << (2*K - 2)) - 1
                py = 0
                u = uf()
                for j in xrange(lens[i] - K + 1):
                    x = a[j][1]
                    y = x >> 2
                    if j > 0:
                        d = ham(py, y)
                        if d == 0:
                            u.union(j-1, j)
                    py = x & m
                udx = {}
                for j in xrange(lens[i] - K + 1):
                    v = u.find(j)
                    if v not in udx:
                        udx[v] = []
                    udx[v].append(j)
                fragments = {}
                fdxLhs = {}
                fdxRhs = {}
                kx = []
                for (jx, js) in udx.iteritems():
                    q = 0
                    for j in js:
                        q += math.log1p(-nm[a[j][2]])
                    if q > math.log(0.05/len(js)):
                        continue
                    print i, jx, len(js), q
                    fragments[jx] = (js[0], a[js[0]][1], js[-1] + 1, a[js[-1]][1])
                    kx.append((len(js), jx))
                    fdxLhs[a[js[0]][1]] = jx
                    fdxRhs[a[js[-1]][1]] = jx
                kx.sort()
                links = {}
                for (_, jx) in kx[::-1]:
                    (jL, xL, jR, xR) = fragments[jx]
                    if jR == lens[i] - K + 1:
                        continue
                    x = xR
                    xs = []
                    lnk = None
                    for k in xrange(100):
                        ys = succ(K, X, x)
                        if len(ys) != 1:
                            break
                        x = ys[0]
                        if x in fdxLhs:
                            lnk = fdxLhs[x]
                            break
                        xs.append(x)
                    if lnk is not None:
                        links[jx] = xs
                        u.union(jx, lnk)
                udx = {}
                for j in [jx for (_, jx) in kx]:
                    v = u.find(j)
                    if v not in udx:
                        udx[v] = []
                    udx[v].append(j)
                res = []
                for (jxx, jxs) in udx.iteritems():
                    fs = [(fragments[jx], jx) for jx in jxs]
                    fs.sort()
                    sxs = []
                    for fj in xrange(len(fs)):
                        (f, jx) = fs[fj]
                        beg = f[0]
                        end = f[2]
                        if fj == 0:
                            for j in xrange(beg):
                                sxs.append((0, 0))
                        xs = links.get(jx, None)
                        for j in xrange(beg, end):
                            x = a[j][1]
                            l = a[j][2]
                            sxs.append((x, l))
                        if xs:
                            for x in xs:
                                sxs.append((x, 27))
                        else:
                            if fj < len(fs) - 1:
                                ff = fs[fj+1][0]
                                nxt = ff[0]
                            else:
                                nxt = lens[i] - K + 1
                            for j in xrange(end, nxt):
                                sxs.append((0, 0))
                    seq = [[0, 0, 0, 0] for j in xrange(len(sxs) + K - 1)]
                    for j in xrange(len(sxs)):
                        (x, l) = sxs[j]
                        p = math.log1p(-nm[l])
                        for k in xrange(K):
                            seq[j + K - k - 1][x&3] += p
                            x >>= 2
                    ax = []
                    p = 0
                    for j in xrange(len(seq)):
                        b = 0
                        for k in xrange(4):
                            if seq[j][k] < qc:
                                b |= 1 << k
                        ax.append(fasta(b))
                        p += log1mexp(min(-1e-100, sum(seq[j])))
                    dst = counts2cdf(cnt[i])
                    (_, kd) = ksDistance2(dst, nm)
                    q = log1mexp(p)
                    res.append((q, kd, ''.join(ax)))
                res.sort()
                for j in xrange(len(res)):
                    if j > 0 and res[j][0] > qc:
                        break
                    print '%d\t%d\t%d\t%g\t%g\t%s\t%s' % (i, lens[i], len(res[j][2]), res[j][1], res[j][0], nms[i], res[j][2])
                    if res[j][0] > qc:
                        break
                sys.stdout.flush()

def add(cmds):
    cmds['scan'] = Scan()
