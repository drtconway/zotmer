from base import Cmd

from pykmer.basics import fasta, kmersWithPos, lcp, rc, render
from pykmer.file import readFasta
import pykmer.kset as kset
import pykmer.kfset as kfset
from pykmer.stats import logChi2Crit, logChi2CDF, counts2cdf, ksDistance2, logLowerGamma, logGamma
from pykmer.uf import uf
from pykmer.exceptions import MismatchedK

import math
import sys

def getK(ins):
    k = None
    for fn in ins:
        k0 = kset.probeK(fn)
        if k is None:
            k = k0
        elif k != k0:
            raise MismatchedK(k, k0)
    return k

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

def rank1(xs, x):
    l = 0
    h = len(xs)
    while h > l:
        m = (h + l) // 2
        y = xs[m]
        if y == x:
            return m + 1
        if y < x:
            l = m + 1
        else:
            h = m - 1
    return l

def rank2(xs, x0, x1):
    assert x0 < x1
    r0 = rank1(xs, x0)
    r1 = r0
    while r1 < len(xs) and xs[r1] < x1:
        r1 += 1
    return (r0, r1)

def succ(K, xs, x):
    m = (1 << (2*K)) - 1
    y0 = (x << 2) & m
    y1 = y0 + 4
    r = rank2(xs, y0, y1)
    return [xs[i] for i in range(r[0], r[1])]

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
        K = getK(opts['<input>'])
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
            (meta, xs) = kset.read(fn)
            sat = meta['acgt'][0] + meta['acgt'][3]
            px = 0
            pxn = 0
            for x in xs:
                xp = (x >> s)
                if px != xp:
                    print >> sys.stderr, '%s\t%d' % (render(J, px), pxn)
                    px = xp
                    pxn = 0
                pxn += 1
                resolveLcp(K, idx, lcps, x)
                m += 1

            g = qat * sat - 0.5 * (qat + sat) + 0.5
            gs = sat * sat - 0.5 * (sat + sat) + 0.5
            print >> sys.stderr, "qat =", qat
            print >> sys.stderr, "sat =", sat
            print >> sys.stderr, "g =", g
            print >> sys.stderr, "gs =", gs
            nm = [null(g, m, j) for j in range(0, K+1)]

            cnt = [[0 for j in range(K+1)] for i in range(len(nms))]
            alleles = [[] for i in range(len(nms))]
            for i in range(len(idx)):
                cnt[idx[i][1]][lcps[i][0]] += 1
                p = 1 - null(g, m, lcps[i][0])
                alleles[idx[i][1]].append((idx[i][2], p, lcps[i][0], lcps[i][1]))

            for i in range(len(nms)):
                dst = counts2cdf(cnt[i])
                (ks1, ks2) = ksDistance2(nm, dst)
                kspv = -2*(ks1**2)*lens[i] / math.log(10)
                alleles[i].sort()

                js = [None for j in range(lens[i] - K + 1)]
                xs = [None for j in range(lens[i] - K + 1)]
                qs = [None for j in range(lens[i])]
                ax = [[1, 1, 1, 1] for j in range(lens[i])]
                dx = [0 for j in range(lens[i])]
                for (p, v, k, x) in alleles[i]:
                    if p < 0:
                        x = rc(K, x)
                        p = -p
                    p -= 1 # positions from 1-offset -> 0-offset
                    if js[p] < k:
                        js[p] = k
                        xs[p] = x
                    for l in range(K):
                        j = p + l
                        qs[j] = max(qs[j], k)
                        dx[j] += 1
                        k = K - l - 1
                        b = (x >> (2*k)) & 3
                        ax[j][b] *= v

                pcrit = 0.01 / (len(xs) - 1)
                ocrit = 0
                for j in range(K, 1, -1):
                    if 1 - null(gs, m, j) > pcrit:
                        break
                    ocrit = j

                suf = indexSuffixes(K, xs)
                u = uf()
                for j in range(len(xs) - 1):
                    (o, k) = maxOverlap(K, xs[j], j, suf)
                    if o >= ocrit:
                        u.union(j, k)

                idx = {}
                for j in range(len(xs)):
                    j0 = u.find(j)
                    if j0 not in idx:
                        idx[j0] = []
                    idx[j0].append(j)

                for (_, ls) in idx.items():
                    a = [0 for j in range(lens[i])]
                    q = []
                    for j in ls:
                        x = xs[j]
                        q0 = null(g, m, js[j])
                        q.append(q0)
                        for l in range(K):
                            k = (K - l - 1)
                            b = (x >> (2*k)) & 3
                            a[j+l] |= 1 << b

                    hasZeros = False
                    for j in range(len(a)):
                        if a[j] == 0:
                            hasZeros = True
                        a[j] = fasta(a[j])

                    if hasZeros:
                        continue
                    else:
                        chi2 = -2*sum([math.log(q0) for q0 in q])
                        print '%s\t%d\t%g\t%s' % (nms[i], lens[i], chi2, ''.join(a))

                #a = []
                #for j in range(len(ax)):
                #    b = 0
                #    for k in range(4):
                #        if ax[j][k] < 1e-12:
                #            b |= 1 << k
                #    a.append(fasta(b))
                #    #print j, dx[j], fasta(b), ax[j]
#
#                chi2 = 0
#                for j in qs:
#                    q0 = null(g, m, j)
#                    chi2 += math.log(q0 + 1e-100)
#                for j in range(len(xs) - 1):
#                    (k, _) = overlap(K, xs, j)
#                    q0 = null(gs, K - k, k)
#                    chi2 += math.log(q0 + 1e-100)
#
#                chi2 *= -2
#                pv = logChi2LowerPval(2*(len(qs) + len(xs) - 1), chi2)/math.log(10)
#
#                print '%d\t%d\t%d\t%g\t%g\t%g\t%g\t%s\t%s' % (i, lens[i], m, ks1, kspv, chi2, pv, nms[i], ''.join(a))
                sys.stdout.flush()

def add(cmds):
    cmds['scan'] = Scan()
