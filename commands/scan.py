from base import Cmd

from pykmer.basics import fasta, kmersWithPos, lcp, render
from pykmer.file import readFasta
import pykmer.kset as kset
import pykmer.kfset as kfset
from pykmer.stats import logChi2Crit, logChi2CDF, counts2cdf, ksDistance2
from pykmer.exceptions import MismatchedK

import math
import sys

def getK(ins):
    k = None
    for fn in ins:
        k0 = kfset.probeK(fn)
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

def gcat(k, x):
    cs = [0, 0, 0, 0]
    for i in range(k):
        cs[x&3] += 1
        x >>= 2
    return (float(cs[0] + cs[3])/k, float(cs[1] + cs[2])/k)

def null(g, s, j):
    if j == 0:
        return 0.0
    return math.exp(s*math.log1p(-math.pow(g,j)))

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
                for (x,p) in kmersWithPos(K, seq):
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
            sat = meta['acgt'][0] + meta['acgt'][3]
            px = 0
            pxn = 0
            for (x,_) in xs:
                xp = (x >> s)
                if px != xp:
                    print >> sys.stderr, '%s\t%d' % (render(J, px), pxn)
                    px = xp
                    pxn = 0
                pxn += 1
                resolveLcp(K, idx, lcps, x)
                m += 1

            g = qat * sat - 0.5 * (qat + sat) + 0.5
            cnt = [[0 for j in range(K+1)] for i in range(len(nms))]
            res = [0 for i in range(len(nms))]
            alleles = [[] for i in range(len(nms))]
            nm = [null(g, m, j) for j in range(0, K+1)]
            for i in range(len(idx)):
                p = math.log1p(-null(g, m, lcps[i][0]))
                res[idx[i][1]] += p
                cnt[idx[i][1]][lcps[i][0]] += 1
                alleles[idx[i][1]].append((idx[i][2], p, lcps[i][0], lcps[i][1]))
            for i in range(len(res)):
                dst = counts2cdf(cnt[i])
                (ks1, ks2) = ksDistance2(nm, dst)
                a = [[0, 0, 0, 0] for j in range(lens[i])]
                alleles[i].sort()
                for (p, v, k, x) in alleles[i]:
                    for k in range(K):
                        j = K - k - 1
                        b = (x >> (2*j)) & 3
                        q = p + k
                        a[q][b] += v
                ax = []
                crit = logChi2Crit(K//2, math.log(0.01))
                print "crit = ", crit
                for acgt in a:
                    print >> sys.stderr, b, v
                    b = 0
                    for j in range(len(acgt)):
                        if -2*acgt[j] > crit:
                            b |= 1 << j
                    ax.append(fasta(b))
                a = ''.join(ax)
                print '%d\t%d\t%d\t%g\t%g\t%g\t%s\t%s' % (i, lens[i], m, ks1, -2*res[i], min(0,logChi2CDF(2*lens[i], -2*res[i]))/math.log(10), nms[i], a)

def add(cmds):
    cmds['scan'] = Scan()
