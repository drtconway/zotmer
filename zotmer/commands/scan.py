"""
Usage:
    zot scan [-AX] <genes> <input>...

Options:
    -A          print all alleles
    -X          create an index
"""

from pykmer.basics import fasta, kmersWithPos, ham, lcp, rc, render
from pykmer.container import container
from pykmer.container.vectors import read64, write64, read32, write32, read32s, write32s
from pykmer.container.std import readKmers
from pykmer.file import readFasta
from pykmer.misc import unionfind
from pykmer.sparse import sparse
from pykmer.stats import counts2cdf, ksDistance2, log1mexp, logAdd, logChoose, logFac

import array
import cPickle
import docopt
import math
import sys

def lev(x, y):
    zx = len(x)
    zy = len(y)
    buf0 = array.array('I', [0 for i in xrange(zy + 1)])
    buf1 = array.array('I', [0 for i in xrange(zy + 1)])

    v0 = buf0
    v1 = buf1

    for i in xrange(zy):
        v0[i] = i

    for i in xrange(zx):
        v1[0] = i + 1
        for j in xrange(zy):
            if x[i] == y[j]:
                c = 0
            else:
                c = 1
            a1 = v0[j] + c
            a2 = v1[j] + 1
            a3 = v0[j + 1] + 1
            v1[j + 1] = min([a1, a2, a3])
        t = v1
        v1 = v0
        v0 = t
    return v0[zy]

def resolveLcp(K, S, L, Y, x):
    z = S.rank(x)
    for k in xrange(z):
        j = z - k - 1
        l = lcp(K, x, S.select(j))
        if l > L[j]:
            L[j] = l
            Y[j] = x
        else:
            break
    for j in xrange(z, S.count()):
        l = lcp(K, x, S.select(j))
        if l > L[j]:
            L[j] = l
            Y[j] = x
        else:
            break

def resolveAll(K, S, L, Y, xs):
    i = 0
    Z = S.count()
    for x in xs:
        while i < Z and S.select(i) <= x:
            y = S.select(i)
            l = lcp(K, x, y)
            if l > L[i]:
                L[i] = l
                Y[i] = x
            i += 1
    while i < Z:
        y = S.select(i)
        l = lcp(K, x, y)
        if l > L[i]:
            L[i] = l
            Y[i] = x
        i += 1

    for i in xrange(1, Z):
        x = Y[i-1]
        y = S.select(i)
        l = lcp(K, x, y)
        if l > L[i]:
            L[i] = l
            Y[i] = x

def succ(K, xs, x):
    m = (1 << (2*K)) - 1
    y0 = (x << 2) & m
    y1 = y0 + 4
    (r0, r1) = xs.rank2(y0, y1)
    return range(r0, r1)

def null(g, s, j):
    if j == 0:
        return 0.0
    return math.exp(s*math.log1p(-math.pow(g,j)))

def logNull(g, s, j):
    return s*math.log1p(-math.pow(g,j))

def logGammaPInt(n, x):
    if x < n:
        lx = math.log(x)
        w = 0
        v = logFac(n)
        s = w - v
        j = 1
        while True:
            w += lx
            v += math.log(j+n)
            u = logAdd(s, w - v)
            if s == u:
                break
            s = u
            j += 1
        return n*lx + s - x
    else:
        lx = math.log(x)
        s = 0
        v = 0
        for k in xrange(1, n):
            v += math.log(k)
            t = k*lx - v
            u = logAdd(t, s)
            s = u
        w = min(s - x, -1e-200)
        return log1mexp(w)

def chi2(n, x):
    d = int(math.ceil(n/2.0))
    c2 = x/2.0
    return logGammaPInt(d, c2)

def uniq(xs):
    px = None
    for x in xs:
        if x != px:
            yield x
            px = x

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    if opts['-X']:
        K = 27
        S = []
        N = 0
        qacgt = [0, 0, 0, 0]
        for fn in opts['<input>']:
            with open(fn) as f:
                for (nm, seq) in readFasta(f):
                    if len(seq) < K:
                        continue
                    for (x,p) in kmersWithPos(K, seq, True):
                        S.append(x) 
                        qacgt[x&3] += 1
                        N += 1
        S.sort()
        qacgt = [float(c)/float(N) for c in qacgt]
        S = sparse(2*K, array.array('L', uniq(S)))
        lens = array.array('I', [])
        nms = []
        seqs = []
        n = 0
        tmp = [[] for i in xrange(S.count())]
        for fn in opts['<input>']:
            with open(fn) as f:
                for (nm, seq) in readFasta(f):
                    if len(seq) < K:
                        print >> sys.stderr, "warning: `%s' skipped" % (nm,)
                        continue
                    nms.append(nm)
                    seqs.append(seq)
                    lens.append(len(seq))
                    for (x,p) in kmersWithPos(K, seq, True):
                        r = S.rank(x)
                        tmp[r].append((n, p))
                    n += 1
        T = array.array('I', [])
        U = array.array('I', [])
        V = array.array('i', [])
        t = 0
        for nps in tmp:
            T.append(t)
            t += len(nps)
            for (n, p) in nps:
                U.append(n)
                V.append(p)
        T.append(t)
        del tmp

        gfn = opts['<genes>']
        with container(gfn, 'w') as z:
            z.meta['K'] = K
            z.meta['S'] = S.count()
            write64(z, S.xs, 'S')
            z.meta['T'] = len(T)
            write64(z, T, 'T')
            z.meta['U'] = len(U)
            write32(z, U, 'U')
            z.meta['V'] = len(V)
            write32s(z, V, 'V')
            z.meta['lens'] = lens
            z.meta['qacgt'] = qacgt
            z.meta['nms'] = nms
            z.meta['seqs'] = seqs

        return

    print >> sys.stderr, "loading..."

    gfn = opts['<genes>']
    with container(gfn, 'r') as z:
        K = z.meta['K']
        S = array.array('L', read64(z, 'S', z.meta['S']))
        S = sparse(2*K, S)
        T = array.array('L', read64(z, 'T', z.meta['T']))
        U = array.array('I', read32(z, 'U', z.meta['U']))
        V = array.array('i', read32s(z, 'V', z.meta['V']))
        lens = z.meta['lens']
        qacgt = z.meta['qacgt']
        nms = z.meta['nms']
        seqs = z.meta['seqs']

    print >> sys.stderr, "done."

    for fn in opts['<input>']:
        L = array.array('B', [0 for i in xrange(S.count())])
        Y = array.array('L', [0 for i in xrange(S.count())])
        with container(fn, 'r') as z:
            sacgt = z.meta['acgt']
            xs = readKmers(z)
            X = array.array('L', xs)
        M = len(X)
        resolveAll(K, S, L, Y, X)
        X = sparse(2*K, X)

        g = sum([qp*sp for (qp, sp) in zip(qacgt, sacgt)])
        print >> sys.stderr, "g =", g
        nm = [null(g, M, j) for j in range(0, K+1)]

        # counts for computing distribution of prefix lengths
        cnt = [[0 for j in xrange(K+1)] for i in xrange(len(nms))]

        # the k-mers that we pulled by lcp from the sample
        # for each position of each query.
        P = [array.array('L', [0 for j in xrange(lens[i] - K + 1)]) for i in xrange(len(lens))]

        # the length of the lcp for each position of each query.
        Q = [array.array('B', [0 for j in xrange(lens[i] - K + 1)]) for i in xrange(len(lens))]

        for i in xrange(S.count()):
            for j in xrange(T[i], T[i+1]):
                n = U[j]
                p = V[j]
                y = Y[i]
                l = L[i]
                cnt[n][l] += 1
                if p > 0:
                    p -= 1
                else:
                    p = -(p + 1)
                    y = rc(K, y)
                if l > Q[n][p]:
                    Q[n][p] = l
                    P[n][p] = y

        for i in xrange(len(nms)):
            # iterate over the queries

            qc = math.log(K*0.05/float(lens[i] - K + 1)/2)

            # Link up "de Bruijn" sequences
            m = (1 << (2*K - 2)) - 1
            py = 0
            u = unionfind()
            for j in xrange(lens[i] - K + 1):
                x = P[i][j]
                y = x >> 2
                if j > 0:
                    d = ham(py, y)
                    if d == 0:
                        u.union(j-1, j)
                py = x & m

            # Gather up the de Bruin fragments
            udx = {}
            for j in xrange(lens[i] - K + 1):
                v = u.find(j)
                if v not in udx:
                    udx[v] = []
                udx[v].append(j)

            # Index the left hand k-mers
            idxLhs = {}
            kx = []
            for (jx, js) in udx.iteritems():
                q = 0
                for j in js:
                    q += math.log1p(-nm[Q[i][j]])
                if q > math.log(0.05/len(js)):
                    continue
                kx.append((-len(js), jx))
                idxLhs[P[i][js[0]]] = jx
            kx.sort()

            # Attempt to link up fragments
            links = {}
            for (_, jx) in kx:
                jR = udx[jx][-1]
                if jR == lens[i] - K + 1:
                    continue
                x = P[i][jR]
                xs = []
                lnk = None
                for k in xrange(100):
                    ys = succ(K, X, x)
                    if len(ys) != 1:
                        break
                    x = ys[0]
                    if x in idxLhs:
                        lnk = idxLhs[x]
                        break
                    xs.append(x)
                if lnk is not None:
                    links[jx] = xs
                    u.union(jx, lnk)

            # Gather up the linked fragments
            vdx = {}
            for j in [jx for (_, jx) in kx]:
                v = u.find(j)
                if v not in vdx:
                    vdx[v] = []
                vdx[v].append(j)

            res = []
            for (jxx, jxs) in vdx.iteritems():
                # Order the gragments by start position
                fs = [(udx[jx][0], jx) for jx in jxs]
                fs.sort()
                sxs = []
                for fj in xrange(len(fs)):
                    (_, jx) = fs[fj]
                    beg = udx[jx][0]
                    end = udx[jx][-1] + 1
                    if fj == 0:
                        for j in xrange(beg):
                            sxs.append((0, 0))
                    xs = links.get(jx, None)
                    for j in xrange(beg, end):
                        x = P[i][j]
                        l = Q[i][j]
                        sxs.append((x, l))
                    if xs:
                        for x in xs:
                            sxs.append((x, 27))
                    else:
                        if fj < len(fs) - 1:
                            nxt = fs[fj+1][0]
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
                p = None
                inf = False
                for j in xrange(len(seq)):
                    b = 0
                    for k in xrange(4):
                        if seq[j][k] < qc:
                            b |= 1 << k
                    ax.append(fasta(b))
                    ssj = sum(seq[j])
                    if p is None:
                        p = ssj
                    else:
                        p = logAdd(p, ssj)
                    if ssj > -1e-300:
                        inf = True
                dst = counts2cdf(cnt[i])
                (_, kd) = ksDistance2(dst, nm)
                df = math.ceil(len(seq)/float(K))
                if inf:
                    q = 1e300
                    pv = 0.0
                else:
                    q = 2*math.exp(p)
                    pv = chi2(df, q)
                res.append((pv, q, kd, ''.join(ax)))

            if len(res) == 0:
                continue

            res.sort()
            if res[0][0] < -2:
                #ed = lev(seqs[i], res[0][2])
                ed = 0
                pv = res[0][0]/math.log(10)
                c2 = res[0][1]
                kd = res[0][2]
                a = res[0][3]
                print '%d\t%d\t%d\t%g\t%g\t%g\t%s\t%s' % (i, lens[i], len(a), kd, c2, pv, nms[i], a)
            sys.stdout.flush()

if __name__ == '__main__':
    main(sys.argv[1:])
