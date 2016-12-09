from base import Cmd

from pykmer.basics import fasta, kmersWithPos, ham, lcp, rc, render
from pykmer.file import readFasta
import pykmer.kset as kset
import pykmer.kfset as kfset
from pykmer.stats import counts2cdf, ksDistance2, log1mexp
from pykmer.uf import uf
from pykmer.sparse import sparse
from pykmer.exceptions import MismatchedK

import cPickle
import array
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

def succ(K, xs, x):
    m = (1 << (2*K)) - 1
    y0 = (x << 2) & m
    y1 = y0 + 4
    r0 = xs.rank(y0)
    r1 = xs.rank(y1)
    return [xs.select(i) for i in xrange(r0, r1)]

def null(g, s, j):
    if j == 0:
        return 0.0
    return math.exp(s*math.log1p(-math.pow(g,j)))

def logNull(g, s, j):
    return s*math.log1p(-math.pow(g,j))

def uniq(xs):
    px = None
    for x in xs:
        if x != px:
            yield x
            px = x

class Scan(Cmd):
    def run(self, opts):
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

            meta = (K, S.count(), len(T), len(U), len(V), qacgt)
            gfn = opts['<genes>']
            with open(gfn + '.meta', 'w') as f:
                cPickle.dump(meta, f)
            with open(gfn + '.kmers', 'w') as f:
                S.xs.tofile(f)
            with open(gfn + '.sums', 'w') as f:
                T.tofile(f)
            with open(gfn + '.nums', 'w') as f:
                U.tofile(f)
            with open(gfn + '.poss', 'w') as f:
                V.tofile(f)
            with open(gfn + '.lens', 'w') as f:
                cPickle.dump(lens, f)
            with open(gfn + '.names', 'w') as f:
                cPickle.dump(nms, f)
            with open(gfn + '.seqs', 'w') as f:
                cPickle.dump(seqs, f)

            return

        print >> sys.stderr, "loading..."

        gfn = opts['<genes>']
        with open(gfn + '.meta') as f:
            meta = cPickle.load(f)
        (K, zS, zT, zU, zV, qacgt) = meta
        with open(gfn + '.kmers') as f:
            S = array.array('L', [])
            S.fromfile(f, zS)
            S = sparse(2*K, S)
        with open(gfn + '.sums') as f:
            T = array.array('I', [])
            T.fromfile(f, zT)
        with open(gfn + '.nums') as f:
            U = array.array('I', [])
            U.fromfile(f, zU)
        with open(gfn + '.poss') as f:
            V = array.array('i', [])
            V.fromfile(f, zV)
        with open(gfn + '.lens') as f:
            lens = cPickle.load(f)
        with open(gfn + '.names') as f:
            nms = cPickle.load(f)
        with open(gfn + '.seqs') as f:
            seqs = cPickle.load(f)

        print >> sys.stderr, "done."

        for fn in opts['<input>']:
            L = array.array('B', [0 for i in xrange(S.count())])
            Y = array.array('L', [0 for i in xrange(S.count())])
            (meta, xs) = kfset.read(fn)
            sacgt = [0, 0, 0, 0]
            X = array.array('L', [])
            M = 0
            for (x,_) in xs:
                X.append(x)
                sacgt[x&3] += 1
                resolveLcp(K, S, L, Y, x)
                M += 1
            sacgt = [float(c)/float(M) for c in sacgt]
            X = sparse(2*K, X)

            g = sum([qp*sp for (qp, sp) in zip(qacgt, sacgt)])
            print >> sys.stderr, "g =", g
            nm = [null(g, M, j) for j in range(0, K+1)]

            cnt = [[0 for j in xrange(K+1)] for i in xrange(len(nms))]
            P = [array.array('L', [0 for j in xrange(lens[i] - K + 1)]) for i in xrange(len(lens))]
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
                qc = math.log(0.05/(lens[i] - K + 1)/2)

                # Link up "de Bruijn" sequences
                m = (1 << (2*K - 2)) - 1
                py = 0
                u = uf()
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
                if res[0][0] < qc:
                    #ed = lev(seqs[i], res[0][2])
                    ed = 0
                    print '%d\t%d\t%d\t%g\t%g\t%d\t%s\t%s' % (i, lens[i], len(res[0][2]), res[0][1], res[0][0], ed, nms[i], res[0][2])
                #for j in xrange(len(res)):
                #    if j > 0 and res[j][0] > qc:
                #        break
                #    print '%d\t%d\t%d\t%g\t%g\t%s\t%s' % (i, lens[i], len(res[j][2]), res[j][1], res[j][0], nms[i], res[j][2])
                #    if res[j][0] > qc:
                #        break
                sys.stdout.flush()

def add(cmds):
    cmds['scan'] = Scan()
