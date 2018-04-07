"""
zot alu-finder - find ALU insertion points

Usage:
    zot alu-finder [options] <regions> <output-prefix> <input>...

Options:
    -k K            value of k to use [default: 25]
    -g PATH         directory of FASTQ reference sequences
    -C INT          coverage cutoff value [default: 5]
    -v              produce verbose output
    -z              compress the output reads
"""

import sys

import docopt
import yaml

from pykmer.basics import ham, kmersWithPosList, kmersWithPosLists, rc, render
from pykmer.file import openFile, readFasta, readFastq
from zotmer.library.hgvs import hg19ToRefSeq, makeHGVS, refSeq2Hg19
from zotmer.library.reads import reads

def outputFile(nm):
    if nm is None:
        return sys.stdout
    return openFile(nm, 'w')

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

def quant(xs, q):
    z = len(xs)
    i = int(q*z)
    return xs[i]

def flatten(xss):
    ys = []
    for xs in xss:
        ys += xs
    return ys

def hits(idx, K, xps):
    hc = 0
    mc = 0
    l = {}
    mx = {}
    hx = {}
    for (x,p) in xps:
        p -= 1
        if x not in idx:
            mc += 1
            if x not in mx:
                mx[x] = []
            mx[x].append(p)
            continue
        hc += 1
        for (z,q) in idx[x]:
            r = q - p
            k = (z,r)
            if k not in l:
                l[k] = 0
            l[k] += 1
            if x not in hx:
                hx[x] = []
            hx[x].append(p)

    if len(l) == 0:
        return None

    pMin = min(flatten(hx.values()))
    pMax = max(flatten(hx.values()))

    pre = []
    aft = []
    for (z,r) in l.keys():
        qL = r + pMin
        qR = r + pMax + K - 1

        for (x, ps) in mx.items():
            for p in ps:
                if p < pMin:
                    pre.append((z, qL, p-pMin, x))
                if p > pMax:
                    aft.append((z, qR, p-pMax, x))
    pre.sort()
    aft.sort()

    if len(pre) == 0 and len(aft) == 0:
        return None

    return (pre, aft)

def accumulateHits(hx, acc):
    if hx is None:
        return

    for vs in hx:
        for k in vs:
            (z,q,p,x) = k
            gk = (z,q)
            if gk not in acc:
                acc[gk] = {}
            if p not in acc[gk]:
                acc[gk][p] = {}
            if x not in acc[gk][p]:
                acc[gk][p][x] = 0
            acc[gk][p][x] += 1

def groupKmers(K, xcs):
    M = 10
    if len(xcs) == 1:
        return xcs
    grps = []
    yns = []
    for (x,c) in xcs.items():
        if c < M:
            yns.append((x,c))
        else:
            grps.append([x,[(x,c)]])
    for (y,n) in yns:
        used = False
        for i in range(len(grps)):
            x = grps[i][0]
            d = ham(y,x)
            if d < 2:
                grps[i][1].append((y,n))
                used = True
                break
        if not used:
            grps.append([y,[(y,n)]])

    res = {}
    for (x,yns) in grps:
        bs = [[0,0,0,0] for j in range(K)]
        for (y,n) in yns:
            for j in range(K):
                bs[j][(y >> (2*j))&3] += n
        m = None
        y = 0
        for j in range(K):
            vs = [(bs[j][b],b) for b in range(4)]
            vs.sort(reverse=True)
            b = vs[0][1]
            y |= b << (2*j)
            if m is None:
                m = vs[0][0]
            else:
                m = min(m, vs[0][0])
        res[y] = m
    return res

def follow(K, x, y):
    K1 = K - 1
    M1 = (1 << (2*K1)) - 1
    return (x & M1) == (y >> 2)

def renderPath(K, xs):
    if len(xs) == 0:
        return ''
    r = [render(K, xs[0])]
    for x in xs[1:]:
        r.append("ACGT"[x&3])
    return ''.join(r)

def contig(K, ps, pxc):
    ctgs = []

    p0 = ps[0]
    for (x,c) in pxc[p0].items():
        ctgs.append(([x],[c]))

    for i in range(1, len(ps)):
        pi = ps[i]
        nctgs = []
        for ctg in ctgs:
            y = ctg[0][-1]
            for (x,c) in pxc[pi].items():
                if follow(K, y, x):
                    nctgs.append((ctg[0] + [x], ctg[1] + [c]))
        ctgs = nctgs
    return ctgs


def output(K, C, z, q, ps, pxc):
    if len(ps) == 0:
        return
    mp = min(ps)

    K1 = K-1

    for (xs,cs) in contig(K, ps, pxc):
        assert len(xs) == len(cs)
        if False: # Trim the ends
            while len(cs) > 0 and cs[-1] <= C:
                xs.pop()
                cs.pop()

            xs = xs[::-1]
            cs = cs[::-1]
            while len(cs) > 0 and cs[-1] <= C:
                xs.pop()
                cs.pop()
            xs = xs[::-1]
            cs = cs[::-1]

            if not all([c > C for c in cs]):
                continue

        cs.sort()
        q1 = quant(cs, 0.1)
        q5 = quant(cs, 0.5)
        q9 = quant(cs, 0.9)

        if q1 <= C:
            continue
        if len(xs) == 0:
            continue

        seq = renderPath(K, xs)

        if mp < 0:
            d = 'before'
            seq = seq[:-K1] + seq[-K1:].lower()
        else:
            d = 'after'
            seq = seq[:K1].lower() + seq[K1:]

        if len(seq) < K:
            continue

        print '%s\t%d\t%s\t%g\t%g\t%g\t%s' % (z, q, d, q1, q5, q9, seq)

    #for p in ps:
    #    w = ' ' * (p - mp)
    #    for (x,c) in sorted(xcs.items()):
    #        if c < C:
    #            continue
    #        print '%s\t%d\t%d\t%d\t%s%s' % (z, q, p, c, w, render(K, x))

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    verbose = opts['-v']

    K = int(opts['-k'])

    C = int(opts['-C'])

    d = "."
    if opts['-g']:
        d = opts['-g']
    sf = SequenceFactory(d)

    with openFile(opts['<regions>']) as f:
        R = readBED(f)

    refIdx = {}
    zoneIdx = {}
    for (acc, zones) in R.items():
        accSeq = sf[acc]
        for (s, e, nm) in zones:
            zoneIdx[nm] = (acc, s, e)
            seq = accSeq[s:e]
            for (x,p) in kmersWithPosList(K, seq, False):
                #p = s + p - 1
                p -= 1
                if x not in refIdx:
                    refIdx[x] = []
                refIdx[x].append((nm,p))

    acc = {}
    pile = {}
    for itm in reads(opts['<input>'], K=K, paired=True, reads=True, kmers=False, verbose=verbose):
        rdL = itm.reads[0]
        zL = len(rdL)
        (fwdL, revL) = kmersWithPosLists(K, rdL[1])
        fwdLHits = hits(refIdx, K, fwdL)
        revLHits = hits(refIdx, K, revL)

        rdR = itm.reads[1]
        zR = len(rdR)
        (fwdR, revR) = kmersWithPosLists(K, rdR[1])
        fwdRHits = hits(refIdx, K, fwdR)
        revRHits = hits(refIdx, K, revR)

        if fwdLHits is not None:
            accumulateHits(fwdLHits, acc)
        if revRHits is not None:
            accumulateHits(revRHits, acc)
        if revLHits is not None:
            accumulateHits(revLHits, acc)
        if fwdRHits is not None:
            accumulateHits(fwdRHits, acc)

    for gk in sorted(acc.keys()):
        (z,q) = gk
        pxc = acc[gk]
        ps = pxc.keys()
        ls = []
        rs = []
        for p in sorted(ps):
            if p > 0:
                rs.append(p)
            else:
                ls.append(p)

        if q > 0:
            output(K, C, z, q, ls, pxc)

        (a,s,e) = zoneIdx[z]
        zl = e - s - K
        if q < zl:
            output(K, C, z, q, rs, pxc)

if __name__ == '__main__':
    main(sys.argv[1:])

