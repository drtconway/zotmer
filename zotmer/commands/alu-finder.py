"""
zot alu-finder - find ALU insertion points

Usage:
    zot alu-finder [options] <regions> <input>...

Options:
    -k K            value of k to use [default: 25]
    -g PATH         directory of FASTQ reference sequences
    -C INT          coverage cutoff value [default: 5]
    -L INT          Minimum length of spurs [default: 29]
    -r              produce raw spurs
    -S INT          Maximum distance to shift spurs [default: 5]
    -V FLOAT        minimum relative frequency for insertions [default: 0.05]
    -v              produce verbose output
"""

import os.path
import sys

import docopt
import yaml

from zotmer.library.basics import ham, kmersWithPosList, kmersWithPosLists, murmer, rc, render
from zotmer.library.file import openFile, readFasta, readFastq
from zotmer.library.hgvs import hg19ToRefSeq, makeHGVS, refSeq2Hg19
from zotmer.library.reads import reads

def succ(K, x):
    M = (1 << (2*K)) - 1
    y0 = (x << 2) & M
    r = []
    for i in range(4):
        r.append(y0 + i)
    return r

def pred(K, x):
    S = 2*(K-1)
    y0 = x >> 2
    r = []
    for i in range(4):
        r.append(y0 + (i << S))
    return r

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

            pth = self.home + '/' + h + '.fa'
            if not os.path.exists(pth):
                pth += '.gz'

            with openFile(pth) as f:
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

def hits(idx, K, xps, hx = None):
    loc = {}
    for (x,p) in xps:
        p -= 1
        if x not in idx:
            continue
        for (z,q) in idx[x]:
            r = q - p
            k = (z,r)
            if k not in loc:
                loc[k] = 0
            loc[k] += 1

    if len(loc) == 0:
        return None

    if hx is None:
        hx = {}
    for k in loc.keys():
        (z,r) = k
        if z not in hx:
            hx[z] = {}
        for (x,p) in xps:
            p -= 1
            q = r + p
            if q not in hx[z]:
                hx[z][q] = {}
            if x not in hx[z][q]:
                hx[z][q][x] = 0
            hx[z][q][x] += 1
    return hx

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

def forwardSpurs(K, ref, Z):
    for p0 in sorted(ref.keys()):
        x0 = ref[p0]
        if p0 not in Z or x0 not in Z[p0]:
            continue
        c0 = Z[p0][x0]

        short = []
        spurs = [[(x0, c0)]]
        p = p0 + 1
        while p in Z and len(spurs) > 0:
            spurs1 = []
            for spur in spurs:
                x = spur[-1][0]
                ext = False
                for (y,c) in Z[p].items():
                    if follow(K, x, y):
                        if p in ref and ref[p] == y:
                            continue
                        spurs1.append(spur + [(y,c)])
                        ext = True
                if not ext:
                    short.append(spur)
            spurs = spurs1
            p += 1
        yield (p0, sorted(short + spurs))

def reverseSpurs(K, ref, Z):
    for p0 in sorted(ref.keys()):
        x0 = ref[p0]
        if p0 not in Z or x0 not in Z[p0]:
            continue
        c0 = Z[p0][x0]

        short = []
        spurs = [[(x0,c0)]]
        p = p0 - 1
        while p in Z and len(spurs) > 0:
            spurs1 = []
            for spur in spurs:
                x = spur[0][0]
                ext = False
                for (y,c) in Z[p].items():
                    if follow(K, y, x):
                        if p in ref and ref[p] == y:
                            continue
                        spurs1.append([(y,c)] + spur)
                        ext = True
                if not ext:
                    short.append(spur)
            spurs = spurs1
            p -= 1
        yield (p0, sorted(short + spurs))

def sig(xs):
    assert len(xs) > 0
    h = xs[0]
    for x in xs[1:]:
        h = murmer(x, h)
    return h

def shiftForwardSpur(ref, Z, S, p, spur):
    i = 0
    while i < S:
        yield (p, spur, i)
        i += 1
        p -= 1
        if p not in ref:
            break
        x = ref[p]
        if p not in Z or x not in Z[p]:
            break
        c = Z[p][x]
        spur = [(x,c)] + spur

def shiftReverseSpur(ref, Z, S, p, spur):
    i = 0
    while i < S:
        yield (p, spur, i)
        i += 1
        p += 1
        if p not in ref:
            break
        x = ref[p]
        if p not in Z or x not in Z[p]:
            break
        c = Z[p][x]
        spur = spur + [(x,c)]

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    verbose = opts['-v']

    K = int(opts['-k'])

    C = int(opts['-C'])

    L = int(opts['-L'])

    raw = opts['-r']

    S = int(opts['-S'])

    V = float(opts['-V'])

    d = "."
    if opts['-g']:
        d = opts['-g']
    sf = SequenceFactory(d)

    with openFile(opts['<regions>']) as f:
        R = readBED(f)

    refTbl = {}
    refIdx = {}
    zoneIdx = {}
    for (acc, zones) in R.items():
        accSeq = sf[acc]
        for (s, e, nm) in zones:
            zoneIdx[nm] = (acc, s, e)
            seq = accSeq[s-1:e]
            if nm not in refTbl:
                refTbl[nm] = {}
            for (x,p) in kmersWithPosList(K, seq, False):
                p -= 1
                p += s
                refTbl[nm][p] = x
                if x not in refIdx:
                    refIdx[x] = []
                refIdx[x].append((nm,p))

    acc = {}
    for itm in reads(opts['<input>'], K=K, paired=True, reads=True, kmers=False, verbose=verbose):
        rdL = itm.reads[0]
        zL = len(rdL)
        (fwdL, revL) = kmersWithPosLists(K, rdL[1])
        fwdLHits = hits(refIdx, K, fwdL, acc)
        revLHits = hits(refIdx, K, revL, acc)

        rdR = itm.reads[1]
        zR = len(rdR)
        (fwdR, revR) = kmersWithPosLists(K, rdR[1])
        fwdRHits = hits(refIdx, K, fwdR, acc)
        revRHits = hits(refIdx, K, revR, acc)

    killZ = set([])
    for z in acc.keys():
        killP = set([])
        for p in acc[z].keys():
            killX = set([])
            vv = {}
            for x in acc[z][p].keys():
                y = x >> 2
                if y not in vv:
                    vv[y] = []
                vv[y].append((x, acc[z][p][x]))
            for vs in vv.values():
                vt = V * sum([c for (x,c) in vs])
                for (x,c) in vs:
                    if c < vt or c < C:
                        killX.add(x)
            for x in killX:
                del acc[z][p][x]
            if len(acc[z][p]) == 0:
                killP.add(p)
        for p in killP:
            del acc[z][p]
        if len(acc[z]) == 0:
            killZ.add(z)
    for z in killZ:
        del acc[z]

    if raw:
        print '\t'.join(['chrom', 'pos', 'side', 'label', 'anchor', 'insSeq'])
    else:
        print '\t'.join(['chrom', 'after', 'before', 'label', 'rhsShift', 'lhsShift', 'lhsAnc', 'rhsAnc', 'lhsSeq', 'rhsSeq'])

    for z in sorted(acc.keys()):
        (ch, st, en) = zoneIdx[z]

        Z = acc[z]
        ref = refTbl[z]
        aft = dict(forwardSpurs(K, ref, Z))
        bef = dict(reverseSpurs(K, ref, Z))

        scoredAft = {}
        for p in sorted(aft.keys()):
            if p+K-1 == en:
                continue

            for spur in aft[p]:

                if len(spur) < L:
                    continue

                if raw:
                    (xs, cs) = zip(*spur)
                    seq = renderPath(K, xs)
                    anc = seq[:K]
                    ins = seq[K:]
                    print '%s\t%d\t%s\t%s\t%s\t%s\t%s' % (ch, p+K-1, 'after', z, anc, ins, ','.join(map(str, cs)))
                    continue

                for (q, xcs, v) in shiftForwardSpur(ref, Z, S, p, spur):
                    q += K-1
                    if q not in scoredAft:
                        scoredAft[q] = []
                    (xs, cs) = zip(*xcs)
                    seq = renderPath(K, xs)
                    anc = seq[:K]
                    ins = seq[K:]
                    scoredAft[q].append((v, anc, ins, cs))

        scoredBef = {}
        for p in sorted(bef.keys()):
            if p == st:
                continue

            for spur in bef[p]:
                if len(spur) < L:
                    continue

                if raw:
                    (xs, cs) = zip(*spur)
                    seq = renderPath(K, xs)
                    anc = seq[-K:]
                    ins = seq[:-K]
                    print '%s\t%d\t%s\t%s\t%s\t%s\t%s' % (ch, p, 'before', z, anc, ins, ','.join(map(str, cs)))
                    continue

                for (q, xcs, v) in shiftReverseSpur(ref, Z, S, p, spur):
                    if q not in scoredBef:
                        scoredBef[q] = []
                    (xs, cs) = zip(*xcs)
                    seq = renderPath(K, xs)
                    anc = seq[-K:]
                    ins = seq[:-K]
                    scoredBef[q].append((v, anc, ins, cs))

        for p0 in sorted(scoredAft.keys()):
            p1 = p0 + 1
            if p1 not in scoredBef:
                continue
            for (aftV, aftAnc, aftIns, aftCov) in scoredAft[p0]:
                for (befV, befAnc, befIns, befCov) in scoredBef[p1]:
                    if befAnc in aftIns or aftAnc in befIns:
                        continue
                    v = aftV + befV
                    print '%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s' % (ch, p0, p1, z, aftV, befV, aftAnc, befAnc, aftIns, befIns)

if __name__ == '__main__':
    main(sys.argv[1:])

