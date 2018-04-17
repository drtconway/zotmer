"""
zot alu-finder - find ALU insertion points

Usage:
    zot alu-finder [options] <regions> <input>...

Options:
    -k K            value of k to use [default: 25]
    -g PATH         directory of FASTQ reference sequences
    -C INT          coverage cutoff value [default: 5]
    -v              produce verbose output
"""

import os.path
import sys

import docopt
import yaml

from pykmer.basics import ham, kmersWithPosList, kmersWithPosLists, murmer, rc, render
from pykmer.file import openFile, readFasta, readFastq
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
        spurs = [[ref[p0]]]
        short = []
        p = p0 + 1
        while p in Z and len(spurs) > 0:
            spurs1 = []
            for spur in spurs:
                x = spur[-1]
                ext = False
                for (y,c) in Z[p].items():
                    if follow(K, x, y):
                        if p in ref and ref[p] == y:
                            continue
                        spurs1.append(spur + [y])
                        ext = True
                if not ext:
                    short.append(spur)
            spurs = spurs1
            p += 1
        yield (p0, sorted(short + spurs))

def reverseSpurs(K, ref, Z):
    for p0 in sorted(ref.keys()):
        spurs = [[ref[p0]]]
        short = []
        p = p0 - 1
        while p in Z and len(spurs) > 0:
            spurs1 = []
            for spur in spurs:
                x = spur[0]
                ext = False
                for (y,c) in Z[p].items():
                    if follow(K, y, x):
                        if p in ref and ref[p] == y:
                            continue
                        spurs1.append([y] + spur)
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

def shiftForwardSpur(ref, p, spur):
    i = 0
    while i < 5:
        yield (p, spur, i)
        i += 1
        p -= 1
        if p not in ref:
            break
        spur = [ref[p]] + spur

def shiftReverseSpur(ref, p, spur):
    i = 0
    while i < 5:
        yield (p, spur, i)
        i += 1
        p += 1
        if p not in ref:
            break
        spur = spur + [ref[p]]

def contigs(K, Z, p0):
    ctgs = []
    for (x,c) in Z[p0].items():
        if c < 10:
            continue
        ctgs.append(([x],[c]))

    short = []
    pi = p0 + 1
    while pi in Z:
        nctgs = []
        for ctg in ctgs:
            y = ctg[0][-1]
            ext = False
            for (x,c) in Z[pi].items():
                if c < 10:
                    continue
                if follow(K, y, x):
                    ext = True
                    nctgs.append((ctg[0] + [x], ctg[1] + [c]))
            if not ext:
                short.append(ctg)
        ctgs = nctgs
        pi += 1
    return ctgs

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

class ResultCollator(object):
    def __init__(self):
        self.res = {}

    def collate(self, K, C, z, q, ps, pxc):
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

            scs = sorted(cs)
            q1 = quant(scs, 0.1)
            q5 = quant(scs, 0.5)
            q9 = quant(scs, 0.9)

            if q1 <= C:
                continue
            if len(xs) == 0:
                continue

            r = {}
            r['cnts'] = cs
            r['q10'] = q1
            r['q50'] = q5
            r['q90'] = q9

            seq = renderPath(K, xs)

            if mp < 0:
                d = -1
                r['side'] = 'before'
                r['seq'] = seq[:-K1]
                r['anc'] = seq[-K1:].lower()
            else:
                d = 1
                r['side'] = 'after'
                r['seq'] = seq[K1:]
                r['anc'] = seq[:K1].lower() 

            if len(r['seq']) < 10:
                continue

            #if len(r['seq']) < K:
            #    continue

            if z not in self.res:
                self.res[z] = {}
            if q not in self.res[z]:
                self.res[z][q] = {}
            if d not in self.res[z][q]:
                self.res[z][q][d] = []
            self.res[z][q][d].append(r)

    def dump(self, out):
        for z in sorted(self.res.keys()):
            Z = self.res[z]
            for q in sorted(Z.keys()):
                Q = Z[q]
                for s in sorted(Q.keys()):
                    for v in Q[s]:
                        print '%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s' % (z, q, v['side'], v['q10'], v['q50'], v['q90'], v['anc'], v['seq'], ','.join(map(str, v['cnts'])))

    def output(self, out):
        hdrShown = False
        for z in sorted(self.res.keys()):
            Z = self.res[z]
            for q0 in sorted(Z.keys()):
                Q0 = Z[q0]
                if 1 not in Q0:
                    continue
                q1 = q0+1
                if q1 not in Z:
                    continue
                Q1 = Z[q1]
                if -1 not in Q1:
                    continue

                for i in range(len(Q0[1])):
                    b = Q0[1][i]
                    for j in range(len(Q1[-1])):
                        a = Q1[-1][j]

                        hdr = ['accession', 'after', 'before']
                        fmt = ['%s', '%d', '%d']
                        val = [z, q0+1, q1+1]

                        hdr += ['lhsQ10', 'lhsQ50', 'lhsQ90']
                        fmt += ['%d', '%d', '%d']
                        val += [b['q10'], b['q50'], b['q90']]

                        hdr += ['rhsQ10', 'rhsQ50', 'rhsQ90']
                        fmt += ['%d', '%d', '%d']
                        val += [a['q10'], a['q50'], a['q90']]

                        hdr += ['lhsAnc', 'rhsAnc']
                        fmt += ['%s', '%s']
                        val += [b['anc'], a['anc']]

                        hdr += ['lhsSeq', 'rhsSeq']
                        fmt += ['%s', '%s']
                        val += [b['seq'], a['seq']]
                        
                        if not hdrShown:
                            print >> out, '\t'.join(hdr)
                            hdrShown = True
                        print >> out, '\t'.join(fmt) % tuple(val)

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
                vt = 0.05 * sum([c for (x,c) in vs])
                for (x,c) in vs:
                    if c < vt or c < C:
                        killX.add(x)
            #for x in acc[z][p].keys():
            #    if acc[z][p][x] < C:
            #        killX.add(x)
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

    print '\t'.join(['accession', 'after', 'before', 'shift', 'lhsAnc', 'rhsAnc', 'lhsSeq', 'rhsSeq'])
    for z in sorted(acc.keys()):
        Z = acc[z]
        ref = refTbl[z]
        aft = dict(forwardSpurs(K, ref, Z))
        bef = dict(reverseSpurs(K, ref, Z))

        print >> sys.stderr, z, len(aft), len(bef)

        scoredAft = {}
        for p in sorted(aft.keys()):
            for spur in aft[p]:
                #if len(spur) < 2*K:
                if len(spur) < 10:
                    continue
                for (q, xs, v) in shiftForwardSpur(ref, p, spur):
                    q += K-1
                    if q not in scoredAft:
                        scoredAft[q] = []
                    seq = renderPath(K, xs)
                    anc = seq[:K]
                    ins = seq[K:]
                    scoredAft[q].append((v, anc, ins))

        scoredBef = {}
        for p in sorted(bef.keys()):
            for spur in bef[p]:
                if len(spur) < K:
                    continue
                for (q, xs, v) in shiftReverseSpur(ref, p, spur):
                    if q not in scoredBef:
                        scoredBef[q] = []
                    seq = renderPath(K, xs)
                    anc = seq[-K:]
                    ins = seq[:-K]
                    scoredBef[q].append((v, anc, ins))

        for p0 in sorted(scoredAft.keys()):
            p1 = p0 + 1
            if p1 not in scoredBef:
                continue
            for (aftV, aftAnc, aftIns) in scoredAft[p0]:
                for (befV, befAnc, befIns) in scoredBef[p1]:
                    v = aftV + befV
                    print '%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s' % (z, p0, p1, v, aftAnc, befAnc, aftIns, befIns)

if __name__ == '__main__':
    main(sys.argv[1:])

