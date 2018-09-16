"""
zot exonic - estimate coverage on transcripts.

Usage:
    zot exonic -X [options] <RefGene>
    zot exonic [options] <index> <input>...

Options:
    -g PATH         directory of FASTQ reference sequences
    -k K            value of k to use [default: 25]
    -v              produce verbose output
    -S SLICE        process a subset of genes given a n/N [default: 1/1]
    -X              extract exons from a BED file
"""

import array
import math
import random
import sys
import zipfile

import docopt
import yaml

from pykmer.basics import ham, kmers, kmersList, kmersWithPosList, murmer, rc, render
from pykmer.file import openFile, readFasta, readFastq
from pykmer.misc import uniq
from zotmer.library.debruijn import debruijn_intersection
from zotmer.library.hgvs import hg19ToRefSeq, refSeq2Hg19
from zotmer.library.reads import reads

def ham1(K, x):
    xs = []
    for i in xrange(K):
        for j in xrange(3):
            xs.append(x ^ ((j + 1) << (2*i)))
    xs.sort()
    return xs

def readRefGene(f):
    res = {}
    for l in f:
        t = l.split()
        tx = t[1]
        c = t[2]
        if c in hg19ToRefSeq:
            c = hg19ToRefSeq[c]
        strand = t[3]
        sts = t[9].split(',')
        ens = t[10].split(',')
        nm = t[12]
        if c not in res:
            res[c] = {}
        if nm not in res[c]:
            res[c][nm] = set([])
        for (st,en) in zip(sts, ens):
            if len(st) == 0:
                continue
            ex = (strand, int(st), int(en))
            res[c][nm].add(ex)
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

def hits(idx, fwd, rev):
    hx = {}
    for x in fwd:
        if x not in idx:
            continue
        for k in idx[x]:
            if k not in hx:
                hx[k] = 0
            hx[k] += 1
    for x in rev:
        if x not in idx:
            continue
        for k in idx[x]:
            if k not in hx:
                hx[k] = 0
            hx[k] += 1
    return hx

def bridge(K, st, en, X):
    K1 = K - 1
    M1 = (1 << (2*K1)) - 1

    xs = sorted(X.keys())

    # Create an index of the leading K-1mers
    #
    xIdx = {}
    for x in xs:
        x1 = x >> 2
        if x1 not in xIdx:
            xIdx[x1] = []
        xIdx[x1].append(x)

    vec = [sorted([x for x in st if x in X])]
    for j in range(K-1):
        # get the trailing suffixes
        #
        ys = [x & M1 for x in vec[-1]]
        ys.sort()
        uniq(ys)

        # use them to find following k-mers
        #
        zs = []
        for y in ys:
            if y in xIdx:
                zs += xIdx[y]
        zs.sort()
        uniq(zs)
        if len(zs) == 0:
            # Short circuit
            return None
        vec.append(zs)
    vec += [sorted([x for x in en if x in X])]

    # Arrange the vector so we cycle forwards then backwards
    #
    fwd = range(1, len(vec))
    rev = fwd[::-1]
    order = [fwd, rev]

    changed = True
    flip = 0
    while changed:
        changed = False
        for i in order[flip]:
            xs = vec[i-1]
            ys = vec[i]
            (vs, ws) = debruijn_intersection(K, xs, ys)
            if len(vs) == 0 or len(ws) == 0:
                return None
            if len(xs) != len(vs):
                vec[i-1] = vs
                changed = True
            if len(ys) != len(ws):
                vec[i] = ws
                changed = True
        flip = 1 - flip
    return vec

def paths(K, vec):
    if vec is None:
        return

    V = len(vec)

    M1 = (1 << (2*(K-1))) - 1

    stk = []
    for x in vec[0]:
        stk.append([x])

    while len(stk):
        xs = stk.pop()
        i = len(xs)
        if i == V:
            yield xs
            continue
        v = xs[-1] & M1
        for x in vec[i]:
            if x >> 2 == v:
                stk.append(xs + [x])

def renderPath(K, xs):
    if len(xs) == 0:
        return ''
    r = [render(K, xs[0])]
    for x in xs[1:]:
        r.append("ACGT"[x&3])
    return ''.join(r)

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    verbose = opts['-v']

    K = int(opts['-k'])

    J = [min(j, K+1-j) for j in range(K+2)]
    W = [1 - math.exp(math.log(0.75)*j) for j in J]
    w0 = sum(W)
    W = [w/w0 for w in W]

    d = "."
    if opts['-g']:
        d = opts['-g']
    sf = SequenceFactory(d)

    if opts['-X']:
        with openFile(opts['<RefGene>']) as f:
            ref = readRefGene(f)

        res = {}
        for g in sorted(ref.keys()):
            seq = sf[g]
            for nm in sorted(ref[g].keys()):
                if nm not in res:
                    res[nm] = []
                for (s,b,e) in sorted(ref[g][nm]):
                    v = seq[b:e].upper()
                    x = (g, s, b, e, v)
                    res[nm].append(x)
        yaml.safe_dump(res, sys.stdout)
        return


    S = 1
    s = 0
    if opts['-S']:
        ss = opts['-S'].split('/')
        assert len(ss) == 2
        s = int(ss[0]) - 1
        S = int(ss[1])
        assert S > 0
        assert s < S

    if verbose:
        print >> sys.stderr, 'processing slice %d/%d' % (s+1,S)

    with openFile(opts['<index>']) as f:
        ref = yaml.load(f)
        for g in ref.keys():
            exons = []
            for t in ref[g]:
                exons.append(tuple(t))
            ref[g] = exons
        keep = sorted(ref.keys())[s::S]
        gx = sorted(set(ref.keys()) - set(keep))
        for g in gx:
            del ref[g]

    idx = {}
    for g in ref.keys():
        i = 0
        for t in ref[g]:
            ex = t[4]
            xs = kmersWithPosList(K, ex)
            if len(xs) == 0:
                print >> sys.stderr, 'warning: %s/%d length %d < K, has no k-mers' % (g, i, len(ex))
                continue
            for (x,p) in xs:
                p -= 1
                k = (g,i,p)
                if x not in idx:
                    idx[x] = []
                idx[x].append(k)
            i += 1

    rn = 0

    acc = {}
    grp = {}
    for itm in reads(opts['<input>'], K=K, paired=True, reads=False, kmers=True, separate=True, verbose=verbose):
        rn += 1
        lhs = itm.kmers[0]
        lhsFwd = lhs[0]
        lhsRev = lhs[1]
        rhs = itm.kmers[1]
        rhsFwd = rhs[0]
        rhsRev = rhs[1]

        idxHits0 = posHits(idx, lhsFwd, rhsRev)
        idxHits1 = posHits(idx, rhsFwd, lhsRev)

        if len(idxHits0) + len(idxHits1) == 0:
            continue

        gs = sorted(set([g for (g,i) in idxHits0.keys() + idxHits1.keys()]))
        for g in gs:
            if g not in acc:
                acc[g] = {}
            for x in lhsFwd + rhsRev:
                if x not in acc[g]:
                    acc[g][x] = 0
                acc[g][x] += 1

            k = tuple(sorted([(h,i) for (h,i) in idxHits0.keys() if h == g]))
            if len(k) > 0:
                if k not in grp:
                    grp[k] = 0
                grp[k] += 1

            k = tuple(sorted([(h,i) for (h,i) in idxHits1.keys() if h == g]))
            if len(k) > 0:
                if k not in grp:
                    grp[k] = 0
                grp[k] += 1

    if False:
        for g in acc.keys():
            hist = {}
            for v in acc[g].itervalues():
                if v not in hist:
                    hist[v] = 0
                hist[v] += 1
            for (c,f) in sorted(hist.items()):
                print '%s\t%d\t%d' % (g, c, f)
        return

    grpNum = 0
    for k in sorted(grp.keys()):
        grpNum += 1
        n = grp[k]
        for (g,i) in k:
            l = len(ref[g][i][4])
            print '%d\t%d\t%s\t%d\t%d' % (grpNum, n, g, i, l)

    for g in sorted(acc.keys()):
        rhsExtra = {}
        for x in sorted(rhsIdx.keys()):
            ys = []
            for y in ham1(K, x):
                if y not in acc[g]:
                    continue
                xc = acc[g].get(x, 0)
                yc = acc[g][y]
                if yc >= xc / 4:
                    if y not in rhsExtra:
                        rhsExtra[y] = set([])
                    rhsExtra[y].add(x)

        lhsExtra = {}
        for x in sorted(lhsIdx.keys()):
            ys = []
            for y in ham1(K, x):
                if y not in acc[g]:
                    continue
                xc = acc[g].get(x, 0)
                yc = acc[g][y]
                if yc >= xc / 4:
                    if y not in lhsExtra:
                        lhsExtra[y] = set([])
                    lhsExtra[y].add(x)

        vec = bridge(K, set(rhsIdx.keys() + rhsExtra.keys()), set(lhsIdx.keys() + lhsExtra.keys()), acc[g])
        res = {}
        for pth in paths(K, vec):
            fm = pth[0]
            if fm in rhsIdx:
                rr = rhsIdx[fm]
            else:
                rr = []
                for y in rhsExtra[fm]:
                    rr += rhsIdx[y]
            rr = sorted(set(rr))

            to = pth[-1]
            if to in lhsIdx:
                ll = lhsIdx[to]
            else:
                ll = []
                for y in lhsExtra[to]:
                    ll += lhsIdx[y]
            ll = sorted(set(ll))

            for r in rr:
                for l in ll:
                    kk = (r,l)
                    if kk not in res:
                        res[kk] = set([])
                    res[kk].add(tuple(pth))

        for kk in sorted(res.keys()):
            (r,l) = kk
            for pth in sorted(res[kk]):
                (hf,jf) = r
                (ht,jt) = l
                c = 0.0
                for i in range(K+1):
                    c += W[i]*acc[g][pth[i]]
                print hf, jf, ht, jt, renderPath(K, pth), c

        sys.stdout.flush()

if __name__ == '__main__':
    main(sys.argv[1:])
