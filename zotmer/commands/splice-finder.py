"""
zot splice-finder - find splice junctions from known exons

Usage:
    zot splice-finder [options] <exons> <input>...

Options:
    -k K            value of k to use [default: 25]
    -F              feature mode
    -g PATH         directory of FASTQ reference sequences
    -v              produce verbose output
"""

import os.path
import sys

import docopt
import yaml

from pykmer.basics import ham, kmersWithPosList, kmersWithPosLists, murmer, rc, render
from pykmer.file import openFile, readFasta, readFastq
from zotmer.library.hgvs import hg19ToRefSeq, refSeq2Hg19
from zotmer.library.reads import reads

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
        d = None
        if len(t) > 4:
            d = t[4]
        v = (s, e, n, d)
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

def updateLinkage(idx, xps, links):
    Z = []
    for (x,p) in xps:
        zqs = idx[x]
        if len(zqs) == 0:
            continue
        if len(zqs) > 1:
            # damn. a tricky case - let's just skip it for now
            # ideally we'd see if there's one that's consistent
            # with the other positions, since that's most likely,
            # and ignore the rest, or pick one somehow.
            continue
        Z.append((zqs[0][0], zqs[0][1], p))
    Z.sort()

    for i in range(1, len(Z)):
        (z0, q0, p0) = Z[i-1]
        (z1, q1, p1) = Z[i]
        if z0 != z1 or q1 - q0 != p1 - p0:
            lnk = ((z0,q0), (z1,q1))
            if lnk not in links:
                links[lnk] = 0
            links[lnk] += 1

def hits(idx, K, xps):
    loc = {}
    for (x,p) in xps:
        p -= 1
        if x not in idx:
            continue
        for (z,q) in idx[x]:
            r = q - p
            if z not in loc:
                loc[z] = {}
            if r not in loc[z]:
                loc[z][r] = 0
            loc[z][r] += 1
    return loc

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    verbose = opts['-v']

    K = int(opts['-k'])

    feature = opts['-F']

    D = 3

    Q = K // 2

    d = "."
    if opts['-g']:
        d = opts['-g']
    sf = SequenceFactory(d)

    with openFile(opts['<exons>']) as f:
        R = readBED(f)

    refIdx = {}
    refTbl = {}
    zoneIdx = {}
    for (acc, zones) in R.items():
        if verbose:
            print >> sys.stderr, 'reading %s' % (acc, )
        accSeq = sf[acc]
        for (s, e, nm, d) in zones:
            refTbl[nm] = {}
            zoneIdx[nm] = (acc, s, e, d)
            seq = accSeq[s-1:e]
            xps = kmersWithPosList(K, seq, False)
            for (x,p) in xps:
                p -= 1
                if x not in refIdx:
                    refIdx[x] = set([])
                refIdx[x].add((nm,p))
                refTbl[nm][p] = x

    counts = {}
    X = {}
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

        if feature:
            for (hs,xps) in [(fwdLHits, fwdL), (revLHits, revL), (fwdRHits, fwdR), (revRHits, revR)]:
                for z in sorted(hs.keys()):
                    for r in sorted(hs[z].keys()):
                        for (x,p) in xps:
                            q = r + p
                            if q % Q != 0:
                                continue
                            if q not in refTbl[z]:
                                continue
                            y = refTbl[z][q]
                            if ham(x, y) > D:
                                continue
                            if z not in X:
                                X[z] = {}
                            if q not in X[z]:
                                X[z][q] = {}
                            if x not in X[z][q]:
                                X[z][q][x] = 0
                            X[z][q][x] += 1
            continue

        if len(fwdLHits) > 0 or len(revRHits) > 0:
            k = '&'.join(sorted(set(fwdLHits.keys()) | set(revRHits.keys())))
            if k not in counts:
                counts[k] = 0
            counts[k] += 1
            if k not in X:
                X[k] = {}
            Xk = X[k]
            for (x,p) in fwdL:
                if x not in Xk:
                    Xk[x] = 0
                Xk[x] += 1
            for (x,p) in revR:
                if x not in Xk:
                    Xk[x] = 0
                Xk[x] += 1

        if len(fwdRHits) > 0 or len(revLHits) > 0:
            k = '&'.join(sorted(set(fwdRHits.keys()) | set(revLHits.keys())))
            if k not in counts:
                counts[k] = 0
            counts[k] += 1
            if k not in X:
                X[k] = {}
            Xk = X[k]
            for (x,p) in fwdR:
                if x not in Xk:
                    Xk[x] = 0
                Xk[x] += 1
            for (x,p) in revL:
                if x not in Xk:
                    Xk[x] = 0
                Xk[x] += 1

    if feature:
        for z in sorted(X.keys()):
            for p in sorted(X[z].keys()):
                for (x,c) in sorted(X[z][p].items()):
                    print '%s\t%d\t%s\t%d' % (z, p, render(K, x), c)
        return

    for (k,c) in counts.items():
        for (x,n) in X[k].items():
            print '%s\t%d\t%d\t%s' % (render(K, x), n, c, k)

if __name__ == '__main__':
    main(sys.argv[1:])

