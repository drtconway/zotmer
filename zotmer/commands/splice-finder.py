"""
zot splice-finder - find splice junctions from known exons

Usage:
    zot splice-finder [options] <reference> <input>...

Options:
    -k K            value of k to use [default: 25]
    -C              show read counts only
    -F              feature mode
    -g PATH         directory of FASTQ reference sequences
    -T FASTA        transcript to capture against.
    -v              produce verbose output
"""

import os.path
import sys

import docopt
import yaml

from pykmer.basics import ham, kmersList, kmersWithPosLists, murmer, rc, render
from pykmer.file import openFile, readFasta, readFastq
from zotmer.library.hgvs import hg19ToRefSeq, refSeq2Hg19
from zotmer.library.index import kmerPosIndex
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

def simpleHits(idx, xps):
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
    return Z

def updateLinkage(Z, links):
    Y = []
    for (z, q, p) in sorted(Z):
        if len(Y) > 0 and Y[-1][0] == z and Y[-1][2] + 1 == q and Y[-1][4] + 1 == p:
            Y[-1] = (Y[-1][0], Y[-1][1], q, Y[-1][3], p)
        else:
            Y.append((z, q, q, p, p))

    # For now, let's assume any Y != 2 is spurious.
    if len(Y) != 2:
        return

    for i in range(1, len(Y)):
        (z0, q0, p0, r0, s0) = Y[i-1]
        (z1, q1, p1, r1, s1) = Y[i]
        if z0 != z1 or q1 - p0 != r1 - s0:
            lnk = ((z0,p0), (z1,q1))
            if lnk not in links:
                links[lnk] = 0
            links[lnk] += 1
    return

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

def accumulateHits(idx, K, xps, acc):
    hs = hits(idx, K, xps)
    for z in hs.keys():
        for r in hs[z].keys():
            for (x,p) in xps:
                p -= 1
                p += r
                if z not in acc:
                    acc[z] = {}
                if p not in acc[z]:
                    acc[z][p] = {}
                if x not in acc[z][p]:
                    acc[z][p][x] = 0
                acc[z][p][x] += 1

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    verbose = opts['-v']

    K = int(opts['-k'])

    D = 3

    Q = K // 2

    d = "."
    if opts['-g']:
        d = opts['-g']
    sf = SequenceFactory(d)

    #with openFile(opts['<exons>']) as f:
    #    R = readBED(f)

    idx = kmerPosIndex(K)
    with openFile(opts['<reference>']) as f:
        for (nm, seq) in readFasta(f):
            idx.addSeq(nm, seq)

    X = {}
    for itm in reads(opts['<input>'], K=K, paired=True, reads=True, kmers=False, verbose=verbose):
        rdL = itm.reads[0]

        rdR = itm.reads[1]

        idx.pairHits(rdL[1], rdR[1], X)

    for z in sorted(X.keys()):
        for p in sorted(X[z].keys()):
            y = idx.atPos(z, p)
            if y is None:
                continue
            xs = set(X[z][p].keys()) - set([y])
            for x in xs:
                if ham(x,y) < 4:
                    continue
                xc = X[z][p][x]
                if xc <= 1:
                    continue
                for (u,q) in idx[x]:
                    #if abs(q - p) > 10000:
                    #    continue
                    print '%s\t%d\t%s\t%s\t%d\t%s' % (z, 1+p, render(K, y), u, 1+q, render(K, x))

if __name__ == '__main__':
    main(sys.argv[1:])

