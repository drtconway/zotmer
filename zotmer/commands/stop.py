"""
zot stop - search for stop codons

Usage:
    zot stop [options] <input>...

Options:
    -k K            k-mer size [default: 27]
    -r FILE         filter candidates against the given transcripts
    -v              produce verbose output
"""

import math
import random
import re
import sys

import docopt
import yaml

from zotmer.library.basics import ham, kmer, kmersList, kmersWithPosList, render
from zotmer.library.file import openFile, readFasta, readFastq
from zotmer.library.prot import nucToAa3

stopStrs = ['TAA', 'TAG', 'TGA']
stops = [kmer(s) for s in stopStrs]

def kmersWithPosLists(K, seq):
    fwd = []
    rev = []
    L = len(seq) - K + 1
    for (x,p) in kmersWithPosList(K, seq, True):
        if p > 0:
            fwd.append((x, p - 1))
        else:
            rev.append((x, L + p))
    return (fwd, rev)

def codons(o, seq):
    L = len(seq)
    res = [' ' * o]
    for i in xrange(o, L, 3):
        if i + 3 <= L:
            c = seq[i:i+3]
            x = kmer(c)
            res.append(nucToAa3(x))
    return ''.join(res)

rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def proj(s, seq):
    if s == 0:
        return seq
    res = []
    for c in seq:
        res.append(rcDict.get(c, 'N'))
    return ''.join(res[::-1])

def nearest(K, xs, x):
    d = K+1
    z = x
    for y in xs:
        d0 = ham(x, y)
        if d0 < d:
            d = d0
            z = y
    return (d, z)

def nearest3(K, xs, x):
    d = K+1
    z = x
    for y in xs:
        d0 = ham(x >> 3, y >> 3)
        if d0 < d:
            d = d0
            z = y
    return (d, z >> 3)

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    random.seed(17)

    K = int(opts['-k'])
    S = 2*(K-3)

    frameAnchors = {}
    knownStops = {}
    sequences = {}
    seqKmers = {}

    if opts['-r']:
        with openFile(opts['-r']) as f:
            for (nm,seq) in readFasta(f):
                sequences[nm] = seq
                # trim polyA tails
                seq = re.sub('AAAAAA*$', '', seq)
                seqKmers[nm] = set([])
                for (x,p1) in kmersWithPosList(K, seq, False):
                    seqKmers[nm].add(x)
                    p = p1 - 1
                    w = p % 3
                    if x not in frameAnchors:
                        frameAnchors[x] = set([])
                    frameAnchors[x].add((nm,p))
                    y = x & 63
                    if w == 0 and y in stops:
                        if x not in knownStops:
                            knownStops[x] = set([])
                        knownStops[x].add(nm)

    rn = 0
    res = {}
    for fn in opts['<input>']:
        with openFile(fn) as f:
            for rd in readFastq(f):
                L = len(rd[1])
                rn += 1
                fwdAndRev = kmersWithPosLists(K, rd[1])
                frames = {}
                possibleStops = {}
                for i in range(2):
                    #print i, sorted([p for (x,p) in fwdAndRev[i]])
                    for (x,p) in fwdAndRev[i]:
                        if x in frameAnchors:
                            for (nm,q) in frameAnchors[x]:
                                o = (q - p)
                                k = (nm, o, i)
                                frames[k] = 1 + frames.get(k, 0)
                if len(frames) == 0:
                    continue
                n = sum(frames.values())
                probs = []
                for ((nm, off, strnd), cnt) in sorted(frames.items()):
                    probs.append((float(cnt)/float(n), cnt, off, strnd, nm))
                v = random.random()
                for (pv, cnt, off, strnd, nm) in probs:
                    if v < pv:
                        #print rd[1]
                        #print proj(strnd, sequences[nm][off:off+len(rd[1])])
                        #print codons(off % 3, rd[1]), off
                        for (x,p) in fwdAndRev[strnd]:
                            if (p + off + K - 3) % 3 == 0 and (x & 63) in stops:
                                if nm not in res:
                                    res[nm] = {}
                                if x not in res[nm]:
                                    res[nm][x] = 0
                                res[nm][x] += 1
                        break
                    v -= pv
    for (nm,stps) in res.iteritems():
        for (x,c) in stps.iteritems():
            (d,y) = nearest3(K, seqKmers[nm], x)
            if x in knownStops:
                k = 'known'
            else:
                k = 'novel'
            print '%s\t%s\t%d\t%d\t%s\t%s' % (k, render(K, x), c, d, render(K, y), nm)

if __name__ == '__main__':
    main(sys.argv[1:])
