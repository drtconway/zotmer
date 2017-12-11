"""
zot spikein - spike in variants

Usage:
    zot spikein [options] <output-prefix> [<variant>...]

Options:
    -b BEDFILE      restrict sampling of reads to regions in the BED file
    -d DEVIATION    standard deviation of insert sizes (proportion of mean) [default: 0.1]
    -e RATE         error rate [default: 0.01]
    -g DIR          directory containing reference sequences [default: .]
    -I SIZE         insert size [default: 300]
    -L LENGTH       read length [default: 100]
    -N NUMBER       number of read pairs [default: 1000000]
    -S SEED         random number generator seed [default: 23]
    -z              gzip the output
    -V VAF          variant allele frequency [default: 1.0]

Variants are spiked in right-to-left from the command line or reverse order when read from a file.
There is no compensation for position in the application of the variants, so they should be specified
in position order, and overlapping variants won't work.
"""

import math
import random
import sys

import docopt

from pykmer.file import openFile, readFasta, readFastq
from zotmer.library.hgvs import hg19ToRefSeq, makeHGVS, refSeq2Hg19
from zotmer.library.rope import rope

class SequenceFactory(object):
    def __init__(self, home):
        self.home = home
        self.prevAcc = None
        self.prevSeq = None

    def __getitem__(self, acc):
        if acc != self.prevAcc:
            if acc not in refSeq2Hg19:
                print >> sys.stderr, "accession %s not available." % (acc)
            assert acc in refSeq2Hg19
            h = refSeq2Hg19[acc]

            with openFile(self.home + "/" + h + ".fa.gz") as f:
                for (nm,seq) in readFasta(f):
                    self.prevAcc = acc
                    self.prevSeq = seq
                    break
        return self.prevSeq

class MultiGen(object):
    def __init__(self, pcs):
        self.cum = []
        t = 0.0
        for (p,c) in pcs:
            t += p
            self.cum.append((t, c))

    def gen(self):
        u = random.random()
        l = 0
        h = len(self.cum)
        while l < h:
            m = (l + h) // 2
            if u < self.cum[m][0]:
                h = m
            else:
                l = m + 1
        return self.cum[l][1]

rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def revComp(s):
    r = []
    for c in s:
        r.append(rcDict.get(c, c))
    return ''.join(r[::-1])

altDict = {}
altDict['A'] = ['C', 'G', 'T']
altDict['C'] = ['A', 'G', 'T']
altDict['G'] = ['A', 'C', 'T']
altDict['T'] = ['A', 'C', 'G']

bases = ['A', 'C', 'G', 'T']

def mutate(E, seq):
    r = []
    e = []
    i = 0
    for c in seq:
        c = c.upper()
        if random.random() < E:
            d = random.choice(altDict.get(c, bases))
            e.append('%d%s>%s' % (i, c, d))
            c = d
        r.append(c)
        i += 1
    return (''.join(r), e)

def readBED(f):
    res = {}
    first = True
    for l in f:
        t = l.split()
        if first:
            first = False
            if t[0] == 'track' or t[0] =='browser':
                continue
        c = hg19ToRefSeq[t[0]]
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

def liftover(vs, p):
    for v in vs[::-1]:
        p = v.liftover(p)
    return p

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    sf = SequenceFactory(opts['-g'])

    I = int(opts['-I'])
    D = float(opts['-d'])
    E = float(opts['-e'])
    L = int(opts['-L'])
    N = int(opts['-N'])
    S = int(opts['-S'])
    V = float(opts['-V'])

    random.seed(S)

    if opts['-b']:
        zones = readBED(openFile(opts['-b']))
    else:
        zones = {}
        for c in refSeq2Hg19.keys():
            s = sf[c]
            v = (1, len(s), refSeq2Hg19[c])

    zps = []
    t = 0
    for c in zones.keys():
        for i in range(len(zones[c])):
            (s, e, _n) = zones[c][i]
            l = e - s + 1
            zps.append((l, (c,i)))
            t += l
    zps = [(float(l)/float(t), c) for (l,c) in zps]
    zgen = MultiGen(zps)
    zcs = dict([(c,0) for (p,c) in zps])
    for n in xrange(N):
        c = zgen.gen()
        zcs[c] += 1

    vs = {}
    for s in opts['<variant>']:
        v = makeHGVS(s, sf)
        if v is None:
            print >> sys.stderr, 'unable to parse variant: %s', (s, )
            continue
        if v.anonymous():
            n = v.size()
            seq = ''.join([random.choice(['A', 'C', 'G', 'T']) for i in range(n)])
            v.setSequence(seq)
        a = v.accession()
        if a not in vs:
            vs[a] = []
        vs[a].append(v)

    pfx = opts['<output-prefix>']
    sfx = ''
    if opts['-z']:
        sfx = '.gz'
    with openFile(pfx + '_1.fastq' + sfx, 'w') as out1, openFile(pfx + '_2.fastq' + sfx, 'w') as out2:
        wtSeq = None
        currC = None
        currSeq = None
        for (c,i) in sorted(zcs.keys()):
            if c != currC:
                print >> sys.stderr, 'loading %s' % (c, )
                currC = c
                currSeq = rope.atom(sf[c])
                for v in vs.get(c, [])[::-1]:
                    print >> sys.stderr, 'applying %s' % (str(v), )
                    currSeq = v.apply(currSeq)
                wtSeq = rope.atom(sf[c])

            n = zcs[(c,i)]
            (zb, ze, nm) = zones[c][i]
            for j in xrange(n):
                u0 = random.randint(zb, ze)
                if random.random() < V:
                    u = liftover(vs.get(c, []), u0)
                    ss = currSeq
                else:
                    u = u0
                    ss = wtSeq

                fl = int(random.gauss(I, I*D))
                fl = max(fl, L + L // 2)

                if random.random() < 0.5:
                    rp1 = u
                    rp2 = u + fl - L
                else:
                    rp1 = u - fl
                    rp2 = u - L

                r1 = ss[rp1:rp1+L]
                r2 = ss[rp2:rp2+L]

                if random.random() < 0.5:
                    s = '+'
                    r2 = revComp(r2)
                else:
                    s = '-'
                    tmp = r2
                    r2 = revComp(r1)
                    r1 = tmp

                (r1, e1) = mutate(E, r1)
                (r2, e2) = mutate(E, r2)

                print >> out1, '>' + ' '.join([c, str(zb), str(ze), nm, str(rp1), str(fl), s] + e1)
                print >> out1, r1
                #print >> out1, '\t'.join([c, str(zb), str(ze), nm, str(rp1), str(fl), s, r1] + e1)

                print >> out2, '>' + ' '.join([c, str(zb), str(ze), nm, str(rp2), str(fl), s] + e2)
                print >> out2, r2
                #print >> out2, '\t'.join([c, str(zb), str(ze), nm, str(rp2), str(fl), s, r2] + e2)

if __name__ == '__main__':
    main(sys.argv[1:])

