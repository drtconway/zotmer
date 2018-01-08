"""
zot spikein - spike in variants

Usage:
    zot spikein [options] <output-prefix> [<variant>...]

Options:
    -b BEDFILE      restrict sampling of reads to regions in the BED file
    -d DEVIATION    standard deviation of insert sizes (proportion of mean) [default: 0.1]
    -e RATE         error rate [default: 0.005]
    -g DIR          directory containing reference sequences [default: .]
    -I SIZE         insert size [default: 300]
    -l FILE         log file for introduced mutations
    -L LENGTH       read length [default: 100]
    -M PROB         introduce extra substitution mutations around the variants with the given per-base probability
    -N NUMBER       number of read pairs [default: 1000000]
    -S SEED         random number generator seed
    -z              gzip the output
    -v              produce verbose progress messages
    -V VAF          variant allele frequency [default: 1.0]

Variants are spiked in right-to-left from the command line or reverse order when read from a file.
There is no compensation for position in the application of the variants, so they should be specified
in position order, and overlapping variants won't work.
"""

import math
import random
import sys

import docopt
from tqdm import tqdm

from pykmer.file import openFile, readFasta, readFastq
from pykmer.misc import uniq
from zotmer.library.hgvs import hg19ToRefSeq, makeHGVS, refSeq2Hg19, Substitution
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

class GeomVarSource(object):
    def __init__(self, p):
        self.l1mp = math.log1p(-p)

    def __call__(self):
        return int(math.log(random.random()) / self.l1mp)

    def make(self, lo, hi):
        r = []
        j = lo
        i = 0
        while True:
            j += i + self()
            if j >= hi:
                break
            i = 1
            r.append(j)
        return r

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

def mutate(egen, seq):
    ps = egen.make(0, len(seq))

    seq = seq.upper()

    if len(ps) == 0:
        return (seq, [])

    r = [seq[i] for i in range(len(seq))]
    e = []
    for p in ps:
        c = r[p]
        d = random.choice(altDict.get(c, bases))
        e.append('%d%s>%s' % (p, c, d))
        r[p] = d
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
    V = float(opts['-V'])

    M = None
    if opts['-M'] is not None:
        M = float(opts['-M'])
        # compute the 99% quantile
        W = int(math.log1p(-0.99)/math.log1p(-M))

    S = None
    if opts['-S'] is not None:
        S = int(opts['-S'])
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
    maxC = 0
    for c in zones.keys():
        maxC = max(maxC, len(c))
        for i in range(len(zones[c])):
            (s, e, _n) = zones[c][i]
            l = e - s + 1
            zps.append((l, (c,i)))
            t += l

    if opts['-v']:
        print >> sys.stderr, 'mean coverage = %g' % (float(N*L)/float(t),)

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

    numOverlaps = 0
    for xs in vs.values():
        xs.sort()
        for i in range(len(xs)):
            for j in range(i + 1, len(xs)):
                if xs[i].overlaps(xs[j]):
                    print >> sys.stderr, "variants overlap: %s <> %s" % (str(xs[i]), str(xs[j]))
                    numOverlaps += 1
    if numOverlaps > 0:
        sys.exit(1)

    prog = None
    if opts['-v']:
        prog = tqdm(total = N, unit='pairs')

    egen = GeomVarSource(E)

    fasta = False
    quals = 'I' * L

    logfile = None
    if opts['-l']:
        logfile = open(opts['-l'], 'w')

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
                if prog is not None:
                    prog.set_description(c.ljust(maxC, ' '))
                    prog.update(0)
                currC = c
                seq = sf[c]
                Z = len(seq)
                currSeq = rope.atom(seq)
                given = vs.get(c, [])
                if M is not None:
                    extra = []
                    mgen = GeomVarSource(M)
                    for v in given:
                        vr = v.range()
                        extra += mgen.make(max(0, vr[0] - W), vr[0])
                        extra += mgen.make(vr[1], min(vr[1] + W, Z))
                    extra.sort()
                    uniq(extra)
                    ok = []
                    for p in extra:
                        r = seq[p].upper()
                        w = Substitution(c, p, r, random.choice(altDict[r]))
                        hit = False
                        for v in given:
                            if v.overlaps(w):
                                hit = True
                                break
                        if not hit:
                            ok.append(w)
                    if logfile is not None:
                        for v in ok:
                            print >> logfile, str(v)
                    given += ok
                    given.sort()
                    if prog is not None and len(given) > 0:
                        prog.write("spiking in the following variants:")
                        for v in given:
                            prog.write('\t' + str(v))
                for v in given[::-1]:
                    currSeq = v.apply(currSeq)
                wtSeq = rope.atom(sf[c])

            n = zcs[(c,i)]
            (zb, ze, nm) = zones[c][i]
            if prog is not None:
                prog.update(n)
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

                (r1, e1) = mutate(egen, r1)
                (r2, e2) = mutate(egen, r2)

                #print >> out1, '\t'.join([c, str(zb), str(ze), nm, str(rp1), str(fl), s, r1] + e1)
                #print >> out2, '\t'.join([c, str(zb), str(ze), nm, str(rp2), str(fl), s, r2] + e2)
                if fasta:
                    print >> out1, '>' + ' '.join([c, str(zb), str(ze), nm, str(rp1), str(fl), s] + e1)
                    print >> out1, r1
                    print >> out2, '>' + ' '.join([c, str(zb), str(ze), nm, str(rp2), str(fl), s] + e2)
                    print >> out2, r2
                else:
                    print >> out1, '@%s:%d-%d %s %d %d %s %s' % (c, zb, ze, nm, rp1, fl, s, ';'.join(e1))
                    print >> out1, r1
                    print >> out1, '+'
                    print >> out1, quals
                    print >> out2, '@%s:%d-%d %s %d %d %s %s' % (c, zb, ze, nm, rp2, fl, s, ';'.join(e2))
                    print >> out2, r2
                    print >> out2, '+'
                    print >> out2, quals

    if prog is not None:
        prog.__exit__(None, None, None)

if __name__ == '__main__':
    main(sys.argv[1:])

