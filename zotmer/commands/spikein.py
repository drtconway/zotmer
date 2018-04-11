"""
zot spikein - spike in variants

Usage:
    zot spikein [options] <output-prefix> [<variant>...]

Options:
    -b BEDFILE      restrict sampling of reads to regions in the BED file
    -d DEVIATION    standard deviation of insert sizes (proportion of mean) [default: 0.1]
    -e RATE         error rate [default: 0.005]
    -f FILE         read variants from a file
    -g DIR          directory containing reference sequences [default: .]
    -I SIZE         insert size [default: 300]
    -l FILE         log file for introduced mutations
    -L LENGTH       read length [default: 100]
    -m FILE         read a file containing <probability, variant> pairs, one per line to introduce stochastically
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
import os.path
import random
import sys

import docopt
from tqdm import tqdm

from pykmer.file import openFile, readFasta, readFastq
from pykmer.misc import uniq
from zotmer.library.hgvs import hg19ToRefSeq, makeHGVS, refSeq2Hg19, Substitution
from zotmer.library.rope import rope

def normalizeAccession(acc):
    if acc in refSeq2Hg19:
        acc = refSeq2Hg19[acc]
    if not acc.startswith('chr'):
        acc = 'chr' + acc
    return acc

class SequenceFactory(object):
    def __init__(self, home):
        self.home = home
        self.prevAcc = None
        self.prevSeq = None

    def __getitem__(self, acc):
        if acc != self.prevAcc:
            acc = normalizeAccession(acc)
            pth = self.home + '/' + acc + '.fa'
            if not os.path.exists(pth):
                pth = pth + '.gz'
            with openFile(pth) as f:
                for (nm,seq) in readFasta(f):
                    self.prevAcc = acc
                    self.prevSeq = seq
                    break
        return self.prevSeq

class MultiGen(object):
    def __init__(self, pcs):
        self.cum = []
        t = 0.0
        for (p,v) in pcs:
            t += p
            self.cum.append((t, v))

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
        self.p = p
        self.l1mp = None
        if p > 0:
            self.l1mp = math.log1p(-p)

    def __call__(self):
        if self.p == 0:
            return 0
        return int(math.log(random.random()) / self.l1mp)

    def make(self, lo, hi):
        if self.p == 0:
            return []
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
        ch = normalizeAccession(t[0])
        s = int(t[1])
        e = int(t[2])
        n = None
        if len(t) > 3:
            n = t[3]
        v = (s, e, n)
        if ch not in res:
            res[ch] = []
        res[ch].append(v)
    for ch in res.keys():
        res[ch].sort()
    return res

def applyBackgroundVariants(ch, given, popVars):
    if len(popVars) == 0:
        return given
    if ch not in popVars:
        return given

    extra = []
    for (m,p) in popVars[ch]:
        if p < random.random():
            continue
        if len(extra) and extra[-1].overlaps(m):
            continue
        ok = True
        for v in given:
            if v.overlaps(m):
                ok = False
                break
        if ok:
            extra.append(m)
    allOfThem = given + extra
    allOfThem.sort()
    return allOfThem

def between(st, en, v):
    """test if a variant lies within a range"""
    r = v.range()
    return all([st <= r[0], r[0] <= en, st <= r[1], r[1] <= en])

def overlaps(st, en, v):
    """test if a variant lies within a range"""
    r = v.range()
    if st <= r[0] and r[0] <= en:
        return True
    if st <= r[1] and r[1] <= en:
        return True
    if r[0] <= st and st <= r[1]:
        return True
    if r[0] <= en and en <= r[1]:
        return True
    return False

class ReadMaker(object):
    def __init__(self, **kwargs):
        # Chromosome and Region of interest
        #
        self.chrom = kwargs['chrom']
        self.zone = kwargs['zone']

        # Read Length
        #
        self.L = kwargs['L']

        # Insert Size
        #
        self.I = kwargs['I']

        # Fragment relative std-dev (as a faction of insert size)
        #
        self.D = kwargs['D']

        # List of *chromosome* variants
        #
        self.variants = kwargs['variants']

        # Fasta or Fastq?
        #
        self.fasta = kwargs['fasta']

        self.seqFac = kwargs['sf']
        self.egen = kwargs['egen']

    def prepareAllele(self):
        (st, en, nm) = self.zone

        # Select the variants that lie within this range.
        #
        ws = [v for v in self.variants if overlaps(st, en, v)]
        ws.sort(reverse=True)

        seq = self.seqFac[self.chrom]
        G = len(seq)

        s0 = st
        e0 = en
        seq = rope.atom(seq)

        for v in ws:
            r = v.range()
            l = r[1] - r[0]
            e0 -= l
            e0 += v.size()
            seq = v.apply(seq)

        self.seq = seq
        self.sourceRange = (s0, e0)
        self.quals = 'I' * self.L

    def makeReadFromZone(self):
        G = len(self.seq)

        (s0, e0) = self.sourceRange

        u = random.randint(s0, e0)

        fl = int(random.gauss(self.I, self.I*self.D))
        fl = max(fl, self.L + self.L // 2)

        if u - fl < 0:
            u = fl
        if u + fl > G:
            u = G - fl

        assert 0 <= u and u + fl <= G

        if random.random() < 0.5:
            rp1 = u
            rp2 = u + fl - self.L
        else:
            rp1 = u - fl
            rp2 = u - self.L

        r1 = self.seq[rp1:rp1+self.L]
        r2 = self.seq[rp2:rp2+self.L]

        if random.random() < 0.5:
            s = '+'
            (r1, r2) = (r1, revComp(r2))
        else:
            s = '-'
            (r1, r2) = (r2, revComp(r1))
            (rp1,rp2) = (rp2,rp1)

        (r1, e1) = mutate(self.egen, r1)
        (r2, e2) = mutate(self.egen, r2)

        ch = self.chrom
        (st, en, nm) = self.zone
        out1 = []
        out2 = []
        if self.fasta:
            out1 += ['>' + ' '.join([ch, str(st), str(en), nm, str(rp1), str(fl), s] + e1)]
            out1 += [r1]

            out2 += ['>' + ' '.join([ch, str(st), str(en), nm, str(rp2), str(fl), s] + e2)]
            out2 += [r2]
        else:
            out1 += ['@%s:%d-%d %s %d %d %s %s' % (ch, st, en, nm, rp1, fl, s, ';'.join(e1))]
            out1 += [r1]
            out1 += ['+']
            out1 += [self.quals]
            out2 += ['@%s:%d-%d %s %d %d %s %s' % (ch, st, en, nm, rp2, fl, s, ';'.join(e2))]
            out2 += [r2]
            out2 += ['+']
            out2 += [self.quals]
        return ('\n'.join(out1), '\n'.join(out2))

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
        for ch in refSeq2Hg19.values():
            ch = normalizeAccession(ch)
            s = sf[ch]
            v = (1, len(s), ch)
            if ch not in zones:
                zones[ch] = []
            zones[ch].append(v)

    popVars = {}
    if opts['-m'] is not None:
        with openFile(opts['-m']) as f:
            for l in f:
                t = l.split()
                p = float(t[0])
                v = makeHGVS(t[1], sf)
                a = normalizeAccession(v.accession())
                if a not in popVars:
                    popVars[a] = []
                popVars[a].append((v,p))
    for ch in popVars.keys():
        popVars[ch].sort()

    t = 0
    chMax = 0
    for ch in zones.keys():
        chMax = max(chMax, len(ch))
        for zone in zones[ch]:
            (s, e, _n) = zone
            l = e - s + 1
            t += l

    if opts['-v']:
        print >> sys.stderr, 'mean coverage = %g' % (float(N*L)/float(t),)

    zoneCounts = {}
    zoneProbs = []
    for ch in zones.keys():
        zoneCounts[ch] = {}
        for zone in zones[ch]:
            zoneCounts[ch][zone] = 0
            (s, e, _n) = zone
            l = e - s + 1
            zoneProbs.append((float(l)/float(t), (ch, zone)))
    zgen = MultiGen(zoneProbs)
    for n in xrange(N):
        (ch,z) = zgen.gen()
        if z not in zoneCounts[ch]:
            zoneCounts[ch][z] = 0
        zoneCounts[ch][z] += 1

    vStrs = opts['<variant>']
    if opts['-f'] is not None:
        with openFile(opts['-f']) as f:
            for l in f:
                s = l.strip()
                vStrs.append(s)

    allVars = {}
    for s in vStrs:
        v = makeHGVS(s, sf)
        if v is None:
            print >> sys.stderr, 'unable to parse variant: %s', (s, )
            continue
        if v.anonymous():
            n = v.size()
            seq = ''.join([random.choice(['A', 'C', 'G', 'T']) for i in range(n)])
            v.setSequence(seq)
        a = normalizeAccession(v.accession())
        if a not in allVars:
            allVars[a] = []
        allVars[a].append(v)

    numOverlaps = 0
    for xs in allVars.values():
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

    logfile = None
    if opts['-l']:
        logfile = open(opts['-l'], 'w')

    pfx = opts['<output-prefix>']
    sfx = ''
    if opts['-z']:
        sfx = '.gz'
    with openFile(pfx + '_1.fastq' + sfx, 'w') as out1, openFile(pfx + '_2.fastq' + sfx, 'w') as out2:
        for ch in zones.keys():
            if prog is not None:
                prog.set_description(ch.ljust(chMax, ' '))
                prog.update(0)
            chVars = []
            if ch in allVars:
                chVars = allVars[ch]

            wtVars = applyBackgroundVariants(ch, [], popVars)
            mutVars = applyBackgroundVariants(ch, chVars, popVars)

            for zone in zones[ch]:
                wtMaker = ReadMaker(chrom=ch, zone=zone, L=L, I=I, D=D, variants=wtVars, fasta=fasta, egen=egen, sf=sf)
                wtMaker.prepareAllele()

                mutMaker = ReadMaker(chrom=ch, zone=zone, L=L, I=I, D=D, variants=mutVars, fasta=fasta, egen=egen, sf=sf)
                mutMaker.prepareAllele()

                for i in xrange(zoneCounts[ch][zone]) :
                    if prog is not None:
                        prog.update(1)
                    u = random.random()
                    if u > V:
                        (rd1,rd2) = wtMaker.makeReadFromZone()
                    else:
                        (rd1,rd2) = mutMaker.makeReadFromZone()
                    print >> out1, rd1
                    print >> out2, rd2

    if prog is not None:
        prog.__exit__(None, None, None)

if __name__ == '__main__':
    main(sys.argv[1:])

