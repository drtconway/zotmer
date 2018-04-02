"""
zot hgvs-gen - generate random variants

Usage:
    zot hgvs-gen [options] <regions>

Options:
    -D LENGTH       mean length for deletion events [default: 1.5]
    -I LENGTH       mean length for insertion events [default: 1.5]
    -U LENGTH       mean length for duplication events [default: 3.0]
    -g PATH         directory of FASTQ reference sequences
    -N NUMBER       number of HGVSg variants to generate [default: 1]
    -S SEED         set a random number seed
    -t TYPES        comma separated list of variant types [default: sub,del,ins,delins,dup]
    -T              produce output in a tabular form
    -v              produce progress messages
    -V              produce verbose output

Variant Types

The following variant types are supported
    sub     substitution
    del     deletion
    ins     insertion
    delins  deletion/insertion
    dup     duplication
    ani     anonymous insertion
    and     anonymous deletion/insertion

If specified as prob:type, the type of variant is generated with a
relative probability of prob, with the default being 1.0.

"""

import math
import random
import sys

import docopt
from tqdm import tqdm

from pykmer.file import openFile, readFasta, readFastq
from zotmer.library import hgvs

def geomvar(p):
    return int(1+math.floor(math.log(random.random())/math.log1p(-p)))

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
        if c in hgvs.hg19ToRefSeq:
            c = hgvs.hg19ToRefSeq[c]
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
            if acc not in hgvs.refSeq2Hg19:
                print >> sys.stderr, "accession %s not available." % (acc)
            assert acc in hgvs.refSeq2Hg19
            h = hgvs.refSeq2Hg19[acc]

            with openFile(self.home + "/" + h + ".fa.gz") as f:
                for (nm,seq) in readFasta(f):
                    self.prevAcc = acc
                    self.prevSeq = seq
                    break
        return self.prevSeq

class MultiGen(object):
    def __init__(self, pcs):
        # Accumulate cumulative relative probabilities
        self.cum = []
        t = 0.0
        for (p,c) in pcs:
            t += p
            self.cum.append((t, c))

        # Normalise
        for i in range(len(self.cum)):
            (p, c) = self.cum[i]
            self.cum[i] = (p/t, c)

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

bases = ['A','C','G','T']
alts = {'A':['C','G','T'], 'C':['A','G','T'], 'G':['A','C','T'], 'T':['A','C','G']}

def genVar(c, seq, s, e, ts, ds):
    assert 0 <= s and s < len(seq)
    assert s <= e and e < len(seq)
    f = ts.gen()
    if f == 'sub':
        p = random.randint(s, e)
        r = seq[p].upper()
        a = random.choice(alts[r])
        return hgvs.Substitution(c, p, r, a)
    if f == 'del':
        p = random.randint(s, e)
        l = geomvar(ds['del'])
        q = min(p+l, e)
        return hgvs.Deletion(c, p, q)
    if f == 'ins':
        p = random.randint(s, e)
        m = geomvar(ds['ins'])
        a = ''.join([random.choice(bases) for j in range(m)])
        return hgvs.Insertion(c, p, a)
    if f == 'delins':
        p = random.randint(s, e)
        l = geomvar(ds['del'])
        q = min(p+l, e)
        m = geomvar(ds['ins'])
        a = ''.join([random.choice(bases) for j in range(m)])
        return hgvs.DeletionInsertion(c, p, q, a)
    if f == 'dup':
        p = random.randint(s, e)
        u = geomvar(ds['dup'])
        q = min(p+u, e)
        return hgvs.Duplication(c, p, q)
    if f == 'ani':
        p = random.randint(s, e)
        m = geomvar(ds['ins'])
        return hgvs.Anonymous(c, p, m)
    if f == 'and':
        p = random.randint(s, e)
        l = geomvar(ds['del'])
        q = min(p+l, e)
        m = geomvar(ds['ins'])
        return hgvs.AnonymousDelIns(c, p, q, m)
    print >> sys.stderr, "unknown variant type", f
    sys.exit(1)

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    d = "."
    if opts['-g']:
        d = opts['-g']
    sf = SequenceFactory(d)

    T = []
    for t0 in opts['-t'].split(','):
        t1 = t0.split(':')
        if len(t1) == 1:
            T.append((1.0, t1[0]))
        elif len(t1) == 2:
            T.append((float(t1[0]), t1[1]))
        else:
            print >> sys.stderr, "unexpected variant category descriptor:", t0
            sys.exit(1)
    Tgen = MultiGen(T)

    N = int(opts['-N'])

    D = float(opts['-D'])
    I = float(opts['-I'])
    U = float(opts['-U'])

    Ds = {}
    Ds['del'] = 1.0/D
    Ds['ins'] = 1.0/I
    Ds['dup'] = 1.0/U

    verbose = opts['-v']

    if opts['-S'] is not None:
        random.seed(int(opts['-S']))

    with openFile(opts['<regions>']) as f:
        zones = readBED(f)

    zps = []
    t = 0
    maxC = 0
    maxM = 0
    for c in zones.keys():
        maxC = max(maxC, len(c))
        for i in range(len(zones[c])):
            (s, e, m) = zones[c][i]
            maxM = max(maxM, len(m))
            l = e - s + 1
            zps.append((l, (c,i)))
            t += l
    zps = [(float(l)/float(t), c) for (l,c) in zps]
    zgen = MultiGen(zps)
    zcs = dict([(c, 0) for (p, c) in zps])
    for n in xrange(N):
        c = zgen.gen()
        zcs[c] += 1
    zcs = [(c,n) for (c,n) in zcs.items() if n > 0]
    zcs.sort()

    prog = None
    progFmt = None
    if verbose:
        prog = tqdm(total = N, unit='vars')
        progFmt = '%-' + str(maxC) + 's : %-' + str(maxM) + 's'
    curC = None
    curSeq = None
    prevM = None
    for ((c, i), n) in zcs:
        (s, e, m) = zones[c][i]
        if curC != c:
            curC = c
            if prog is not None:
                prog.set_description(c.ljust(maxC, ' '))
                prog.update(0)
            curSeq = sf[c]
        assert s < len(curSeq)
        assert e <= len(curSeq)
        if prog is not None:
            #prog.set_description(progFmt % (c, m))
            prog.update(n)
        if opts['-V'] and m != prevM:
            prevM = m
            print '# %s : %s' % (c, m)
        for j in xrange(n):
            v = genVar(c, curSeq, s, e, Tgen, Ds)
            fmts = []
            vals = []
            if opts['-T']:
                v.setSequenceFactory(sf)
                wt = curSeq[v.range()[0]:v.range()[1]].upper()
                mut = v.sequence()
                if mut is None:
                    mut = '*'
                fmts += ['%s', '%d', '%d', '%d', '%s', '%s']
                vals += [v.accession(), v.range()[0], v.range()[1], v.size(), wt, mut]
            fmts += ['%s']
            vals += [str(v)]
            print '\t'.join(fmts) % tuple(vals)

if __name__ == '__main__':
    main(sys.argv[1:])
