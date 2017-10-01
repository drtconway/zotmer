import re

from zotmer.library.rope import rope

hg19Data = [("chr1",   "CM000663.1",  "NC_000001.10"),
            ("chr2",   "CM000664.1",  "NC_000002.11"),
            ("chr3",   "CM000665.1",  "NC_000003.11"),
            ("chr4",   "CM000666.1",  "NC_000004.11"),
            ("chr5",   "CM000667.1",  "NC_000005.9"),
            ("chr6",   "CM000668.1",  "NC_000006.11"),
            ("chr7",   "CM000669.1",  "NC_000007.13"),
            ("chr8",   "CM000670.1",  "NC_000008.10"),
            ("chr9",   "CM000671.1",  "NC_000009.11"),
            ("chr10",  "CM000672.1",  "NC_000010.10"),
            ("chr11",  "CM000673.1",  "NC_000011.9"),
            ("chr12",  "CM000674.1",  "NC_000012.11"),
            ("chr13",  "CM000675.1",  "NC_000013.10"),
            ("chr14",  "CM000676.1",  "NC_000014.8"),
            ("chr15",  "CM000677.1",  "NC_000015.9"),
            ("chr16",  "CM000678.1",  "NC_000016.9"),
            ("chr17",  "CM000679.1",  "NC_000017.10"),
            ("chr18",  "CM000680.1",  "NC_000018.9"),
            ("chr19",  "CM000681.1",  "NC_000019.9"),
            ("chr20",  "CM000682.1",  "NC_000020.10"),
            ("chr21",  "CM000683.1",  "NC_000021.8"),
            ("chr22",  "CM000684.1",  "NC_000022.10"),
            ("chrX",   "CM000685.1",  "NC_000023.10"),
            ("chrY",   "CM000686.1",  "NC_000024.9")]

refSeq2Hg19 = dict([(r, h) for (h, g, r) in hg19Data])
hg19ToRefSeq= dict([(h, r) for (h, g, r) in hg19Data])

class HGVS(object):
    def __init__(self, acc):
        self.acc = acc
        self.seqFac = None

    def accession(self):
        return self.acc

    def setSequenceFactory(self, seqFac):
        self.seqFac = seqFac

    def loadAccession(self):
        assert self.seqFac is not None
        return self.seqFac[self.accession()]

    def apply(self, s):
        (p,q) = self.range()
        v = self.sequence()
        a = rope.substr(s, 0, p)
        b = rope.substr(s, q, len(s))
        if len(v) > 0:
            c = rope.atom(v)
            return rope.join([a, c, b])
        else:
            return rope.concat(a, b)

class Substitution(HGVS):
    def __init__(self, acc, pos, ref, var):
        super(Substitution, self).__init__(acc)
        self.pos = pos
        self.ref = ref
        self.var = var

    def range(self):
        return (self.pos, self.pos+1)

    def sequence(self):
        return self.var

    def __str__(self):
        return '%s:g.%d%s>%s' % (self.acc, self.pos, self.ref, self.var)

class Deletion(HGVS):
    def __init__(self, acc, fst, lst):
        super(Deletion, self).__init__(acc)
        self.fst = fst
        self.lst = lst

    def range(self):
        return (self.fst, self.lst+1)

    def sequence(self):
        return ''

    def __str__(self):
        return '%s:g.%d_%ddel' % (self.acc, self.fst, self.lst)

class Insertion(HGVS):
    def __init__(self, acc, aft, bef, seq):
        super(Insertion, self).__init__(acc)
        self.aft = aft
        assert self.aft + 1 == self.bef
        self.seq = seq

    def range(self):
        return (self.aft, self.aft)

    def sequence(self):
        return self.seq

    def __str__(self):
        return '%s:g.%d_%dins%s' % (self.acc, self.aft, self.aft+1, self.seq)

class DeletionInsertion(HGVS):
    def __init__(self, acc, fst, lst, seq):
        super(DeletionInsertion, self).__init__(acc)
        self.fst = fst
        self.lst = lst
        self.seq = seq

    def range(self):
        return (self.fst, self.lst+1)

    def sequence(self):
        return self.seq

    def __str__(self):
        return '%s:g.%d_%ddelins%s' % (self.acc, self.fst, self.lst, self.seq)

class Repeat(HGVS):
    def __init__(self, acc, fst, lst, cnt):
        super(Repeat, self).__init__(acc)
        self.fst = fst
        self.lst = lst
        self.cnt = cnt
        self.seq = None

    def range(self):
        return (self.fst, self.lst+1)

    def sequence(self):
        if not self.seq:
            big = self.loadAccession()
            self.seq = big[self.fst-1:self.lst]
        return self.cnt * self.seq

    def __str__(self):
        return '%s:g.%d_%d[%d]' % (self.acc, self.fst, self.lst, self.cnt)

class Duplication(HGVS):
    def __init__(self, acc, fst, lst):
        super(Duplication, self).__init__(acc)
        self.fst = fst
        self.lst = lst
        self.seq = None

    def range(self):
        return (self.fst, self.lst+1)

    def sequence(self):
        if not self.seq:
            big = self.loadAccession()
            self.seq = big[self.fst-1:self.lst]
        return 2*self.seq

    def __str__(self):
        return '%s:g.%d_%ddup' % (self.acc, self.fst, self.lst)

gvarSub = re.compile("(\w+\.\d+):g\.(\d+)([ACGT])>([ACGT])$")
gvarRpt = re.compile("(\w+\.\d+):g\.(\d+)_(\d+)([ACGT]*)\[(\d+)\]$")
gvarIns = re.compile("(\w+\.\d+):g\.(\d+)_(\d+)ins([ACGT]+)$")
gvarDel0 = re.compile("(\w+\.\d+):g\.(\d+)del([ACGT])?$")
gvarDel1 = re.compile("(\w+\.\d+):g\.(\d+)_(\d+)del([ACGT]+)?$")
gvarDel2 = re.compile("(\w+\.\d+):g\.(\d+)_(\d+)del(\d+)?$")
gvarDup = re.compile("(\w+\.\d+):g\.(\d+)_(\d+)dup$")
gvarDin = re.compile("(\w+\.\d+):g\.(\d+)_(\d+)delins([ACGT]+)$")

def parseHGVS(s):
    m = gvarSub.match(s)
    if m:
        r = {}
        r['type'] = 'substitution'
        r['accession'] = m.group(1)
        r['position'] = int(m.group(2))
        r['reference'] = m.group(3)
        r['variant'] = m.group(4)
        return r

    m = gvarRpt.match(s)
    if m:
        r = {}
        r['type'] = 'repeat'
        r['accession'] = m.group(1)
        r['start-position'] = int(m.group(2))
        r['stop-position'] = int(m.group(3))
        r['sequence'] = m.group(4)
        r['count'] = int(m.group(5))
        return r

    m = gvarIns.match(s)
    if m:
        r = {}
        r['type'] = 'insertion'
        r['accession'] = m.group(1)
        r['after-position'] = int(m.group(2))
        r['before-position'] = int(m.group(3))
        r['sequence'] = m.group(4)
        return r

    m = gvarDel0.match(s)
    if m:
        r = {}
        r['type'] = 'deletion'
        r['accession'] = m.group(1)
        r['first-position'] = int(m.group(2))
        r['last-position'] = int(m.group(2))
        return r

    m = gvarDel1.match(s)
    if m:
        r = {}
        r['type'] = 'deletion'
        r['accession'] = m.group(1)
        r['first-position'] = int(m.group(2))
        r['last-position'] = int(m.group(3))
        return r

    m = gvarDel2.match(s)
    if m:
        r = {}
        r['type'] = 'deletion'
        r['accession'] = m.group(1)
        r['first-position'] = int(m.group(2))
        r['last-position'] = int(m.group(3))
        return r

    m = gvarDup.match(s)
    if m:
        r = {}
        r['type'] = 'duplication'
        r['accession'] = m.group(1)
        r['first-position'] = int(m.group(2))
        r['last-position'] = int(m.group(3))
        return r

    m = gvarDin.match(s)
    if m:
        r = {}
        r['type'] = 'deletion-insertion'
        r['accession'] = m.group(1)
        r['first-position'] = int(m.group(2))
        r['last-position'] = int(m.group(3))
        r['insertion'] = m.group(4)
        return r

    return None

def makeHGVS(s):
    v = parseHGVS(s)
    if v['type'] == 'substitution':
        return Substitution(v['accession'], v['position'], v['reference'], v['variant'])
    if v['type'] == 'repeat':
        return Repeat(v['accession'], v['start-position'], v['stop-position'], v['count'])
    if v['type'] == 'insertion':
        return Insertion(v['accession'], v['after-position'], v['before-position'], v['sequence'])
    if v['type'] == 'deletion':
        return Deletion(v['accession'], v['first-position'], v['last-position'])
    if v['type'] == 'duplication':
        return Duplication(v['accession'], v['first-position'], v['last-position'])
    if v['type'] == 'deletion-insertion':
        return DeletionInsertion(v['accession'], v['first-position'], v['last-position'], v['sequence'])

def applyVariant(s, v):
    if v['type'] == 'substitution':
        p = v['position'] - 1
        l = rope.substr(s, 0, p)
        m = rope.atom(v['variant'])
        r = rope.substr(s, p + 1, len(s))
        return rope.join([l, m, r])
    elif v['type'] == 'insertion':
        p = v['after-position']
        q = v['before-position'] - 1
        assert p == q
        l = rope.substr(s, 0, p)
        m = rope.atom(v['sequence'])
        r = rope.substr(s, p, len(s))
        return rope.join([l, m, r])
    elif v['type'] == 'deletion':
        p = v['first-position'] - 1
        q = v['last-position']
        l = rope.substr(s, 0, p)
        r = rope.substr(s, q, len(s))
        return rope.concat(l, r)
    return None

