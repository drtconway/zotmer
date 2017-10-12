import re

from zotmer.library.rope import rope
from zotmer.library.align import revComp

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

    def anonymous(self):
        return False

    def apply(self, s):
        (p,q) = self.range()
        p -= 1
        q -= 1
        v = self.sequence()
        a = rope.substr(s, 0, p)
        b = rope.substr(s, q, len(s))
        if len(v) > 0:
            c = rope.atom(v)
            return rope.join([a, c, b])
        else:
            return rope.concat(a, b)

    def context(self, w):
        wt = self.range()
        wt = (wt[0] - 1, wt[1] - 1) # convert to zero offset
        n = self.size()
        d = wt[0] - wt[1] + n
        mut = (wt[0], wt[1] + d)

        seq = self.loadAccession()
        s = rope.atom(seq)

        r = self.apply(s)

        wt = (wt[0]-w+1, wt[1]+w-1)
        wt = (max(0, wt[0]), min(wt[1], len(s)))

        mut = (mut[0]-w+1, mut[1]+w-1)
        mut = (max(0, mut[0]), min(mut[1], len(r)))

        wtSeq = s[wt[0]:wt[1]]
        mutSeq = r[mut[0]:mut[1]]

        return (wtSeq, mutSeq)

class Substitution(HGVS):
    def __init__(self, acc, pos, ref, var):
        super(Substitution, self).__init__(acc)
        self.pos = pos
        self.ref = ref
        self.var = var

    def size(self):
        return 1

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

    def size(self):
        return 0

    def range(self):
        return (self.fst, self.lst+1)

    def sequence(self):
        return ''

    def __str__(self):
        if self.fst == self.lst:
            return '%s:g.%ddel' % (self.acc, self.fst)
        return '%s:g.%d_%ddel' % (self.acc, self.fst, self.lst)

class Insertion(HGVS):
    def __init__(self, acc, aft, bef, seq):
        super(Insertion, self).__init__(acc)
        self.aft = aft
        assert self.aft + 1 == bef
        self.seq = seq

    def size(self):
        return len(self.seq)

    def range(self):
        return (self.aft, self.aft)

    def sequence(self):
        return self.seq

    def __str__(self):
        return '%s:g.%d_%dins%s' % (self.acc, self.aft, self.aft+1, self.seq)

class Anonymous(HGVS):
    def __init__(self, acc, aft, bef, iln):
        super(Anonymous, self).__init__(acc)
        self.aft = aft
        assert self.aft + 1 == bef
        self.iln = iln
        self.seq = None

    def anonymous(self):
        return True

    def setSequence(self, seq):
        assert len(seq) == self.iln
        self.seq = seq

    def size(self):
        return self.iln

    def range(self):
        return (self.aft, self.aft)

    def sequence(self):
        if self.seq is not None:
            return self.seq
        return self.iln * 'N'

    def __str__(self):
        return '%s:g.%d_%dins%d' % (self.acc, self.aft, self.aft+1, self.iln)

class DeletionInsertion(HGVS):
    def __init__(self, acc, fst, lst, seq):
        super(DeletionInsertion, self).__init__(acc)
        self.fst = fst
        self.lst = lst
        self.seq = seq

    def size(self):
        return len(self.seq)

    def range(self):
        return (self.fst, self.lst+1)

    def sequence(self):
        return self.seq

    def __str__(self):
        if self.fst == self.lst:
            return '%s:g.%ddelins%s' % (self.acc, self.fst, self.seq)
        return '%s:g.%d_%ddelins%s' % (self.acc, self.fst, self.lst, self.seq)

class AnonymousDelIns(HGVS):
    def __init__(self, acc, fst, lst, iln):
        super(AnonymousDelIns, self).__init__(acc)
        self.fst = fst
        self.lst = lst
        self.iln = iln
        self.seq = None

    def anonymous(self):
        return True

    def setSequence(self, seq):
        assert len(seq) == self.iln
        self.seq = seq

    def size(self):
        return self.iln

    def range(self):
        return (self.fst, self.lst+1)

    def sequence(self):
        if self.seq is not None:
            return self.seq
        return self.iln * 'N'

    def __str__(self):
        if self.fst == self.lst:
            return '%s:g.%ddelins%d' % (self.acc, self.fst, self.iln)
        return '%s:g.%d_%ddelins%d' % (self.acc, self.fst, self.lst, self.iln)

class Repeat(HGVS):
    def __init__(self, acc, fst, lst, cnt):
        super(Repeat, self).__init__(acc)
        self.fst = fst
        self.lst = lst
        self.cnt = cnt
        self.seq = None

    def size(self):
        return self.cnt * (self.lst - self.fst + 1)

    def range(self):
        return (self.fst, self.lst+1)

    def sequence(self):
        if not self.seq:
            big = self.loadAccession()
            self.seq = big[self.fst-1:self.lst]
        return self.cnt * self.seq

    def __str__(self):
        if self.fst == self.lst:
            return '%s:g.%d[%d]' % (self.acc, self.fst, self.cnt)
        return '%s:g.%d_%d[%d]' % (self.acc, self.fst, self.lst, self.cnt)

class Duplication(HGVS):
    def __init__(self, acc, fst, lst):
        super(Duplication, self).__init__(acc)
        self.fst = fst
        self.lst = lst
        self.seq = None

    def size(self):
        return 2 * (self.lst - self.fst + 1)

    def range(self):
        return (self.fst, self.lst+1)

    def sequence(self):
        if not self.seq:
            big = self.loadAccession()
            self.seq = big[self.fst-1:self.lst]
        return 2*self.seq

    def __str__(self):
        if self.fst == self.lst:
            return '%s:g.%ddup' % (self.acc, self.fst)
        return '%s:g.%d_%ddup' % (self.acc, self.fst, self.lst)

class Inversion(HGVS):
    def __init__(self, acc, fst, lst):
        super(Inversion, self).__init__(acc)
        self.fst = fst
        self.lst = lst
        self.seq = None

    def size(self):
        return self.lst - self.fst + 1

    def range(self):
        return (self.fst, self.lst+1)

    def sequence(self):
        if not self.seq:
            big = self.loadAccession()
            self.seq = revComp(big[self.fst-1:self.lst].upper())
        return self.seq

    def __str__(self):
        return '%s:g.%d_%dinv' % (self.acc, self.fst, self.lst)

gvarPfx = re.compile('([^:]+):g\.(.*)$')
gvarPos = re.compile('([0-9]+)(_([0-9]+))?(.*)$')
gvarSub = re.compile('([ACGT])>([ACGT])$')
gvarDel = re.compile('del([ACGT]+)?$')
gvarDup = re.compile('dup([ACGT]+)?$')
gvarIns = re.compile('ins([ACGT]+)$')
gvarAnI = re.compile('ins([0-9]+)$')
gvarInv = re.compile('inv$')
gvarDIn = re.compile('delins([ACGT]+)$')
gvarAnD = re.compile('delins([0-9]+)$')
gvarRpt = re.compile('([ACGT]+)?\[([0-9]+)\]$')

def makeHGVS(txt):
    s = txt
    mpfx = gvarPfx.match(s)
    if mpfx is None:
        return None
    ss = mpfx.groups()
    acc = ss[0]
    s = ss[1]
    mpos = gvarPos.match(s)
    if mpos is None:
        return None
    ss = mpos.groups()
    pos0 = int(ss[0])
    pos1 = ss[2]
    s = ss[3]

    m = gvarSub.match(s)
    if m:
        if pos1 is not None:
            return None
        ss = m.groups()
        ref = ss[0]
        alt = ss[1]
        return Substitution(acc, pos0, ref, alt)

    m = gvarDel.match(s)
    if m:
        if pos1 is None:
            pos1 = pos0
        else:
            pos1 = int(pos1)
        return Deletion(acc, pos0, pos1)

    m = gvarDup.match(s)
    if m:
        if pos1 is None:
            pos1 = pos0
        else:
            pos1 = int(pos1)
        return Duplication(acc, pos0, pos1)

    m = gvarIns.match(s)
    if m:
        if pos1 is None:
            return None
        pos1 = int(pos1)
        ss = m.groups()
        return Insertion(acc, pos0, pos1, ss[0])

    m = gvarAnI.match(s)
    if m:
        if pos1 is None:
            return None
        pos1 = int(pos1)
        ss = m.groups()
        l = int(ss[0])
        return Anonymous(acc, pos0, pos1, l)

    m = gvarInv.match(s)
    if m:
        if pos1 is None:
            return None
        pos1 = int(pos1)
        return Inversion(acc, pos0, pos1)

    m = gvarDIn.match(s)
    if m:
        if pos1 is None:
            pos1 = pos0
        else:
            pos1 = int(pos1)
        ss = m.groups()
        return DeletionInsertion(acc, pos0, pos1, ss[0])

    m = gvarAnD.match(s)
    if m:
        if pos1 is None:
            pos1 = pos0
        else:
            pos1 = int(pos1)
        ss = m.groups()
        l = int(ss[0])
        return AnonymousDelIns(acc, pos0, pos1, l)

    m = gvarRpt.match(s)
    if m:
        ss = m.groups()
        if ss[0] is not None and pos1 is None:
            pos1 = pos0 + len(ss[0]) - 1
        elif pos1 is None:
            pos1 = pos0
        else:
            pos1 = int(pos1)
        cnt = int(ss[1])
        return Repeat(acc, pos0, pos1, cnt)

    return None

