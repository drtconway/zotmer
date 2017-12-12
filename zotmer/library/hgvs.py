import re

import pytest

from zotmer.library.rope import rope
from zotmer.library.align import revComp

# For testing NM_004119.2
#
_flt3 = \
"ACCTGCAGCGCGAGGCGCGCCGCTCCAGGCGGCATCGCAGGGCTGGGCCGGCGCGGCCTGGGGACCCCGG" + \
"GCTCCGGAGGCCATGCCGGCGTTGGCGCGCGACGGCGGCCAGCTGCCGCTGCTCGTTGTTTTTTCTGCAA" + \
"TGATATTTGGGACTATTACAAATCAAGATCTGCCTGTGATCAAGTGTGTTTTAATCAATCATAAGAACAA" + \
"TGATTCATCAGTGGGGAAGTCATCATCATATCCCATGGTATCAGAATCCCCGGAAGACCTCGGGTGTGCG" + \
"TTGAGACCCCAGAGCTCAGGGACAGTGTACGAAGCTGCCGCTGTGGAAGTGGATGTATCTGCTTCCATCA" + \
"CACTGCAAGTGCTGGTCGACGCCCCAGGGAACATTTCCTGTCTCTGGGTCTTTAAGCACAGCTCCCTGAA" + \
"TTGCCAGCCACATTTTGATTTACAAAACAGAGGAGTTGTTTCCATGGTCATTTTGAAAATGACAGAAACC" + \
"CAAGCTGGAGAATACCTACTTTTTATTCAGAGTGAAGCTACCAATTACACAATATTGTTTACAGTGAGTA" + \
"TAAGAAATACCCTGCTTTACACATTAAGAAGACCTTACTTTAGAAAAATGGAAAACCAGGACGCCCTGGT" + \
"CTGCATATCTGAGAGCGTTCCAGAGCCGATCGTGGAATGGGTGCTTTGCGATTCACAGGGGGAAAGCTGT" + \
"AAAGAAGAAAGTCCAGCTGTTGTTAAAAAGGAGGAAAAAGTGCTTCATGAATTATTTGGGACGGACATAA" + \
"GGTGCTGTGCCAGAAATGAACTGGGCAGGGAATGCACCAGGCTGTTCACAATAGATCTAAATCAAACTCC" + \
"TCAGACCACATTGCCACAATTATTTCTTAAAGTAGGGGAACCCTTATGGATAAGGTGCAAAGCTGTTCAT" + \
"GTGAACCATGGATTCGGGCTCACCTGGGAATTAGAAAACAAAGCACTCGAGGAGGGCAACTACTTTGAGA" + \
"TGAGTACCTATTCAACAAACAGAACTATGATACGGATTCTGTTTGCTTTTGTATCATCAGTGGCAAGAAA" + \
"CGACACCGGATACTACACTTGTTCCTCTTCAAAGCATCCCAGTCAATCAGCTTTGGTTACCATCGTAGAA" + \
"AAGGGATTTATAAATGCTACCAATTCAAGTGAAGATTATGAAATTGACCAATATGAAGAGTTTTGTTTTT" + \
"CTGTCAGGTTTAAAGCCTACCCACAAATCAGATGTACGTGGACCTTCTCTCGAAAATCATTTCCTTGTGA" + \
"GCAAAAGGGTCTTGATAACGGATACAGCATATCCAAGTTTTGCAATCATAAGCACCAGCCAGGAGAATAT" + \
"ATATTCCATGCAGAAAATGATGATGCCCAATTTACCAAAATGTTCACGCTGAATATAAGAAGGAAACCTC" + \
"AAGTGCTCGCAGAAGCATCGGCAAGTCAGGCGTCCTGTTTCTCGGATGGATACCCATTACCATCTTGGAC" + \
"CTGGAAGAAGTGTTCAGACAAGTCTCCCAACTGCACAGAAGAGATCACAGAAGGAGTCTGGAATAGAAAG" + \
"GCTAACAGAAAAGTGTTTGGACAGTGGGTGTCGAGCAGTACTCTAAACATGAGTGAAGCCATAAAAGGGT" + \
"TCCTGGTCAAGTGCTGTGCATACAATTCCCTTGGCACATCTTGTGAGACGATCCTTTTAAACTCTCCAGG" + \
"CCCCTTCCCTTTCATCCAAGACAACATCTCATTCTATGCAACAATTGGTGTTTGTCTCCTCTTCATTGTC" + \
"GTTTTAACCCTGCTAATTTGTCACAAGTACAAAAAGCAATTTAGGTATGAAAGCCAGCTACAGATGGTAC" + \
"AGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATG" + \
"GGAGTTTCCAAGAGAAAATTTAGAGTTTGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAAC" + \
"GCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGGTTGCCGTCAAAATGCTGAAAGAAAAAG" + \
"CAGACAGCTCTGAAAGAGAGGCACTCATGTCAGAACTCAAGATGATGACCCAGCTGGGAAGCCACGAGAA" + \
"TATTGTGAACCTGCTGGGGGCGTGCACACTGTCAGGACCAATTTACTTGATTTTTGAATACTGTTGCTAT" + \
"GGTGATCTTCTCAACTATCTAAGAAGTAAAAGAGAAAAATTTCACAGGACTTGGACAGAGATTTTCAAGG" + \
"AACACAATTTCAGTTTTTACCCCACTTTCCAATCACATCCAAATTCCAGCATGCCTGGTTCAAGAGAAGT" + \
"TCAGATACACCCGGACTCGGATCAAATCTCAGGGCTTCATGGGAATTCATTTCACTCTGAAGATGAAATT" + \
"GAATATGAAAACCAAAAAAGGCTGGAAGAAGAGGAGGACTTGAATGTGCTTACATTTGAAGATCTTCTTT" + \
"GCTTTGCATATCAAGTTGCCAAAGGAATGGAATTTCTGGAATTTAAGTCGTGTGTTCACAGAGACCTGGC" + \
"CGCCAGGAACGTGCTTGTCACCCACGGGAAAGTGGTGAAGATATGTGACTTTGGATTGGCTCGAGATATC" + \
"ATGAGTGATTCCAACTATGTTGTCAGGGGCAATGCCCGTCTGCCTGTAAAATGGATGGCCCCCGAAAGCC" + \
"TGTTTGAAGGCATCTACACCATTAAGAGTGATGTCTGGTCATATGGAATATTACTGTGGGAAATCTTCTC" + \
"ACTTGGTGTGAATCCTTACCCTGGCATTCCGGTTGATGCTAACTTCTACAAACTGATTCAAAATGGATTT" + \
"AAAATGGATCAGCCATTTTATGCTACAGAAGAAATATACATTATAATGCAATCCTGCTGGGCTTTTGACT" + \
"CAAGGAAACGGCCATCCTTCCCTAATTTGACTTCGTTTTTAGGATGTCAGCTGGCAGATGCAGAAGAAGC" + \
"GATGTATCAGAATGTGGATGGCCGTGTTTCGGAATGTCCTCACACCTACCAAAACAGGCGACCTTTCAGC" + \
"AGAGAGATGGATTTGGGGCTACTCTCTCCGCAGGCTCAGGTCGAAGATTCGTAGAGGAACAATTTAGTTT" + \
"TAAGGACTTCATCCCTCCACCTATCCCTAACAGGCTGTAGATTACCAAAACAAGATTAATTTCATCACTA" + \
"AAAGAAAATCTATTATCAACTGCTGCTTCACCAGACTTTTCTCTAGAAGCTGTCTGCGTTTACTCTTGTT" + \
"TTCAAAGGGACTTTTGTAAAATCAAATCATCCTGTCACAAGGCAGGAGGAGCTGATAATGAACTTTATTG" + \
"GAGCATTGATCTGCATCCAAGGCCTTCTCAGGCTGGCTTGAGTGAATTGTGTACCTGAAGTACAGTATAT" + \
"TCTTGTAAATACATAAAACAAAAGCATTTTGCTAAGGAGAAGCTAATATGATTTTTTAAGTCTATGTTTT" + \
"AAAATAATATGTAAATTTTTCAGCTATTTAGTGATATATTTTATGGGTGGGAATAAAATTTCTACTACAG" + \
"AATTGCCCATTATTGAATTATTTACATGGTATAATTAGGGCAAGTCTTAACTGGAGTTCACGAACCCCCT" + \
"GAAATTGTGCACCCATAGCCACCTACACATTCCTTCCAGAGCACGTGTGCTTTTACCCCAAGATACAAGG" + \
"AATGTGTAGGCAGCTATGGTTGTCACAGCCTAAGATTTCTGCAACAACAGGGGTTGTATTGGGGGAAGTT" + \
"TATAATGAATAGGTGTTCTACCATAAAGAGTAATACATCACCTAGACACTTTGGCGGCCTTCCCAGACTC" + \
"AGGGCCAGTCAGAAGTAACATGGAGGATTAGTATTTTCAATAAAGTTACTCTTGTCCCCACAAAAAAA"

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
    """
    An HGVS object describes a variation.

    Unlike the textual notation used in the specification,
    internally positions are 0 based not 1 based.

    Each variant may be described in terms of a section of
    sequence that is removed from the reference - an interval
    returned by the method range(); and a section of sequence
    that replaces it - the actual nucleotides are returned by
    the method sequence() and the method size() returns the
    length of the inserted sequence.

    This representation is used to map the high level notion
    of the variant to the rope representation used underneath.
    """
    def __init__(self, acc, seqFac):
        self.acc = acc
        self.seqFac = seqFac

    def accession(self):
        return self.acc

    def setSequenceFactory(self, seqFac):
        self.seqFac = seqFac

    def loadAccession(self):
        assert self.seqFac is not None
        return self.seqFac[self.accession()]

    def cryptic(self):
        return None

    def anonymous(self):
        return False

    def liftover(self, p):
        (b, e) = self.range()
        d = self.size() - (e - b)

        if p < b:
            return p
        q = p + d
        if q < b:
            return b
        else:
            return q

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

    def context(self, w):
        seq = self.loadAccession()
        s = rope.atom(seq)

        rng = self.range()

        wt = (rng[0]-w+1, rng[1]+w-1)
        wt = (max(0, wt[0]), min(wt[1], len(s)))

        wtSeq = s[wt[0]:wt[1]]

        n = self.size()
        if n is not None:
            r = self.apply(s)
            d = rng[0] - rng[1] + n
            mut = (rng[0], rng[1] + d)

            mut = (mut[0]-w+1, mut[1]+w-1)
            mut = (max(0, mut[0]), min(mut[1], len(r)))

            mutSeq = r[mut[0]:mut[1]]
        else:
            # must have been cryptic
            mutSeq = None

        return (wtSeq, mutSeq)

    def __cmp__(self, other):
        if self.accession() != other.accession():
            return cmp(self.accession(), other.accession())
        xr = self.range()
        yr = other.range()
        if xr == yr:
            return 0
        if xr[1] < yr[0]:
            return -1
        if yr[1] < xr[0]:
            return 1
        return NotImplemented

    def overlaps(self, other):
        if self.accession() != other.accession():
            return False

        r = cmp(self, other)
        if r == NotImplemented or r == 0:
            return True
        return False

class Substitution(HGVS):
    def __init__(self, acc, pos, ref, var, sf = None):
        super(Substitution, self).__init__(acc, sf)
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
        return '%s:g.%d%s>%s' % (self.acc, self.pos+1, self.ref, self.var)

def test_substitution0():
    sf = {'NM_004119.2': _flt3}
    assert _flt3[49] == 'G'
    v = Substitution('NM_004119.2', 49, 'G', 'T', sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGCTGGGCCGGCGCGGCCT'
    assert ctx[1] == 'GGCTGGGCCTGCGCGGCCT'

def test_substitution1():
    sf = {'NM_004119.2': _flt3}
    h = 'NM_004119.2:g.50G>T'
    v = makeHGVS(h, sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGCTGGGCCGGCGCGGCCT'
    assert ctx[1] == 'GGCTGGGCCTGCGCGGCCT'
    assert str(v) == h

class Deletion(HGVS):
    def __init__(self, acc, fst, lst, sf = None):
        super(Deletion, self).__init__(acc, sf)
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
            return '%s:g.%ddel' % (self.acc, self.fst+1)
        return '%s:g.%d_%ddel' % (self.acc, self.fst+1, self.lst+1)

def test_deletion0():
    sf = {'NM_004119.2': _flt3}
    assert _flt3[49] == 'G'
    v = Deletion('NM_004119.2', 47, 51, sf)
    ctx = v.context(10)
    assert ctx[0] == 'AGGGCTGGGCCGGCGCGGCCTGG'
    assert ctx[1] == 'AGGGCTGGGGCGGCCTGG'

def test_deletion1():
    sf = {'NM_004119.2': _flt3}
    h = 'NM_004119.2:g.48_52del'
    v = makeHGVS(h, sf)
    ctx = v.context(10)
    assert ctx[0] == 'AGGGCTGGGCCGGCGCGGCCTGG'
    assert ctx[1] == 'AGGGCTGGGGCGGCCTGG'
    assert str(v) == h

class Insertion(HGVS):
    def __init__(self, acc, bef, seq, sf = None):
        super(Insertion, self).__init__(acc, sf)
        self.bef = bef
        self.seq = seq

    def size(self):
        return len(self.seq)

    def range(self):
        return (self.bef, self.bef)

    def sequence(self):
        return self.seq

    def __str__(self):
        return '%s:g.%d_%dins%s' % (self.acc, self.bef, self.bef+1, self.seq)

def test_insertion0():
    sf = {'NM_004119.2': _flt3}
    v = Insertion('NM_004119.2', 49, 'TTATT', sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGCTGGGCCGGCGCGGCC'
    assert ctx[1] == 'GGCTGGGCCTTATTGGCGCGGCC'

def test_insertion1():
    sf = {'NM_004119.2': _flt3}
    h = 'NM_004119.2:g.49_50insTTATT'
    v = makeHGVS(h, sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGCTGGGCCGGCGCGGCC'
    assert ctx[1] == 'GGCTGGGCCTTATTGGCGCGGCC'
    assert str(v) == h

class Anonymous(HGVS):
    def __init__(self, acc, bef, iln, sf = None):
        super(Anonymous, self).__init__(acc, sf)
        self.bef = bef
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
        return (self.bef, self.bef)

    def sequence(self):
        if self.seq is not None:
            return self.seq
        return self.iln * 'N'

    def __str__(self):
        return '%s:g.%d_%dins%d' % (self.acc, self.bef, self.bef+1, self.iln)

def test_anonymous0():
    sf = {'NM_004119.2': _flt3}
    v = Anonymous('NM_004119.2', 49, 5, sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGCTGGGCCGGCGCGGCC'
    assert ctx[1] == 'GGCTGGGCCNNNNNGGCGCGGCC'

def test_anonymous1():
    sf = {'NM_004119.2': _flt3}
    v = Anonymous('NM_004119.2', 49, 5, sf)
    v.setSequence('TTATT')
    ctx = v.context(10)
    assert ctx[0] == 'GGCTGGGCCGGCGCGGCC'
    assert ctx[1] == 'GGCTGGGCCTTATTGGCGCGGCC'

def test_anonymous2():
    sf = {'NM_004119.2': _flt3}
    h = 'NM_004119.2:g.49_50ins5'
    v = makeHGVS(h, sf)
    v.setSequence('TTATT')
    ctx = v.context(10)
    assert ctx[0] == 'GGCTGGGCCGGCGCGGCC'
    assert ctx[1] == 'GGCTGGGCCTTATTGGCGCGGCC'
    assert str(v) == h

class DeletionInsertion(HGVS):
    def __init__(self, acc, fst, lst, seq, sf = None):
        super(DeletionInsertion, self).__init__(acc, sf)
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
            return '%s:g.%ddelins%s' % (self.acc, self.fst+1, self.seq)
        return '%s:g.%d_%ddelins%s' % (self.acc, self.fst+1, self.lst+1, self.seq)

def test_delins0():
    sf = {'NM_004119.2': _flt3}
    v = DeletionInsertion('NM_004119.2', 48, 49, 'TTATT', sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGGCTGGGCCGGCGCGGCCT'
    assert ctx[1] == 'GGGCTGGGCTTATTGCGCGGCCT'

def test_delins1():
    sf = {'NM_004119.2': _flt3}
    h = 'NM_004119.2:g.49_50delinsTTATT'
    v = makeHGVS(h, sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGGCTGGGCCGGCGCGGCCT'
    assert ctx[1] == 'GGGCTGGGCTTATTGCGCGGCCT'
    assert str(v) == h

class AnonymousDelIns(HGVS):
    def __init__(self, acc, fst, lst, iln, sf = None):
        super(AnonymousDelIns, self).__init__(acc, sf)
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
            return '%s:g.%ddelins%d' % (self.acc, self.fst+1, self.iln)
        return '%s:g.%d_%ddelins%d' % (self.acc, self.fst+1, self.lst+1, self.iln)

def test_anondelins0():
    sf = {'NM_004119.2': _flt3}
    v = AnonymousDelIns('NM_004119.2', 48, 49, 5, sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGGCTGGGCCGGCGCGGCCT'
    assert ctx[1] == 'GGGCTGGGCNNNNNGCGCGGCCT'

def test_anondelins1():
    sf = {'NM_004119.2': _flt3}
    v = AnonymousDelIns('NM_004119.2', 48, 49, 5, sf)
    v.setSequence('TTATT')
    ctx = v.context(10)
    assert ctx[0] == 'GGGCTGGGCCGGCGCGGCCT'
    assert ctx[1] == 'GGGCTGGGCTTATTGCGCGGCCT'

def test_anondelins2():
    sf = {'NM_004119.2': _flt3}
    h = 'NM_004119.2:g.49_50delins5'
    v = makeHGVS(h, sf)
    v.setSequence('TTATT')
    ctx = v.context(10)
    assert ctx[0] == 'GGGCTGGGCCGGCGCGGCCT'
    assert ctx[1] == 'GGGCTGGGCTTATTGCGCGGCCT'
    assert str(v) == h

class Repeat(HGVS):
    def __init__(self, acc, fst, lst, cnt, sf = None):
        super(Repeat, self).__init__(acc, sf)
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
            self.seq = big[self.fst:self.lst+1].upper()
        return self.cnt * self.seq

    def __str__(self):
        if self.fst == self.lst:
            return '%s:g.%d[%d]' % (self.acc, self.fst+1, self.cnt)
        return '%s:g.%d_%d[%d]' % (self.acc, self.fst+1, self.lst+1, self.cnt)

def test_repeat0():
    sf = {'NM_004119.2': _flt3}
    v = Repeat('NM_004119.2', 48, 49, 5, sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGGCTGGGCCGGCGCGGCCT'
    assert ctx[1] == 'GGGCTGGGCCGCGCGCGCGGCGCGGCCT'

def test_repeat1():
    sf = {'NM_004119.2': _flt3}
    h = 'NM_004119.2:g.49_50[5]'
    v = makeHGVS(h, sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGGCTGGGCCGGCGCGGCCT'
    assert ctx[1] == 'GGGCTGGGCCGCGCGCGCGGCGCGGCCT'
    assert str(v) == h

class Duplication(HGVS):
    def __init__(self, acc, fst, lst, sf = None):
        super(Duplication, self).__init__(acc, sf)
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
            self.seq = big[self.fst:self.lst+1].upper()
        return 2*self.seq

    def __str__(self):
        if self.fst == self.lst:
            return '%s:g.%ddup' % (self.acc, self.fst+1)
        return '%s:g.%d_%ddup' % (self.acc, self.fst+1, self.lst+1)

def test_duplication0():
    sf = {'NM_004119.2': _flt3}
    v = Duplication('NM_004119.2', 48, 49, sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGGCTGGGCCGGCGCGGCCT'
    assert ctx[1] == 'GGGCTGGGCCGCGGCGCGGCCT'

def test_duplication1():
    sf = {'NM_004119.2': _flt3}
    h = 'NM_004119.2:g.49_50dup'
    v = makeHGVS(h, sf)
    ctx = v.context(10)
    assert ctx[0] == 'GGGCTGGGCCGGCGCGGCCT'
    assert ctx[1] == 'GGGCTGGGCCGCGGCGCGGCCT'
    assert str(v) == h

class Inversion(HGVS):
    def __init__(self, acc, fst, lst, sf = None):
        super(Inversion, self).__init__(acc, sf)
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
            self.seq = revComp(big[self.fst:self.lst+1].upper())
        return self.seq

    def __str__(self):
        return '%s:g.%d_%dinv' % (self.acc, self.fst+1, self.lst+1)

def test_inversion0():
    sf = {'NM_004119.2': _flt3}
    v = Inversion('NM_004119.2', 50, 61, sf)
    ctx = v.context(10)
    assert ctx[0] == 'GCTGGGCCGGCGCGGCCTGGGGACCCCGGG'
    assert ctx[1] == 'GCTGGGCCGCCCAGGCCGCGCGACCCCGGG'

def test_inversion1():
    sf = {'NM_004119.2': _flt3}
    h = 'NM_004119.2:g.51_62inv'
    v = makeHGVS(h, sf)
    ctx = v.context(10)
    assert ctx[0] == 'GCTGGGCCGGCGCGGCCTGGGGACCCCGGG'
    assert ctx[1] == 'GCTGGGCCGCCCAGGCCGCGCGACCCCGGG'
    assert str(v) == h

class AluInsertion(HGVS):
    def __init__(self, acc, bef, sf = None):
        super(AluInsertion, self).__init__(acc, sf)
        self.bef = bef

    def cryptic(self):
        return 'Alu'

    def size(self):
        return None

    def range(self):
        return (self.bef, self.bef)

    def sequence(self):
        return None

    def __str__(self):
        return '%s:g.%d_%dinsALU' % (self.acc, self.bef, self.bef+1)

gvarPfx = re.compile('([^:]+):g\.(.*)$')
gvarPos = re.compile('([0-9]+)(_([0-9]+))?(.*)$')
gvarSub = re.compile('([ACGT])>([ACGT])$')
gvarDel = re.compile('del([ACGT]+)?$')
gvarDup = re.compile('dup([ACGT]+)?$')
gvarIns = re.compile('ins([ACGT]+)$')
gvarAnI = re.compile('ins([0-9]+)$')
gvarALU = re.compile('insALU$')
gvarInv = re.compile('inv$')
gvarDIn = re.compile('delins([ACGT]+)$')
gvarAnD = re.compile('delins([0-9]+)$')
gvarRpt = re.compile('([ACGT]+)?\[([0-9]+)\]$')

def makeHGVS(txt, sf = None):
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
        return Substitution(acc, pos0-1, ref, alt, sf)

    m = gvarDel.match(s)
    if m:
        if pos1 is None:
            pos1 = pos0
        else:
            pos1 = int(pos1)
        return Deletion(acc, pos0-1, pos1-1, sf)

    m = gvarDup.match(s)
    if m:
        if pos1 is None:
            pos1 = pos0
        else:
            pos1 = int(pos1)
        return Duplication(acc, pos0-1, pos1-1, sf)

    m = gvarIns.match(s)
    if m:
        if pos1 is None:
            return None
        pos1 = int(pos1)
        ss = m.groups()
        return Insertion(acc, pos1-1, ss[0], sf)

    m = gvarAnI.match(s)
    if m:
        if pos1 is None:
            return None
        pos1 = int(pos1)
        ss = m.groups()
        l = int(ss[0])
        return Anonymous(acc, pos1-1, l, sf)

    m = gvarALU.match(s)
    if m:
        if pos1 is None:
            return None
        pos1 = int(pos1)
        return AluInsertion(acc, pos1-1, sf)

    m = gvarInv.match(s)
    if m:
        if pos1 is None:
            return None
        pos1 = int(pos1)
        return Inversion(acc, pos0-1, pos1-1, sf)

    m = gvarDIn.match(s)
    if m:
        if pos1 is None:
            pos1 = pos0
        else:
            pos1 = int(pos1)
        ss = m.groups()
        return DeletionInsertion(acc, pos0-1, pos1-1, ss[0], sf)

    m = gvarAnD.match(s)
    if m:
        if pos1 is None:
            pos1 = pos0
        else:
            pos1 = int(pos1)
        ss = m.groups()
        l = int(ss[0])
        return AnonymousDelIns(acc, pos0-1, pos1-1, l, sf)

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
        return Repeat(acc, pos0-1, pos1-1, cnt, sf)

    return None

