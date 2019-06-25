"""
zot hgvs-find - search for known variants in read data.

Usage:
    zot fizzle -P [options] <gene-list> <refgene> <bedfile>
    zot fizzle -X [options] <index> <must-have>
    zot fizzle -T [options] <index>
    zot fizzle [options] <index> <input>...

Options:
    -P              prepare a bed file
    -X              index exons
    -T              test the index for crosstalk
    -g PATH         directory of FASTA reference sequences
    -k K            value of k to use [default: 25]
    -M NUM          minimum number of reads to report gene group [default: 1]
    -m NUM          minimum number of reads to report exon group [default: 1]
    -R NUM          minimum number of k-mers/read to report gene group [default: 1]
    -r NUM          minimum number of k-mers/read to report exon group [default: 1]
    -t              when preparing the bedfile, read genes and transcripts, not just genes
    -Z MIN:MAX      minimum:maximum number of genes to report [default: 0:999999]
    -z MIN:MAX      minimum:maximum number of exons to report [default: 0:999999]
    -v              produce verbose output
"""

import math
import sys
import docopt
import yaml

from zotmer.library.basics import kmersWithPos, kmersList, kmersLists, kmersWithPosList, kmersWithPosLists, murmer, render
from zotmer.library.file import openFile, readFasta, readFastq
from zotmer.library.hgvs import hg19ToRefSeq, refSeq2Hg19
from zotmer.library.misc import unionfind, uniq
from zotmer.library.reads import reads

verbose = False

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

            with openFile(self.home + "/" + h + ".fa.gz") as f:
                for (nm,seq) in readFasta(f):
                    self.prevAcc = acc
                    self.prevSeq = seq
                    break
        return self.prevSeq

class Summary(object):
    def __init__(self):
        self.n = 0
        self.s = 0
        self.s2 = 0

    def add(self, x):
        self.n += 1
        self.s += x
        self.s2 += x * x

    def __len__(self):
        return self.n

    def mean(self):
        return float(self.s)/float(self.n)

    def var(self):
        m = self.mean()
        return float(self.s2)/float(self.n) - m*m

    def sd(self):
        return math.sqrt(self.var())

codeDict = { 'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0 }
base64Idx = { }
base64Jdx = [None for i in range(64)]
for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
    i = len(base64Idx)
    base64Idx[c] = i
    base64Jdx[i] = c
for c in 'abcdefghijklmnopqrstuvwxyz':
    i = len(base64Idx)
    base64Idx[c] = i
    base64Jdx[i] = c
for c in '0123456789+/':
    i = len(base64Idx)
    base64Idx[c] = i
    base64Jdx[i] = c
assert len(base64Idx) == 64
assert len(base64Jdx) == 64

def compressInt6(j, r):
    x = []
    while True:
        x.append(base64Jdx[j & 63])
        j = j >> 6
        if j == 0:
            break
    n = len(x)
    if n > 1:
        r.append('.'*(n-1))
    r += x[::-1]

def compressInt12(j, r):
    x = []
    while True:
        x.append(base64Jdx[j & 63])
        j = j >> 6
        x.append(base64Jdx[j & 63])
        j = j >> 6
        if j == 0:
            break
    n = len(x) // 2
    if n > 1:
        r.append('.'*(n-1))
    r += x[::-1]

def compressRead(seq):
    global codeDict
    r = []
    compressInt12(len(seq), r)
    r.append('=')
    q = []
    b = 0
    c = 0
    i = 0
    for v in seq:
        c = (c << 2) | codeDict[v]
        b += 2
        if b == 6:
            r.append(base64Jdx[c])
            c = 0
            b = 0
        if v == 'N':
            q.append(i)
        i += 1

    if b != 0:
        r.append(base64Jdx[c])
    w = []
    if len(q) > 0:
        r += '='
        p = 0
        for i in q:
            compressInt6(i - p, r)
            p = i

    r = ''.join(r)
    print len(seq), len(r), seq, r
    return r

def prepareBedFileGene(geneListFn, refGeneFn, outFn):
    genes = {}
    with openFile(geneListFn) as f:
        for l in f:
            t = l.split()
            genes[t[0]] = set([])

    loci = {}
    strands = {}
    with openFile(refGeneFn) as f:
        for l in f:
            t = l.split()
            tx = t[1]
            ch = t[2]
            st = t[3]
            g = t[12]

            # Ignore chr6_cox_hap2 and friends.
            if '_' in ch:
                continue

            if g not in genes:
                continue

            if g in strands:
                if (ch,st) != strands[g]:
                    print >> sys.stderr, 'gene tx seen in multiple chromosomes: %s [%s] (%s%s, %s%s)' % (g, tx, strands[g][0], strands[g][1], ch, st)
                    continue
            else:
                strands[g] = (ch,st)

            ss = [int(s) for s in t[9].split(',') if len(s) > 0]
            ee = [int(e) for e in t[10].split(',') if len(e) > 0]
            exs = zip(ss, ee)
            for i in range(len(exs)):
                genes[g].add((exs[i][0], exs[i][1]))

    for g in sorted(genes.keys()):
            if len(genes[g]) == 0:
                print >> sys.stderr, 'no annotated transcripts for gene: %s' % (g,)
                continue
        
            (ch, st) = strands[g]
            exs0 = sorted(genes[g])

            exs1 = [exs0[0]]
            for i in range(1, len(exs0)):
                if exs1[-1][1] >= exs0[i][0]:
                    exs1[-1] = (min(exs1[-1][0], exs0[i][0]), max(exs1[-1][1], exs0[i][1]))
                else:
                    exs1.append(exs0[i])

            if len(exs0) != len(exs1):
                genes[g] = exs1
                #print >> sys.stderr, 'merged exons for gene %s (%d -> %d)' % (g, len(exs0), len(exs1))

    with openFile(outFn, 'w') as out:
        for g in sorted(genes.keys()):
            if len(genes[g]) == 0:
                continue

            (ch, st) = strands[g]
            exs = sorted(genes[g])
            if st == '-':
                exs = exs[::-1]

            for i in range(len(exs)):
                j = i + 1
                print >> out, '%s\t%d\t%d\t%s/%02d\t%s' % (ch, exs[i][0], exs[i][1], g, j, st)

def prepareBedFileGeneTx(geneListFn, refGeneFn, outFn):
    genes = set([])
    transcripts = {}
    with openFile(geneListFn) as f:
        for l in f:
            t = l.split()
            genes.add(t[0])
            transcripts[t[1]] = t[0]

    found = {}
    with openFile(refGeneFn) as f, openFile(outFn, 'w') as out:
        for l in f:
            t = l.split()
            tx = t[1]
            ch = t[2]
            st = t[3]
            g = t[12]
            if tx not in transcripts:
                continue
            found[tx] = g
            ss = [int(s) for s in t[9].split(',') if len(s) > 0]
            ee = [int(e) for e in t[10].split(',') if len(e) > 0]
            exs = zip(ss, ee)
            if st == '-':
                exs = exs[::-1]
            for i in range(len(exs)):
                j = i + 1
                print >> out, '%s\t%d\t%d\t%s/%02d\t%s' % (ch, exs[i][0], exs[i][1], g, j, st)

    for tx in sorted(set(transcripts.keys()) - set(found.keys())):
        print >> sys.stderr, 'transcipt not found: %s\t(%s)' % (tx, transcripts[tx])

rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
def revComp(seq):
    global rcDict
    r = []
    for c in seq:
        c = rcDict.get(c, c)
        r.append(c)
    return ''.join(r[::-1])

def indexBedFiles(bedFn, sf):
    global verbose

    idx = {}
    with openFile(bedFn) as f:
        for l in f:
            t = l.split()
            ch = t[0]
            s = int(t[1])
            e = int(t[2])
            gex = t[3]
            st = t[4]
            if ch not in idx:
                idx[ch] = []
            g, ex = gex.split('/')
            idx[ch].append((s, e, st, g, ex))
    for ch in sorted(idx.keys()):
        if verbose:
            print >> sys.stderr, 'processing %s' % (ch, )

        seq = sf[ch]
        idx[ch].sort()
        for (s, e, st, g, ex) in idx[ch]:
            exSeq = seq[s:e].upper()
            if st == '-':
                exSeq = revComp(exSeq)
            itm = {}
            itm['gene'] = g
            itm['exon'] = ex
            itm['chr'] = ch
            itm['st'] = s
            itm['en'] = e
            itm['strand'] = st
            itm['seq'] = exSeq
            yield itm

def sig(xs):
    h = 0
    for x in xs:
        h = murmer(x, h)
    return h

class ExonIndex(object):
    def __init__(self, K, refList):
        self.K = K
        self.exonLengths = [0 for i in range(len(refList))]

        idx = {}
        for i in range(len(refList)):
            itm = refList[i]
            for x in kmersList(self.K, itm['seq'], False):
                if x not in idx:
                    idx[x] = set([])
                idx[x].add(i)
                self.exonLengths[i] += 1

        self.idxUpper = {}
        self.idxLower = {}
        for x in idx.iterkeys():
            k = tuple(sorted(idx[x]))
            h = sig(k)
            if h not in self.idxLower:
                self.idxLower[h] = k
            self.idxUpper[x] = h

    def getKeyHash(self, x):
        if x in self.idxUpper:
            return self.idxUpper[x]
        return None

    def readHash(self, xs):
        ys = set([])
        for x in xs:
            y = self.getKeyHash(x)
            if y is not None:
                ys.add(y)

        if len(ys) == 0:
            return None

        return (sig(ys), ys)

def recHits(idx, xps):
    xs = []
    for (x,p) in xps:
        if x not in idx:
            continue
        xs.append(x)
    
    if len(xs) == 0:
        return ([], 0)

    S = set([])
    sdx = {}
    for x in xs:
        s = idx[x]
        if s not in sdx:
            sdx[s] = 0
        sdx[s] += 1
        S.add(s)
    S = sorted(S)
    print sig(S), S

    xSet = set(sdx.values())

    while True:

        u = unionfind()
        for i in range(len(S)):
            s0 = set(S[i])
            for j in range(i+1, len(S)):
                s1 = set(S[j])
                if len(s0 & s1) > 0:
                    u.union(i, j)
        grps = {}
        cnts = {}
        for i in range(len(S)):
            j = u.find(i)
            if j not in grps:
                grps[j] = set(S[i])
                cnts[j] = sdx[S[i]]
            else:
                grps[j] &= set(S[i])
                cnts[j] += sdx[S[i]]

        res = {}
        for j in grps.keys():
            if len(grps[j]) == 0:
                continue
            res[tuple(sorted(grps[j]))] = cnts[j]

        if len(res) > 0:
            return (res, len(xs))

        xMin = min(xSet)
        xSet.remove(xMin)

        if len(xSet) == 0:
            return (res, len(xs))

        xMin = min(xSet)
        S = [s for s in S if sdx[s] >= xMin]
        first = False

def main(argv):
    global verbose

    opts = docopt.docopt(__doc__, argv)

    verbose = opts['-v']

    genomeDir = '.'
    if opts['-g']:
        genomeDir = opts['-g']
    sf = SequenceFactory(genomeDir)

    if opts['-P']:
        if opts['-t']:
            prepareBedFileGeneTx(opts['<gene-list>'], opts['<refgene>'], opts['<bedfile>'])
        else:
            prepareBedFileGene(opts['<gene-list>'], opts['<refgene>'], opts['<bedfile>'])
        return

    if opts['-X']:
        with openFile(opts['<index>'], 'w') as out:
            yaml.safe_dump_all(indexBedFiles(opts['<must-have>'], sf), out, default_flow_style=False)
        return

    K = int(opts['-k'])
    minGeneReads = int(opts['-M'])
    minExonReads = int(opts['-m'])
    minGeneRate = float(opts['-R'])
    minExonRate = float(opts['-r'])
    (minGeneCount, maxGeneCount) = map(int, opts['-Z'].split(':'))
    (minExonCount, maxExonCount) = map(int, opts['-z'].split(':'))

    with openFile(opts['<index>']) as f:
        ref = list(yaml.load_all(f, Loader=yaml.BaseLoader))

    if True:
        # Test the double-layer index
        idx = ExonIndex(K, ref)

        acc = {}
        toc = {}
        rn = 0
        for itm in reads(opts['<input>'], K=K, paired=True, reads=True, kmers=False, both=True, verbose=verbose):
            rn += 1
            (lhsFwd, lhsRev) = kmersLists(K, itm.reads[0][1])
            (rhsFwd, rhsRev) = kmersLists(K, itm.reads[1][1])
            xs0 = lhsFwd + rhsRev
            rh0 = idx.readHash(xs0)
            if rh0 is not None:
                (h0, ys0) = rh0
                if h0 not in acc:
                    acc[h0] = []
                    toc[h0] = ys0
                acc[h0].append((compressRead(itm.reads[0][1]), compressRead(itm.reads[1][1])))

            xs1 = lhsRev + rhsFwd
            rh1 = idx.readHash(xs1)
            if rh1 is not None:
                (h1, ys1) = rh1
                if h1 not in acc:
                    acc[h1] = []
                    toc[h1] = ys1
                acc[h1].append((compressRead(itm.reads[0][1]), compressRead(itm.reads[1][1])))

        nx = 0
        for h in sorted(acc.keys()):
            for (x,c) in sorted(acc[h].items()):
                nx += 1
                if c <= 1:
                    continue
                print '%016x\t%s\t%d' % (h, render(K, x), c)

        print >> sys.stderr, 'nx =', nx
        return

    if False:
        # Position index
        idx = {}
        for i in range(len(ref)):
            itm = ref[i]
            for (x,p) in kmersWithPosList(K, itm['seq'], False):
                p -= 1
                if x not in idx:
                    idx[x] = []
                idx[x].append((i,p))

    if True:
        # Exon tuple index
        idx = {}
        lens = [0 for i in range(len(ref))]
        for i in range(len(ref)):
            itm = ref[i]
            for (x,p) in kmersWithPosList(K, itm['seq'], False):
                if x not in idx:
                    idx[x] = set([])
                idx[x].add(i)
                lens[i] += 1
        for x in idx.iterkeys():
            idx[x] = tuple(sorted(idx[x]))

    if opts['-T']:
        ak = {}
        for x in sorted(idx.iterkeys()):
            if len(idx[x]) == 1:
                continue
            xStr = render(K, x)
            ak[xStr] = []
            for i in idx[x]:
                itm = ref[i]
                k = '%s/%s' % (itm['gene'], itm['exon'])
                ak[xStr].append(k)
            ak[xStr].sort()
        rep = {}
        rep['aliasing-within'] = ak
        chrs = set([])
        for i in range(len(ref)):
            itm = ref[i]
            chrs.add(itm['chr'])
        counts = [0 for i in range(len(ref))]
        for ch in sorted(chrs):
            if verbose:
                print >> sys.stderr, 'processing %s' % (ch, )
            seq = sf[ch]
            for (x,p) in kmersWithPos(K, seq, True):
                if x not in idx:
                    continue
                for i in idx[x]:
                    counts[i] += 1
        gk = {}
        for i in range(len(ref)):
            if lens[i] == counts[i]:
                continue
            itm = ref[i]
            k = '%s/%s' % (itm['gene'], itm['exon'])
            gk[k] = {'indexed':lens[i], 'genomic':counts[i]}
        rep['aliasing-genomic'] = gk
        yaml.safe_dump(rep, sys.stdout, default_flow_style=False)
        return

    acc = {}
    rn = 0
    hitStats = Summary()
    hitHist = [0 for i in range(1000)]
    for itm in reads(opts['<input>'], K=K, paired=True, reads=True, kmers=False, both=True, verbose=verbose):
        rn += 1
        (lhsFwd, lhsRev) = kmersWithPosLists(K, itm.reads[0][1])
        (rhsFwd, rhsRev) = kmersWithPosLists(K, itm.reads[1][1])
        (hits0, hitCount0) = recHits(idx, lhsFwd + rhsRev)
        (hits1, hitCount1) = recHits(idx, lhsRev + rhsFwd)
        if len(hits0) > 0:
            k = tuple(sorted(hits0.keys()))
            v = sum(hits0.values())
            if k not in acc:
                acc[k] = [0, 0]
            acc[k][0] += 1
            acc[k][1] += v
            hitStats.add(hitCount0)
            hitHist[hitCount0] += 1

        if len(hits1) > 0:
            k = tuple(sorted(hits1.keys()))
            v = sum(hits1.values())
            if k not in acc:
                acc[k] = [0, 0]
            acc[k][0] += 1
            acc[k][1] += v
            hitStats.add(hitCount1)
            hitHist[hitCount1] += 1
    
    if verbose:
        print >> sys.stderr, 'total read hits: %d' % (len(hitStats), )
        print >> sys.stderr, 'total hits per read: %g (%g)' % (hitStats.mean(), hitStats.sd())
        print >> sys.stderr, 'total reads: %d' % (rn, )
        for i in range(len(hitHist)):
            if hitHist[i] > 0:
                print >> sys.stderr, '\t%d\t%d' % (i, hitHist[i])

    def gex(s) :
        r = []
        for n in s:
            itm = ref[n]
            r.append('%s/%s' % (itm['gene'], itm['exon']))
        return '|'.join(r)

    def fmtKey(k):
        nex = len(k)
        gx = set([])
        kStrParts = []
        for s in k:
            kStrParts.append(gex(s))
            gx |= set([ref[i]['gene'] for i in s])
        kStr = '--'.join(sorted(kStrParts))
        return (nex, gx, kStr)

    gxCounts = {}
    for k in acc.keys():
        gx = set([])
        ex = set([])
        for s in k:
            gx |= set([ref[i]['gene'] for i in s])
            ex |= set(s)
        gx = tuple(sorted(gx))
        if gx not in gxCounts:
            gxCounts[gx] = [0, 0]
        gxCounts[gx][0] += acc[k][0]
        gxCounts[gx][1] += acc[k][1]

    hdr = ['numReads', 'numKmers', 'kmersPerRead']
    hdr += ['ggNumReads', 'ggNumKmers', 'ggKmersPerRead']
    hdr += ['numExons', 'numGenes', 'geneGroup', 'exonGroup']
    print '\t'.join(hdr)
    for k in acc.keys():
        (nex, gx, kStr) = fmtKey(k)
        gx = tuple(sorted(gx))
        if len(gx) < minGeneCount or len(gx) > maxGeneCount:
            continue
        if len(ex) < minExonCount or len(ex) > maxExonCount:
            continue
        if gxCounts[gx][0] < minGeneReads:
            continue
        if acc[k][0] < minExonReads:
            continue
        gxRate = float(gxCounts[gx][1])/float(gxCounts[gx][0])
        if gxRate < minGeneRate:
            continue
        exRate = float(acc[k][1])/float(acc[k][0])
        if exRate < minExonRate:
            continue
        gxStr = ':'.join(gx)

        print '%d\t%d\t%g\t%d\t%d\t%g\t%d\t%d\t%s\t%s' % (acc[k][0], acc[k][1], exRate, gxCounts[gx][0], gxCounts[gx][1], gxRate, nex, len(gx), gxStr, kStr)

if __name__ == '__main__':
    main(sys.argv[1:])
