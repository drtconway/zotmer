"""
zot hgvs-find - search for known variants in read data.

Usage:
    zot fizzle -P [options] <gene-tx-list> <refgene> <bedfile>
    zot fizzle -X [options] <index> <must-have>
    zot fizzle [options] <index> <input>...

Options:
    -P              prepare a bed file
    -X              index exons
    -g PATH         directory of FASTA reference sequences
    -k K            value of k to use [default: 25]
    -v              produce verbose output
"""

import sys
import docopt
import yaml

from zotmer.library.basics import kmersWithPosList, kmersWithPosLists, render
from zotmer.library.file import openFile, readFasta, readFastq
from zotmer.library.hgvs import hg19ToRefSeq, refSeq2Hg19
from zotmer.library.misc import unionfind
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

def prepareBedFile(geneListFn, refGeneFn, outFn):
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

def recHits(idx, xps):
    xs = []
    for (x,p) in xps:
        if x not in idx:
            continue
        xs.append(x)
    
    if len(xs) == 0:
        return []

    S = set([])
    sdx = {}
    for x in xs:
        s = idx[x]
        if s not in sdx:
            sdx[s] = 0
        sdx[s] += 1
        S.add(s)
    S = sorted(S)

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
            return res

        xMin = min(xSet)
        xSet.remove(xMin)

        if len(xSet) == 0:
            return res

        xMin = min(xSet)
        S = [s for s in S if sdx[s] >= xMin]
        first = False

    return res


def main(argv):
    global verbose

    opts = docopt.docopt(__doc__, argv)

    verbose = opts['-v']

    genomeDir = '.'
    if opts['-g']:
        genomeDir = opts['-g']
    sf = SequenceFactory(genomeDir)

    if opts['-P']:
        prepareBedFile(opts['<gene-tx-list>'], opts['<refgene>'], opts['<bedfile>'])
        return

    if opts['-X']:
        with openFile(opts['<index>'], 'w') as out:
            yaml.safe_dump_all(indexBedFiles(opts['<must-have>'], sf), out, default_flow_style=False)
        return

    K = int(opts['-k'])

    with openFile(opts['<index>']) as f:
        ref = list(yaml.load_all(f, Loader=yaml.BaseLoader))

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
        for i in range(len(ref)):
            itm = ref[i]
            for (x,p) in kmersWithPosList(K, itm['seq'], False):
                p -= 1
                if x not in idx:
                    idx[x] = set([])
                idx[x].add(i)
        for x in idx.iterkeys():
            idx[x] = tuple(sorted(idx[x]))

    if False:
        for x in sorted(idx.keys()):
            if len(idx[x]) > 1:
                xx = set([])
                for (i,p) in idx[x]:
                    itm = ref[i]
                    k = '%s/%s' % (itm['gene'], itm['exon'])
                    xx.add(k)
                xx = sorted(xx)
                print >> sys.stderr, 'ambiguous k-mer: %s\t%d\t%s' % (render(K, x), len(idx[x]), ', '.join(xx))
    if False:
        gex = {}
        for px in idx.itervalues():
            l = len(px)
            for (i,p) in px:
                v = ref[i]['gene'] + '/' + ref[i]['exon']
                if v not in gex:
                    gex[v] = {}
                if l not in gex[v]:
                    gex[v][l] = 0
                gex[v][l] += 1

        for v in sorted(gex.keys()):
            for (l,c) in sorted(gex[v].items()):
                print '%s\t%d\t%d' % (v, l, c)

    acc = {}
    rn = 0
    for itm in reads(opts['<input>'], K=K, paired=True, reads=True, kmers=False, both=True, verbose=verbose):
        rn += 1
        (lhsFwd, lhsRev) = kmersWithPosLists(K, itm.reads[0][1])
        (rhsFwd, rhsRev) = kmersWithPosLists(K, itm.reads[1][1])
        hits0 = recHits(idx, lhsFwd + rhsRev)
        hits1 = recHits(idx, lhsRev + rhsFwd)
        if len(hits0) > 0:
            k = tuple(sorted(hits0.keys()))
            v = sum(hits0.values())
            if k not in acc:
                acc[k] = [0, 0]
            acc[k][0] += 1
            acc[k][1] += v

        if len(hits1) > 0:
            k = tuple(sorted(hits1.keys()))
            v = sum(hits1.values())
            if k not in acc:
                acc[k] = [0, 0]
            acc[k][0] += 1
            acc[k][1] += v

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
        for s in k:
            gx |= set([ref[i]['gene'] for i in s])
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
        gxStr = ':'.join(gx)
        print '%d\t%d\t%g\t%d\t%d\t%g\t%d\t%d\t%s\t%s' % (acc[k][0], acc[k][1], float(acc[k][1])/float(acc[k][0]), gxCounts[gx][0], gxCounts[gx][1], float(gxCounts[gx][1])/float(gxCounts[gx][0]), nex, len(gx), gxStr, kStr)

if __name__ == '__main__':
    main(sys.argv[1:])
