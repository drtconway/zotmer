"""
zot hgvs-find - search for known variants in read data.

Usage:
    zot fizzle -P [options] <gene-tx-list> <refgene> <bedfile>
    zot fizzle -X [options] <index> <must-have> [<may-have>]
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

from zotmer.library.basics import kmersWithPosList, render
from zotmer.library.file import openFile, readFasta, readFastq
from zotmer.library.hgvs import hg19ToRefSeq, refSeq2Hg19

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
                print >> out, '%s\t%d\t%d\t%s/ex%02d\t%s' % (ch, exs[i][0], exs[i][1], g, j, st)

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

def indexBedFiles(fnAndMust, sf):
    global verbose

    for (fn,mustHave) in fnAndMust:
        idx = {}
        with openFile(fn) as f:
            for l in f:
                t = l.split()
                ch = t[0]
                s = int(t[1])
                e = int(t[2])
                gex = t[3]
                st = t[4]
                if ch not in idx:
                    idx[ch] = []
                idx[ch].append((s, e, st, gex))
        for ch in sorted(idx.keys()):
            if verbose:
                print >> sys.stderr, 'processing %s' % (ch, )

            seq = sf[ch]
            idx[ch].sort()
            for (s, e, st, gex) in idx[ch]:
                exSeq = seq[s:e].upper()
                if st == '-':
                    exSeq = revComp(exSeq)
                itm = {}
                itm['exon'] = gex
                itm['reqd'] = mustHave
                itm['chr'] = ch
                itm['st'] = s
                itm['en'] = e
                itm['strand'] = st
                itm['seq'] = exSeq
                yield itm

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
        inputs = [(opts['<must-have>'],True)]
        if opts['<may-have>']:
            inputs.append((opts['<may-have>'],False))
        with openFile(opts['<index>'], 'w') as out:
            yaml.safe_dump_all(indexBedFiles(inputs, sf), out)
        return

    K = int(opts['-k'])

    with openFile(opts['<index>']) as f:
        ref = list(yaml.load_all(f, Loader=yaml.BaseLoader))

    idx = {}
    for i in range(len(ref)):
        itm = ref[i]
        for (x,p) in kmersWithPosList(K, itm['seq'], False):
            p -= 1
            if x not in idx:
                idx[x] = []
            idx[x].append((i,p))
    gex = {}
    for px in idx.itervalues():
        l = len(px)
        if l <= 10:
            continue
        for (i,p) in px:
            v = ref[i]['exon']
            if v not in gex:
                gex[v] = 0
            gex[v] += 1
    print sorted(gex.items())

if __name__ == '__main__':
    main(sys.argv[1:])
