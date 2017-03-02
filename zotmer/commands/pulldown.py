"""
Usage:
    zot pulldown [-p] [-U FASTA] <baits> <output> <input>...

Pull down reads (or pairs) that have homology to bait sequences.
The output is the name of a zipfile in to which the reads for each
of the baits will be placed.

Options:
    -p          Treat reads as paired end. Requires an even number of inputs.
    -U FASTA    A FASTA file of sequences to "pushup".
"""

import os
import sys
import zipfile

import docopt

from pykmer.basics import kmersList, render
from pykmer.file import openFile, readFasta, readFastq, tmpfile

from zotmer.library.kmers import kmers
from zotmer.library.files import writeKmers

def pairs(xs):
    assert len(xs) & 1 == 0
    i = 0
    while i + 1 < len(xs):
        yield (xs[i], xs[i + 1])
        i += 2

def both(xs, ys):
    while True:
        try:
            x = xs.next()
            y = ys.next()
            yield(x, y)
        except StopIteration:
            return

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = 25

    nms = []
    idx = {}
    for (nm, seq) in readFasta(openFile(opts['<baits>'])):
        n = len(nms)
        nms.append(nm)
        for x in kmersList(K, seq, True):
            if x not in idx:
                idx[x] = set([])
            idx[x].add(n)

    for x in idx.keys():
        idx[x] = list(idx[x])
        idx[x].sort()

    anti = set([])
    if opts['-U']:
        with openFile(opts['-U']) as f:
            for (nm, seq) in readFasta(f):
                for x in kmersList(K, seq, True):
                    anti.add(x)

    rn = 0
    if opts['-p']:

        hist = {}
        for (fn1, fn2) in pairs(opts['<input>']):
            tmps = [(tmpfile('_1.fastq'), tmpfile('_2.fastq')) for i in xrange(len(nms))]
            cache = [[[], []] for i in xrange(len(nms))]
            counts = [0 for i in xrange(len(nms))]
            with openFile(fn1) as f1, openFile(fn2) as f2:
                for fq1, fq2 in both(readFastq(f1), readFastq(f2)):
                    hits = set([])
                    pushup = False
                    for x in kmersList(K, fq1[1]):
                        if x in anti:
                            pushup = True
                            break
                        for i in idx.get(x, []):
                            hits.add(i)
                    for x in kmersList(K, fq2[1]):
                        if x in anti:
                            pushup = True
                            break
                        for i in idx.get(x, []):
                            hits.add(i)
                    if pushup:
                        continue
                    n = len(hits)
                    hist[n] = 1 + hist.get(n, 0)
                    for i in hits:
                        counts[i] += 1
                        cache[i][0].append(fq1)
                        cache[i][1].append(fq2)
                        if len(cache[i][0]) >= 1024:
                            with open(tmps[i][0], 'a') as f:
                                for rd in cache[i][0]:
                                    print >> f, rd[0]
                                    print >> f, rd[1]
                                    print >> f, rd[2]
                                    print >> f, rd[3]
                            with open(tmps[i][1], 'a') as f:
                                for rd in cache[i][1]:
                                    print >> f, rd[0]
                                    print >> f, rd[1]
                                    print >> f, rd[2]
                                    print >> f, rd[3]
                            cache[i][0] = []
                            cache[i][1] = []
            for i in xrange(len(cache)):
                if len(cache[i][0]) > 0:
                    with open(tmps[i][0], 'a') as f:
                        for rd in cache[i][0]:
                            print >> f, rd[0]
                            print >> f, rd[1]
                            print >> f, rd[2]
                            print >> f, rd[3]
                    with open(tmps[i][1], 'a') as f:
                        for rd in cache[i][1]:
                            print >> f, rd[0]
                            print >> f, rd[1]
                            print >> f, rd[2]
                            print >> f, rd[3]
                    cache[i][0] = []
                    cache[i][1] = []
            with zipfile.ZipFile(opts['<output>'], 'w', zipfile.ZIP_DEFLATED) as z:
                for i in xrange(len(nms)):
                    if counts[i] > 0:
                        pth = '/'.join(nms[i].split())
                        z.write(tmps[i][0], pth + '/' + fn1)
                        os.remove(tmps[i][0])
                        z.write(tmps[i][1], pth + '/' + fn2)
                        os.remove(tmps[i][1])
        hist = hist.items()
        hist.sort()
        for (n,f) in hist:
            print '%d\t%d' % (n, f)
    else:
        raise "not implemented"

if __name__ == '__main__':

    main(sys.argv[1:])
