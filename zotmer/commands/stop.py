"""
zot stop - search for stop codons

Usage:
    zot stop [options] <input>...

Options:
    -k K            k-mer size [default: 27]
    -r FILE         filter candidates against the given transcripts
    -v              produce verbose output
"""

import math
import sys

import docopt
import yaml

from pykmer.basics import kmer, kmersList, kmersWithPosList, render
from pykmer.file import openFile, readFasta, readFastq

stopStrs = ['TAA', 'TAG', 'TGA']
stops = [kmer(s) for s in stopStrs]

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = int(opts['-k'])
    K3 = K - 3

    ref = {}
    if opts['-r']:
        with openFile(opts['-r']) as f:
            for (nm,seq) in readFasta(f):
                for (x,p) in kmersWithPosList(K3, seq, False):
                    if (p % 3) == 1:
                        if x not in ref:
                            ref[x] = set([])
                        ref[x].add(nm)

        kx = {}
        for fn in opts['<input>']:
            with openFile(fn) as f:
                for rd in readFastq(f):
                    xs = kmersList(K, rd[1], True)
                    for x in xs:
                        # mask off low order codon.
                        c = x & 63
                        if c in stops:
                            y = x >> 6
                            if y in ref:
                                if x not in kx:
                                    kx[x] = 0
                                kx[x] += 1

        for (x,c) in kx.iteritems():
            y = x >> 6
            for nm in ref[y]:
                print '%s\t%s\t%d' % (render(K, x), nm, c)
    else:
        kx = {}
        for fn in opts['<input>']:
            with openFile(fn) as f:
                for rd in readFastq(f):
                    xs = kmersList(K, rd[1], True)
                    for x in xs:
                        # mask off low order codon.
                        c = x & 63
                        if c in stops:
                            if x not in kx:
                                kx[x] = 0
                            kx[x] += 1

        #for (x,c) in kx.iteritems():
        #    print '%s\t%d' % (render(K, x), c)

        h = {}
        for (x,c) in kx.iteritems():
            if c not in h:
                h[c] = 0
            h[c] += 1
        for (c,f) in h.iteritems():
            print '%d\t%d' % (c, f)

if __name__ == '__main__':
    main(sys.argv[1:])
