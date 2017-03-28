"""
Usage:
    zot sample [-DS SEED] [-P PROBABILITY] <output> <input>

Options:
    -D              use deterministic sampling
    -P PROBABILITY  the proportion of samples to include in the output.
                    default: 0.01
    -S SEED         use the given seed for the sampling
"""

import random
import sys

import docopt

from pykmer.basics import murmer
from zotmer.library.kmers import kmers
from zotmer.library.files import readKmersAndCounts, writeKmersAndCounts

def sampleR(p, xs, h):
    for x in xs:
        if random.random() < p:
            h[x[1]] = 1 + h.get(x[1], 0)
            yield x

def sampleD(p, s, xs, h):
    M = 0xFFFFFFFFFF
    for x in xs:
        y = x[0]
        u = float(murmer(y, s)&M)/float(M)
        if u < p:
            h[x[1]] = 1 + h.get(x[1], 0)
            yield x

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    p = 0.01
    if opts['-P'] is not None:
        p = float(opts['-P'])
    inp = opts['<input>']
    out = opts['<output>']
    with kmers(out, 'w') as z:
        h = {}
        with kmers(inp, 'r') as z0:
            K = z0.meta['K']
            z.meta = z0.meta.copy()
            del z.meta['kmers']
            del z.meta['counts']
            xs = readKmersAndCounts(z0)
            S = 0
            if opts['-D'] is None:
                if opts['-S']:
                    S = long(opts['-S'])
                    random.seed(S)
                writeKmersAndCounts(z, sampleR(p, xs, h))
            else:
                if opts['-S']:
                    S = long(opts['-S'])
                writeKmersAndCounts(z, sampleD(p, S, xs, h))
        z.meta['K'] = K
        z.meta['kmers'] = 'kmers'
        z.meta['counts'] = 'counts'
        z.meta['hist'] = h

if __name__ == '__main__':
    main(sys.argv[1:])
