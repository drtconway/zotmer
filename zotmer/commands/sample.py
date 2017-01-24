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
from pykmer.container import container
from pykmer.container.std import readKmersAndCounts, writeKmersAndCounts

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
    with container(out, 'w') as z:
        h = {}
        with container(inp, 'r') as z0:
            K = z0.meta['K']
            z.meta = z0.meta.copy()
            del z.meta['kmers']
            del z.meta['counts']
            xs = readKmersAndCounts(z0)
            if opts['-D'] is None:
                if opts['-S'] is not None:
                    S = long(opts['-S'])
                    random.seed(S)
                    writeKmersAndCounts(K, sampleR(p, xs, h), z)
            else:
                S = 0
                if opts['-S'] is not None:
                    S = long(opts['-S'])
                    writeKmersAndCounts(K, sampleD(p, S, xs, h), z)
        z.meta['hist'] = h

if __name__ == '__main__':
    main(sys.argv[1:])
