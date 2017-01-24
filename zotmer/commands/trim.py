"""
Usage:
    zot trim [-sc CUTOFF] <output> <input>

Options:
    -c CUTOFF   discard k-mers with frequency less than CUTOFF. A
                cutoff of 0 (the default) indicates that cutoff
                inference should be used.
    -s          force output to be a k-mer set
"""

import math
import sys

import docopt

from pykmer.container import container
from pykmer.container.std import readKmersAndCounts, writeKmersAndCounts

def sqr(x):
    return x * x

def gauss(x0, x, n):
    return math.exp(-sqr(x0 - x)/(2.0*n*n))

def smooth(h):
    h1 = {}
    for x0 in h.keys():
        num = 0
        den = 0
        for (xi,yi) in h.items():
            z = gauss(x0, xi, 1.5)
            num += z*yi
            den += z
        y0 = num/den
        h1[x0] = y0
    return h1

def infer(K, h):
    h1 = smooth(h)
    h1 = h1.items()
    h1.sort()
    v = h1[0]
    for v1 in h1[1:]:
        if v1[1] < v[1]:
            v = v1
        else:
            break
    return v[0]

def trim(xs, c):
    for (x,f) in xs:
        if f >= c:
            yield (x, f)

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    inp = opts['<input>']
    out = opts['<output>']

    c = 0
    if opts['-c'] is not None:
        c = int(opts['-c'])

    with container(inp, 'r') as z:
        K = z.meta['K']
        h = z.meta['hist']
        if c == 0:
            c = infer(K, h)
            print >> sys.stderr, 'inferred cutoff:', c
        xs = readKmersAndCounts(z)
        with container(out, 'w') as w:
            w.meta = z.meta.copy()
            del w.meta['kmers']
            del w.meta['counts']
            writeKmersAndCounts(K, trim(xs, c), w)

if __name__ == '__main__':
    main(sys.argv[1:])
