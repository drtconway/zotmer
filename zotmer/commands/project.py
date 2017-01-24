"""
Usage:
    zot project <ref> <output> <input>

Project one or more inputs on to a reference set. For each k-mer in <ref>,
a whitespace separated 0 or a 1 is printed indicating whether that k-mer
was present in the input, with a separate line for each input k-mer set.
"""

import array
import sys

import docopt

from pykmer.container import container
from pykmer.container.std import readKmers, writeKmers, readKmersAndCounts, writeKmersAndCounts

def project1(xs, ys):
    i = 0
    Z = len(xs)
    for y in ys:
        while i < Z and xs[i] < y:
            i += 1
        if i == Z:
            break
        if xs[i] == y:
            yield y
            i += 1

def project2(xs, ys):
    i = 0
    Z = len(xs)
    for y in ys:
        while i < Z and xs[i] < y[0]:
            i += 1
        if i == Z:
            break
        if xs[i] == y[0]:
            yield y
            i += 1

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    with container(opts['<ref>'], 'r') as z:
        K = z.meta['K']
        xs = array.array('L', readKmers(z))
    Z = len(xs)

    with container(opts['<input>'], 'r') as z0:
        K0 = z0.meta['K']
        if K0 != K:
            print >> sys.stderr, "mismatched K (%d)" % (K0, )
            sys.exit(1)

        with container(opts['<output>'], 'w') as z:
            if 'counts' in z0.meta:
                ys = readKmersAndCounts(z0)
                writeKmersAndCounts(K, project2(xs, ys), z)
            else:
                ys = readKmers(z0)
                writeKmers(K, project1(xs, ys), z)

if __name__ == '__main__':
    main(sys.argv[1:])
