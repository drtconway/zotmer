"""
Usage:
    zot project <ref> <input>...

Project one or more inputs on to a reference set. For each k-mer in <ref>,
a whitespace separated 0 or a 1 is printed indicating whether that k-mer
was present in the input, with a separate line for each input k-mer set.
"""

from pykmer.adaptors import kf2k
from pykmer.container import probe
import pykmer.kset as kset
import pykmer.kfset as kfset

import array
import docopt
import sys

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    refFn = opts['<ref>']
    (meta, xs) = kset.read(refFn)
    xs = array.array('L', xs)
    Z = len(xs)

    K = meta['K']

    for inp in opts['<input>']:
        (m, _) = probe(inp)
        K0 = m['K']
        if K0 != K:
            print >> sys.stderr, "input %s has mismatched K (%d)" % (inp, K0)
            sys.exit(1)

        if m['type'] == 'k-mer set':
            (_, ys) = kset.read(inp)
        else:
            (_, ys) = kfset.read(inp)
            ys = kf2k(ys)

        i = 0
        row = []
        for y in ys:
            while i < Z and xs[i] < y:
                row.append('0')
                i += 1
            if i < Z and xs[i] == y:
                row.append('1')
                i += 1
        while i < Z:
            row.append('0')
        print '\t'.join(row)

if __name__ == '__main__':
    main(sys.argv[1:])
