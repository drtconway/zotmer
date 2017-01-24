"""
Usage:
    zot mlst [-XK K] <alleles> <input>...

Determine which, out of a FASTA database of alleles exist in the
input set of k-mers.

Options:
    -X          Create an index
    -K K        If creating an index, use this value of K (default 27).
"""

import array
import cPickle
import sys

import docopt

from pykmer.container import container
from pykmer.container.std import readKmers
from pykmer.index import buildIndex, index

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    if opts['-X']:
        K = 27
        if opts['-K']:
            K = int(opts['-K'])

        buildIndex(K, opts['<input>'], opts['<alleles>'])

        return

    idx = index(opts['<alleles>'])

    for inp in opts['<input>']:
        with container(inp, 'r') as z:
            K0 = z.meta['K']
            if K0 != idx.K:
                print >> sys.stderr, 'Input "%d" has different k to index' % (inp,)
                sys.exit(1)
            xs = readKmers(z)

            cs = [idx.lens[i] for i in xrange(len(idx.lens))]

            for x in xs:
                for j in idx[x]:
                    cs[j] -= 1

            for j in xrange(len(idx.lens)):
                assert cs[j] >= 0
                if cs[j] == 0:
                    print '%s\t%d\t%s' % (inp, j, idx.names[j])

if __name__ == '__main__':
    main(sys.argv[1:])
