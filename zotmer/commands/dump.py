"""
Usage:
    zot dump <input>
"""

from pykmer.basics import render
from pykmer.container import container
from pykmer.container.std import readKmers, readKmersAndCounts

import docopt
import sys

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    inp = opts['<input>']
    with container(inp, 'r') as z:
        K = z.meta['K']
        if 'kmers' not in z.meta:
            print >> sys.stderr, 'cannot dump "%s" as it contains no k-mers' % (inp,)
            return
        if 'counts' in z.meta:
            xs = readKmersAndCounts(z)
            for (x, c) in xs:
                print '%s\t%d' % (render(K, x), c)
        else:
            xs = readKmers(z)
            for x in xs:
                print render(K, x)

if __name__ == '__main__':
    main(sys.argv[1:])
