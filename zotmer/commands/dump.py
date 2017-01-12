"""
Usage:
    zot dump <input>...
"""

from pykmer.basics import render
from pykmer.container import container
from pykmer.container.std import readKmers, readCounts, readKmersAndCounts

import docopt
import sys

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    for inp in opts['<input>']:
        with container(inp) as z:
            for nm in z.find(kind=['k-mers', 'k-mers and counts']):
                meta = manifest[nm]
                assert 'kind' in meta
                if meta['kind'] == 'k-mers':
                    (m, xs) = readKmers(z, nm)
                    K = m['K']
                    for x in xs:
                        print render(K, x)
                if meta['kind'] == 'k-mers and counts':
                    (m, xs) = readKmers(z, nm)
                    K = m['K']
                    for (x, c) in xs:
                        print '%s\t%d' % (render(K, x), c)

if __name__ == '__main__':
    main(sys.argv[1:])
