"""
Usage:
    zot hist <input>...

Options:
    -u              update the input container to include the histogram
"""

from zotmer.library.kmers import kmers

import docopt
import sys

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    for inp in opts['<input>']:
        with kmers(inp, 'r') as z:
            if 'hist' in z.meta:
                h = z.meta['hist'].items()
                h = [(int(f), c) for (f, c) in h]
                h.sort()
                for (f,c) in h:
                    print '%s\t%d\t%d' % (inp, f, c)

if __name__ == '__main__':
    main(sys.argv[1:])
