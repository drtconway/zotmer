"""
Usage:
    zot info <input>...
"""

from zotmer.library.kmers import kmers

import docopt

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    for inp in opts['<input>']:
        with kmers(inp, 'r') as z:
            itms = z.meta.items()
            itms.sort()
            for (k, v) in itms:
                print k, v

if __name__ == '__main__':
    main(sys.argv[1:])
