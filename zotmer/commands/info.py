"""
Usage:
    zot info <input>...
"""

from pykmer.container import probe

import docopt

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    for inp in opts['<input>']:
        (m, _) = probe(inp)
        print m

if __name__ == '__main__':
    main(sys.argv[1:])
