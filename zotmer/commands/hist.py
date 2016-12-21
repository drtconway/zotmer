"""
Usage:
    zot hist <input>...
"""

import pykmer.kfset as kfset

import docopt

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    h = {}
    for inp in opts['<input>']:
        (m, xs) = kfset.read(inp)
        for (_, f) in xs:
            c = h.get(f, 0)
            h[f] = c + 1
    h = h.items()
    h.sort()
    for (f,c) in h:
        print '%d\t%d' % (f, c)

if __name__ == '__main__':
    main(sys.argv[1:])
