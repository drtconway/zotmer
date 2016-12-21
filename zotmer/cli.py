"""
Usage:
    zot [options] <command> [<args>...]

options:
    --version   print version information
"""

import docopt
import pkgutil
import os
import sys
import importlib

import commands

def main():
    #if 'ZOTMER_PLUGINS' in os.environ:
    #    nms = os.environ['ZOTMER_PLUGINS'].split(':')
    #    for nm in nms:
    #        m = __import__(nm)
    #        m.add(cmds)
    args = docopt.docopt(__doc__, version='Zotmer k-mer toolkit 0.1')
    modname = commands.__name__ + '.' + args['<command>']
    try:
        argv = [args['<command>']] + args['<args>']
        m = importlib.import_module(modname)
        return m.main(argv)
    except ImportError:
        print >> sys.stderr, "unable to load command `%s', use `zot help` for help." % (args['<command>'], )
        return 1


if __name__ == '__main__':
    main()

