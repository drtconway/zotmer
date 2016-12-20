"""Almost python k-mer toolkit.

Usage:
    almost kmerize [-sm MEM] <k> <out> <input>...
    almost merge <output> <input>...
    almost trim <cutoff> <output> <kmers>
    almost dist [-M measure]... <k> <input>...
    almost hist <input>...
    almost info <input>...
    almost scan [-AX] <genes> <input>...
    almost vars [-r ref] <input>...
    almost sample [-S SEED] <prob> <output> <kmers>
    almost dump <input>...

Options:
    -A          print alleles
    -c          cumulative scanning: combine LCPs across ksets
    -m MEM      in-memory buffer size
    -M measure  use "measure" for the distance between k-mer frequency sets.
                Use "-M list" to get a list of available measures.
    -r ref      use a reference sequence/k-mers
    -S SEED     set the seed for the RNG
    -s          generate a k-mer set rather than a k-mer frequency set
    -X          create an index
"""

from inspect import getmembers, isclass
from docopt import docopt

import commands

if __name__ == '__main__':
    arguments = docopt(__doc__, version='Almost k-mer toolkit 0.1')
    for (k,v) in arguments.items():
        if k in commands.cmds and v:
            commands.cmds[k].run(arguments)
