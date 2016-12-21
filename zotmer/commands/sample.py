"""
Usage:
    zot sample [-DS SEED] <probability> <output> <input>

Options:
    -D          use deterministic sampling
    -S SEED     use the given seed for the sampling
"""

import pykmer.kset as kset
import pykmer.kfset as kfset
from pykmer.container import probe

import docopt
import random

def sample(p, xs):
    for x in xs:
        if random.random() < p:
            yield x

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    p = float(opts['<probability>'])
    if opts['-S'] is not None:
        S = long(opts['-S'])
        random.seed(S)
    inp = opts['<input>']
    out = opts['<output>']
    (m, _) = probe(inp)
    if m['type'] == 'k-mer set':
        (m, xs) = kset.read(inp)
        K = m['K']
        kset.write(K, sample(p, xs), out, m)
    else:
        (m, xs) = kfset.read(inp)
        K = m['K']
        kfset.write(K, sample(p, xs), out, m)

if __name__ == '__main__':
    main(sys.argv[1:])
