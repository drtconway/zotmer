"""
Usage:
    zot sample [-DS SEED] <probability> <output> <input>

Options:
    -D          use deterministic sampling
    -S SEED     use the given seed for the sampling
"""

from pykmer.basics import fnv
from pykmer.container import probe
import pykmer.kset as kset
import pykmer.kfset as kfset

import docopt
import random
import sys

def sampleR(p, xs):
    for x in xs:
        if random.random() < p:
            yield x

def sampleD(p, s, xs, accessor):
    M = 0xFFFFFFFFFF
    for x in xs:
        y = accessor(x)
        u = float(fnv(y, s)&M)/float(M)
        if u < p:
            yield x

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    p = float(opts['<probability>'])
    inp = opts['<input>']
    out = opts['<output>']
    (m, _) = probe(inp)
    if opts['-D'] is None:
        if opts['-S'] is not None:
            S = long(opts['-S'])
            random.seed(S)
        if m['type'] == 'k-mer set':
            (m, xs) = kset.read(inp)
            K = m['K']
            kset.write(K, sampleR(p, xs), out, m)
        else:
            (m, xs) = kfset.read(inp)
            K = m['K']
            kfset.write(K, sampleR(p, xs), out, m)
    else:
        S = 0
        if opts['-S'] is not None:
            S = long(opts['-S'])
        if m['type'] == 'k-mer set':
            (m, xs) = kset.read(inp)
            K = m['K']
            kset.write(K, sampleD(p, S, xs, lambda x: x), out, m)
        else:
            (m, xs) = kfset.read(inp)
            K = m['K']
            kfset.write(K, sampleD(p, S, xs, lambda x: x[0]), out, m)

if __name__ == '__main__':
    main(sys.argv[1:])
