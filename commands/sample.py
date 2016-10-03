from base import Cmd

import pykmer.kset as kset
import pykmer.kfset as kfset
from pykmer.container import probe

import random

def sample(p, xs):
    for x in xs:
        if random.random() < p:
            yield x

class Sample(Cmd):
    def run(self, opts):
        p = float(opts['<prob>'])
        if opts['-S'] is not None:
            S = long(opts['-S'])
            random.seed(S)
        inp = opts['<kmers>']
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

def add(cmds):
    cmds['sample'] = Sample()

