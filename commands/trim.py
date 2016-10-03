from base import Cmd

import pykmer.kfset as kfset
from pykmer.container import probe

def trim(xs, c):
    for (x,f) in xs:
        if f >= c:
            yield(x, f)

class Trim(Cmd):
    def run(self, opts):
        inp = opts['<kmers>']
        out = opts['<output>']
        c = int(opts['<cutoff>'])
        (m, xs) = kfset.read(inp)
        K = m['K']
        kfset.write(K, trim(xs, c), out, m)

def add(cmds):
    cmds['trim'] = Trim()

