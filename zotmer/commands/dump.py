from base import Cmd

from pykmer.basics import render
import pykmer.kset as kset
import pykmer.kfset as kfset
from pykmer.container import probe

class Dump(Cmd):
    def run(self, opts):
        for inp in opts['<input>']:
            (m, _) = probe(inp)
            if m['type'] == 'k-mer set':
                (m, xs) = kset.read(inp)
                K = m['K']
                for x in xs:
                    print render(K, x)
            else:
                (m, xs) = kfset.read(inp)
                K = m['K']
                for (x,f) in xs:
                    print '%s\t%d' % (render(K, x), f)

def add(cmds):
    cmds['dump'] = Dump()

