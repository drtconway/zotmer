from base import Cmd

import pykmer.kfset as kfset

class Hist(Cmd):
    def run(self, opts):
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

def add(cmds):
    cmds['hist'] = Hist()

