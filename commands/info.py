from base import Cmd

from pykmer.container import probe

class Info(Cmd):
    def run(self, opts):
        for inp in opts['<input>']:
            (m, _) = probe(inp)
            print m

def add(cmds):
    cmds['info'] = Info()

