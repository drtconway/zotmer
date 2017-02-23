"""
"""

import json

from pykmer.container.casket import casket

class kmers(casket):
    def __init__(self, fn, mode):
        super(kmers, self).__init__(fn, mode)

        self.meta = {}
        if mode == 'r':
            f = self.open('__meta__')
            s = f.read()
            self.meta = json.loads(s)

    def close(self):
        if self.mode == 'w':
            self.add_content('__meta__', json.dumps(self.meta))
        super(kmers, self).close()

