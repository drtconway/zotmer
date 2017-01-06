"""
Usage:
    zot mlst [-XK K] <alleles> <input>...

Determine which, out of a FASTA database of alleles exist in the
input set of k-mers.

Options:
    -X          Create an index
    -K K        If creating an index, use this value of K (default 27).
"""

from pykmer.adaptors import kf2k
from pykmer.basics import kmers
from pykmer.container import probe
from pykmer.file import readFasta
from pykmer.misc import uniq
from pykmer.sparse import sparse
import pykmer.kfset as kfset
import pykmer.kset as kset

import array
import cPickle
import docopt
import sys

class index:
    def __init__(self):
        pass

    def create(self, K, alleles):
        self.K = K
        self.nms = []
        self.lens = array.array('I', [])

        kxs = []
        zs = []
        for (nm, seq) in alleles:
            xs = list(kmers(K, seq))
            self.nms.append(nm)
            kxs.append(xs)
            self.lens.append(len(xs))
            zs += xs

        zs.sort()
        self.zs = sparse(2*K, array.array('L', zs))

        self.cs = array.array('I', [0 for i in xrange(self.zs.count()+1)])
        for xs in kxs:
            for x in xs:
                r = self.zs.access(x)
                assert r is not None
                self.cs[r] += 1

        t = 0
        ts = []
        for i in xrange(len(self.cs)):
            t0 = t
            t += self.cs[i]
            self.cs[i] = t0
            ts.append(self.cs[i])

        self.idx = array.array('I', [0 for i in xrange(sum(self.lens))])
        for i in xrange(len(kxs)):
            xs = kxs[i]
            for x in xs:
                r = self.zs.access(x)
                self.idx[ts[r]] = i
                ts[r] += 1

    def save(self, nm):
        meta = {}
        meta['K'] = self.K
        meta['nms'] = self.nms
        meta['lens'] = len(self.lens)
        meta['zs'] = self.zs.count()
        meta['cs'] = len(self.cs)
        meta['idx'] = len(self.idx)
        with open(nm + ".meta", 'w') as f:
            cPickle.dump(meta, f)
        with open(nm + ".data", 'w') as f:
            self.lens.tofile(f)
            self.zs.xs.tofile(f)
            self.cs.tofile(f)
            self.idx.tofile(f)

    def load(self, nm):
        meta = {}
        with open(nm + ".meta") as f:
            meta = cPickle.load(f)

        self.K = meta['K']
        self.nms = meta['nms']

        with open(nm + ".data") as f:
            self.lens = array.array('I', [])
            self.lens.fromfile(f, meta['lens'])
            zs = array.array('L', [])
            zs.fromfile(f, meta['zs'])
            self.zs = sparse(2*self.K, zs)
            self.cs = array.array('I', [])
            self.cs.fromfile(f, meta['cs'])
            self.idx = array.array('I', [])
            self.idx.fromfile(f, meta['idx'])

    def init(self):
        self.seen = array.array('I', [self.lens[i] for i in xrange(len(self.lens))])

    def add(self, x):
        r = self.zs.access(x)
        if r is None:
            return
        for i in xrange(self.cs[r], self.cs[r+1]):
            a = self.idx[i]
            self.seen[a] -= 1

    def end(self):
        for i in xrange(len(self.lens)):
            if self.seen[i] == 0:
                yield (i, self.nms[i])

def inputs(fns):
    for fn in fns:
        for (nm, seq) in readFasta(open(fn)):
            yield (nm, seq)

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    if opts['-X']:
        K = 27
        if opts['-K']:
            K = int(opts['-K'])

        idx = index()
        idx.create(K, inputs(opts['<input>']))
        idx.save(opts['<alleles>'])

        return

    idx = index()
    idx.load(opts['<alleles>'])

    for inp in opts['<input>']:
        (m, _) = probe(inp)
        if m['type'] == 'k-mer set':
            (m, xs) = kset.read(inp)
            K = m['K']
        else:
            (m, xs) = kfset.read(inp)
            K = m['K']
            xs = kf2k(xs)

        idx.init()
        for x in xs:
            idx.add(x)
        for (n, nm) in idx.end():
            print '%s\t%d\t%s' % (inp, n, nm)

if __name__ == '__main__':
    main(sys.argv[1:])
