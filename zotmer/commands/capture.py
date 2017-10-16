
"""
zot capture - capture reads matching bait sequences

Usage:
    zot capture [options] <sequences> <input>...

Options:
    -k K            value of K to use - must be even. [default: 24]
    -P PREFIX       prefix for output files [default: .]
    -p              treat inputs as paired reads
    -v              produce verbose output
"""

import sys
import zipfile

import docopt
import tqdm
import yaml

from pykmer.basics import kmersList
from pykmer.file import openFile, readFasta, readFastq
from zotmer.library.reads import reads

class ReadCache(object):
    def __init__(self, pfx, nm, paired):
        self.B = 1024
        self.nr = 0
        if paired:
            nm1 = '%s/%s_1.fastq' % (pfx, nm)
            nm2 = '%s/%s_2.fastq' % (pfx, nm)
            self.names = [nm1, nm2]
            self.buffers = [[], []]
            self.N = 2
        else:
            nm = '%s/%s.fastq' % (pfx, nm)
            self.names = [nm]
            self.buffers = [[]]
            self.N = 1

    def add(self, rds):
        assert len(rds) == self.N
        self.nr += 1
        for i in xrange(self.N):
            self.buffers[i].append(rds[i])
        if len(self.buffers[0]) >= self.B:
            self.flush()

    def flush(self):
        for i in range(self.N):
            if len(self.buffers[i]) == 0:
                continue
            with open(self.names[i], 'a') as f:
                for rd in self.buffers[i]:
                    print >> f, rd[0]
                    print >> f, rd[1]
                    print >> f, rd[2]
                    print >> f, rd[3]
            self.buffers[i] = []

    def end(self):
        self.flush()
        print >> sys.stderr, '%s: %d' % (self.names[0], self.nr)

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = int(opts['-k'])
    if (K & 1) != 0:
        print >> sys.stderr, "K must be even."
        return

    paired = opts['-p']

    verbose = opts['-v']

    names = []
    seqs = []
    baits = {}
    with openFile(opts['<sequences>']) as f:
        for (nm, seq) in readFasta(f):
            n = len(names)
            names.append(nm)
            seqs.append(seq)
            for x in kmersList(K, seq, True):
                if x not in baits:
                    baits[x] = set([])
                baits[x].add(n)

    N = len(names)

    caches = [ReadCache(opts['-P'], names[n], paired)  for n in range(N)]

    nr = 0
    nh = 0
    for itm in reads(opts['<input>'], reads=True, kmers=True, fwdOnly=True, paired=paired, verbose=verbose):
        nr += 1
        # itm[0] == reads
        # item[0][i] == read end i
        # itm[1] == kmers
        # itm[1][i] == kmers end i
        # itm[1][i][s] == kmers end i, strand s
        E = len(itm[1])
        hits = set([])
        for i in xrange(E):
            fwd = itm[1][i][0]
            for x in fwd:
                if x in baits:
                    hits |= baits[x]
        for n in hits:
            caches[n].add(itm[0])

        if len(hits) > 0:
            nh += 1

    for n in xrange(N):
        caches[n].end()

if __name__ == '__main__':
    main(sys.argv[1:])
