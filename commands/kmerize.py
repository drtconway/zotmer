from base import Cmd

from pykmer.basics import kmers
from pykmer.file import readFasta, readFastq
import pykmer.kset as kset
import pykmer.kfset as kfset

import gzip
import sys

def openFile(fn):
    if fn == "-":
        return sys.stdin
    if fn.endswith(".gz"):
        return gzip.open(fn)
    return open(fn, 'rb')

def stripCompressionSuffix(nm):
    if nm.endswith('.gz'):
        return nm[:-3]
    return nm

def isFasta(nm):
    bnm = stripCompressionSuffix(nm)
    if bnm.endswith(".fa"):
        return True
    if bnm.endswith(".fasta"):
        return True
    if bnm.endswith(".fas"):
        return True
    if bnm.endswith(".fna"):
        return True

def mkParser(fn):
    if isFasta(fn):
        for (nm, seq) in readFasta(openFile(fn)):
            yield (nm, seq)
    else:
        for grp in readFastq(openFile(fn)):
            yield (grp[0], grp[1])

class Kmerize(Cmd):
    def run(self, opts):
        K = int(opts['<k>'])
        out = opts['<out>']
        s = opts['-s']
        Z = 1024*1024*256
        if opts['-m'] is not None:
            Z = 1024*1024*int(opts['-m'])
        buf = {}
        acgt = [0, 0, 0, 0]
        m = 0
        for fn in opts['<input>']:
            for (nm, seq) in mkParser(fn):
                for x in kmers(K, seq, True):
                    c = buf.get(x, 0)
                    buf[x] = c + 1
                    acgt[x&3] += 1
                    m += 1

        n = float(sum(acgt))
        acgt = tuple([c/n for c in acgt])
        h = {}
        for (_, c) in buf.iteritems():
            n = h.get(c, 0)
            h[c] = n + 1

        meta = {
            'total' : m,
            'distinct' : len(buf),
            'acgt' : acgt,
            'hist' : h
        }

        if s:
            buf = buf.keys()
            buf.sort()
            kset.write(K, buf, out, meta)
        else:
            buf = buf.items()
            buf.sort()
            kfset.write(K, buf, out, meta)

def add(cmds):
    cmds['kmerize'] = Kmerize()
