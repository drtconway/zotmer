from base import Cmd
from merge import merge1, merge2

from pykmer.basics import kmers
from pykmer.file import readFasta, readFastq
import pykmer.kset as kset
import pykmer.kfset as kfset

import gzip
import os
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

def mkPairs(xs):
    p = 0
    n = 0
    for x in xs:
        if x != p:
            if n > 0:
                yield (p, n)
            p = x
            n = 0
        n += 1
    if n > 0:
        yield (p, n)

def mkSet(xs):
    p = 0
    n = 0
    for x in xs:
        if x != p:
            if n > 0:
                yield p
            p = x
            n = 0
        n += 1
    if n > 0:
        yield p

class Kmerize(Cmd):
    def run(self, opts):
        K = int(opts['<k>'])
        out = opts['<out>']
        s = opts['-s']
        Z = 1024*1024*32
        if opts['-m'] is not None:
            Z = 1024*1024*int(opts['-m'])
        buf = []
        n = 0
        tmps = []
        acgt = [0, 0, 0, 0]
        m = 0
        for fn in opts['<input>']:
            for (nm, seq) in mkParser(fn):
                for x in kmers(K, seq, True):
                    buf.append(x)
                    acgt[x&3] += 1
                    m += 1
                    n += 1
                if n >= Z:
                    buf.sort()
                    fn = 'tmps-%d.k%s%d' % (len(tmps), ('' if s else 'f'), K)
                    tmps.append(fn)
                    if s:
                        kset.write(K, mkSet(buf), fn)
                    else:
                        kfset.write(K, mkPairs(buf), fn)
                    buf = []
                    n = 0

        if len(tmps):
            if len(buf):
                buf.sort()
                fn = 'tmps-%d.k%s%d' % (len(tmps), ('' if s else 'f'), K)
                tmps.append(fn)
                if s:
                    kset.write(K, mkSet(buf), fn)
                else:
                    kfset.write(K, mkPairs(buf), fn)
                buf = []

            zs = None
            if s:
                for fn in tmps:
                    (_, xs) = kset.read(fn)
                    if zs is None:
                        zs = xs
                    else:
                        zs = merge1(K, zs, xs)
            else:
                for fn in tmps:
                    (_, xs) = kfset.read(fn)
                    if zs is None:
                        zs = xs
                    else:
                        zs = merge2(K, zs, xs)
        else:
            buf.sort()
            if s:
                zs = mkSet(buf)
            else:
                zs = mkPairs(buf)

        n = float(sum(acgt))
        acgt = tuple([c/n for c in acgt])

        meta = {
            'total' : m,
            'distinct' : len(buf),
            'acgt' : acgt
        }

        if s:
            kset.write(K, zs, out, meta)
        else:
            kfset.write(K, zs, out, meta)

        for fn in tmps:
            os.remove(fn)

def add(cmds):
    cmds['kmerize'] = Kmerize()
