"""
Usage:
    zot kmerize [-m MEM] [-D FRAC [-S SEED]] <k> <output> <input>...

Kmerize FASTA or FASTQ inputs to produce either a k-mer set or a k-mer
frequency set. If neither -e nor -s are given, a k-mer frequency set
is generated.

Arguments:
    <k>         the length of the k-mers. Recommended values: 10-30
    <output>    the name of the output file.
                recommended naming convention
                    - mykmers.k25 for a k-mer set of 25-mers
                    - mykmers.kf25 for a k-mer frequency set of 25-mers
                    - mykmers.e25 for an expanded k-mer set of 25-mers

Options:
    -m MEM      in-memory buffer size
    -D FRAC     subsample k-mers, using FRAC proportion of k-mers
    -S SEED     if -D is given, give a seed for determining the
                subspace (defaults to 0).
"""

from pykmer.basics import kmers, render, sub
from pykmer.file import readFasta, readFastq
from pykmer.container import container
from pykmer.container.std import readKmersAndCounts, writeKmersAndCounts

from merge import merge

import docopt
import array
import gzip
import os
import subprocess
import sys
import time

class KmerAccumulator:
    def __init__(self):
        self.toc = {}
        self.z = 0

    def __len__(self):
        return self.z

    def clear(self):
        self.toc = {}
        self.z = 0

    def add(self, x):
        xh = x >> 32
        xl = x & 0xFFFFFFFF
        if xh not in self.toc:
            self.toc[xh] = array.array('I')
        self.toc[xh].append(xl)
        self.z += 1

    def kmers(self):
        xhs = self.toc.keys()
        xhs.sort()
        for xh in xhs:
            x0 = xh << 32
            xls = self.toc[xh].tolist()
            xls.sort()
            for xl in xls:
                x = x0 | xl
                yield x

def openFile(fn):
    if fn == "-":
        return sys.stdin
    if fn.endswith(".gz"):
        p = subprocess.Popen(['gunzip', '-c', fn], stdout=subprocess.PIPE)
        return p.stdout
    if fn.endswith(".bz2"):
        p = subprocess.Popen(['bunzip2', '-c', fn], stdout=subprocess.PIPE)
        return p.stdout
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
    m = 0
    p = 0
    n = 0
    for x in xs:
        if x != p:
            if n > 0:
                m += 1
                yield (p, n)
            p = x
            n = 0
        n += 1
    if n > 0:
        m += 1
        yield (p, n)
    #print >> sys.stderr, 'wrote %d k-mers' % (m,)

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

def hist(zs, h):
    for z in zs:
        h[z[1]] = 1 + h.get(z[1], 0)
        yield z

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = int(opts['<k>'])
    out = opts['<output>']
    Z = 1024*1024*32
    if opts['-m'] is not None:
        Z = 1024*1024*int(opts['-m'])

    buf = KmerAccumulator()
    n = 0
    tmps = []
    acgt = [0, 0, 0, 0]
    m = 0

    d = None
    if opts['-D'] is not None:
        d = float(opts['-D'])

        S = 0
        if opts['-S'] is not None:
            S = int(opts['-S'])

        cacheYes = set([])
        cacheNo = set([])

    with container('tmp-' + out, 'w') as z:
        pass

    nr = 0
    t0 = time.time()
    for fn in opts['<input>']:
        for (nm, seq) in mkParser(fn):
            nr += 1
            if nr & (1024*1024 - 1) == 0:
                t1 = time.time()
                print >> sys.stderr, 'reads processed:', nr, (1024*1024)/(t1 - t0), 'reads/second'
                t0 = t1
            if d is None:
                for x in kmers(K, seq, True):
                    buf.add(x)
                    acgt[x&3] += 1
                    m += 1
                    n += 1
            else:
                for x in kmers(K, seq, True):
                    if x in cacheNo:
                        continue
                    if x not in cacheYes:
                        if not sub(S, d, x):
                            cacheNo.add(x)
                            continue
                        cacheYes.add(x)
                    buf.add(x)
                    acgt[x&3] += 1
                    m += 1
                    n += 1
                if len(cacheYes) > 1000000:
                    cacheYes = set([])
                if len(cacheNo) > 1000000:
                    cacheNo = set([])
            if 8*n >= Z:
                fn = 'tmps-%d' % (len(tmps),)
                tmps.append(fn)
                with container('tmp-' + out, 'a') as z:
                    writeKmersAndCounts(K, mkPairs(buf.kmers()), z, fn)
                buf.clear()
                n = 0

    if len(tmps):
        if len(buf):
            fn = 'tmps-%d' % (len(tmps),)
            tmps.append(fn)
            with container('tmp-' + out, 'a') as z:
                writeKmersAndCounts(K, mkPairs(buf.kmers()), z, fn)
            buf = []

    with container(out, 'w') as z:
        h = {}
        if len(tmps):
            with container('tmp-' + out, 'r') as z0:
                zs = None
                for fn in tmps:
                    xs = readKmersAndCounts(z0, fn)
                    if zs is None:
                        zs = xs
                    else:
                        zs = merge(zs, xs)
                zs = hist(zs, h)
                writeKmersAndCounts(K, zs, z)
        else:
            zs = hist(mkPairs(buf.kmers()), h)
            writeKmersAndCounts(K, zs, z)
        n = float(sum(acgt))
        acgt = [c/n for c in acgt]
        z.meta['hist'] = h
        z.meta['acgt'] = acgt
        z.meta['reads'] = nr
    os.remove('tmp-' + out)

if __name__ == '__main__':
    main(sys.argv[1:])
