"""
Usage:
    zot kmerize [-m MEM] [-D FRAC [-S SEED]] <k> <output> <input>...

Kmerize FASTA or FASTQ inputs to produce a standard container object.

Arguments:
    <k>         the length of the k-mers. Recommended values: 10-30
    <output>    the name of the output file.
                recommended naming convention
                    - mykmers.k25 for a k-mer set of 25-mers
                    - mykmers.kf25 for a k-mer frequency set of 25-mers
                    - mykmers.e25 for an expanded k-mer set of 25-mers

Options:
    -m MEM      in-memory buffer size (in MB)
    -D FRAC     subsample k-mers, using FRAC proportion of k-mers
    -S SEED     if -D is given, give a seed for determining the
                subspace (defaults to 0).
"""

import array
import gzip
import os
import subprocess
import sys
import time

import docopt

from pykmer.basics import kmersList, render, sub
from pykmer.container import container
from pykmer.container.std import readKmersAndCounts, writeKmersAndCounts
from pykmer.file import openFile, readFasta, readFastqBlock, tmpfile

from zotmer.library.merge import merge2

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

    def addList(self, xs):
        for x in xs:
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

    def stat(self):
        h = {}
        for (xh, xls) in self.toc.iteritems():
            l = len(xls)
            h[l] = 1 + h.get(l, 0)

        h = h.items()
        h.sort()
        for (l,c) in h:
            print '%d\t%d' % (l, c)

def pairs(xs):
    i = 0
    res = []
    while i + 1 < len(xs):
        res.append((xs[i], xs[i + 1]))
        i += 2
    if i < len(xs):
        res.append((xs[i], ))
    return res

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
            yield [(nm, seq)]
    else:
        for grps in readFastqBlock(openFile(fn)):
            yield [(grp[0], grp[1]) for grp in grps]

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

    tmpnm = tmpfile('.pmc')
    with container(tmpnm, 'w') as z:
        pass

    PN = 1024*1024

    nr = 0
    t0 = time.time()
    for fn in opts['<input>']:
        for rds in mkParser(fn):
            for (nm, seq) in rds:
                nr += 1
                if nr & (PN - 1) == 0:
                    t1 = time.time()
                    print >> sys.stderr, 'reads processed:', nr, (PN)/(t1 - t0), 'reads/second'
                    t0 = t1
                    #buf.stat()
                xs = kmersList(K, seq, True)
                if d is None:
                    buf.addList(xs)
                    for x in xs:
                        acgt[x&3] += 1
                        m += 1
                        n += 1
                else:
                    for x in xs:
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
                    #print >> sys.stderr, "writing " + fn + "\t" + tmpnm
                    tmps.append(fn)
                    with container(tmpnm, 'a') as z:
                        writeKmersAndCounts(K, mkPairs(buf.kmers()), z, fn)
                    buf.clear()
                    n = 0

    t1 = time.time()
    print >> sys.stderr, 'reads processed:', nr, (nr % PN)/(t1 - t0), 'reads/second'

    if len(tmps) and len(buf):
        fn = 'tmps-%d' % (len(tmps),)
        #print >> sys.stderr, "writing " + fn + "\t" + tmpnm
        tmps.append(fn)
        with container(tmpnm, 'a') as z:
            writeKmersAndCounts(K, mkPairs(buf.kmers()), z, fn)
        buf = []

    while len(tmps) > 2:
        tmpnm2 = tmpfile('.pmc')
        tmps2 = []
        with container(tmpnm, 'r') as z0, container(tmpnm2, 'w') as z:
            ps = pairs(tmps)
            for p in ps:
                fn = 'tmps-%d' % (len(tmps2),)
                tmps2.append(fn)
                if len(p) == 1:
                    writeKmersAndCounts(K, readKmersAndCounts(z0, p[0]), z, fn)
                    continue
                h = {}
                merge2(z, K, readKmersAndCounts(z0, p[0]), readKmersAndCounts(z0, p[1]), h, fn)
        os.remove(tmpnm)
        tmpnm = tmpnm2
        tmps = tmps2

    with container(out, 'w') as z:
        h = {}
        if len(tmps) == 0:
            zs = hist(mkPairs(buf.kmers()), h)
            writeKmersAndCounts(K, zs, z)
        elif len(tmps) == 1:
            with container(tmpnm, 'r') as z0:
                writeKmersAndCounts(K, hist(readKmersAndCounts(z0, tmps[0]), h), z)
        else:
            assert len(tmps) == 2
            with container(tmpnm, 'r') as z0:
                merge2(z, K, readKmersAndCounts(z0, tmps[0]), readKmersAndCounts(z0, tmps[1]), h)
        n = float(sum(acgt))
        acgt = [c/n for c in acgt]
        z.meta['hist'] = h
        z.meta['acgt'] = acgt
        z.meta['reads'] = nr
    #os.remove(tmpnm)

if __name__ == '__main__':
    main(sys.argv[1:])
