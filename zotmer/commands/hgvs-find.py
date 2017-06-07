"""
zot hgvs-find - search for known variants in read data.

Usage:
    zot hgvs-find -X [options] <index> [<variant>...]
    zot hgvs-find [options] <index> <input>...

Options:
    -X              index HGVS variants
    -t              test the index for aliasing
    -f FILENAME     read variants from a file
    -k K            value of k to use [default: 25]
    -g PATH         directory of FASTQ reference sequences
    -w WIDTH        context width [default: 1000]
"""

import sys

import docopt
import yaml

from pykmer.basics import fasta, ham, kmer, kmersList, render
from pykmer.file import openFile, readFasta, readFastq
from zotmer.library.hgvs import parseHGVS, refSeq2Hg19
from zotmer.library.kmers import kmers
from zotmer.library.files import readKmersAndCounts

def pairs(xs):
    assert len(xs) & 1 == 0
    i = 0
    while i + 1 < len(xs):
        yield (xs[i], xs[i + 1])
        i += 2

def both(xs, ys):
    while True:
        try:
            x = xs.next()
            y = ys.next()
            yield(x, y)
        except StopIteration:
            return

def neigh(K, x, d):
    if d == 0:
        return

    for j in xrange(K):
        for y in [1,2,3]:
            z = x ^ (y << (2*j))
            if d == 1:
                yield z
            else:
                for w in neigh(K, z, d - 1):
                    if ham(x, w) == d:
                        yield w


def ball(K, xs, d):
    ys = set(xs)
    i = 0
    while i < d:
        i += 1
        for x in xs:
            ys |= set(neigh(K, x, i))
    return ys

def main(argv):
    opts = docopt.docopt(__doc__, argv)


    K = int(opts['-k'])
    W = int(opts['-w'])

    if opts['-X']:
        variants = opts['<variant>']
        if opts['-f']:
            with openFile(opts['-f']) as f:
                variants += f.read().split()
        vx = {}
        for v in variants:
            x = parseHGVS(v)
            if x is None:
                print >> sys.stderr, "unable to parse %s" % (v,)
                continue
            if x['type'] != 'substitution':
                print >> sys.stderr, "only substitutions are supported at this time: %s" % (v,)
                continue
            x['hgvs'] = v
            acc = x['accession']
            if acc not in vx:
                vx[acc] = []
            vx[acc].append(x)

        d = "."
        if opts['-g']:
            d = opts['-g']

        rs = []
        for (acc, vs) in vx.iteritems():
            if acc not in refSeq2Hg19:
                print >> sys.stderr, "accession %s not supported" % (acc,)
                continue
            h = refSeq2Hg19[acc]
            with openFile(d + "/" + h + ".fa.gz") as f:
                for (nm,seq) in readFasta(f):
                    for v in vs:
                        p = v['position'] - 1
                        wt = seq[p-K:p + K]
                        mut = wt[:K] + v['variant'] + wt[K+1:]
                        #print "wt:\t" + wt
                        #print "mut:\t" + mut
                        wtXs = set(kmersList(K, wt, True))
                        mutXs = set(kmersList(K, mut, True))
                        wtYs = list(wtXs - mutXs)
                        wtYs.sort()
                        mutYs = list(mutXs - wtXs)
                        mutYs.sort()

                        ctx = seq[p - W:p + W]
                        ctxXs = set(kmersList(K, ctx, True))
                        ctxYs = ctxXs - (set(wtYs) | set(mutYs))
                        ctxYs = list(ctxYs)
                        ctxYs.sort()

                        r = {}
                        r['hgvs'] = v['hgvs']
                        r['pos'] = [render(K, y) for y in mutYs]
                        r['neg'] = [render(K, y) for y in wtYs]
                        r['ctx'] = [render(K, y) for y in ctxYs]
                        rs.append(r)
        with open(opts['<index>'], 'w') as f:
            yaml.safe_dump(rs, f)
        return

    with open(opts['<index>']) as f:
        itms = yaml.load(f)

    posIdx = {}
    posRes = {}
    negIdx = {}
    negRes = {}
    ctxIdx = {}
    ctxRes = {}
    hs = []
    for itm in itms:
        h = itm['hgvs']
        hs.append(h)
        posRes[h] = {}
        negRes[h] = {}
        ctxRes[h] = {}
        for x in itm['pos']:
            y = kmer(x)
            if y not in posIdx:
                posIdx[y] = []
            posIdx[y].append(h)
            posRes[h][y] = 0
        for x in itm['neg']:
            y = kmer(x)
            if y not in negIdx:
                negIdx[y] = []
            negIdx[y].append(h)
            negRes[h][y] = 0
        for x in itm['ctx']:
            y = kmer(x)
            if y not in ctxIdx:
                ctxIdx[y] = []
            ctxIdx[y].append(h)
            ctxRes[h][y] = 0
    hs.sort()

    if opts['-t']:
        for fn in opts['<input>']:
            print >> sys.stderr, fn
            with kmers(fn, 'r') as z:
                xs = readKmersAndCounts(z)
                for (x, c) in xs:
                    if x in posIdx:
                        for h in posIdx[x]:
                            posRes[h][x] += c
                    if x in negIdx:
                        for h in negIdx[x]:
                            negRes[h][x] += c
                    if x in ctxIdx:
                        for h in ctxIdx[x]:
                            ctxRes[h][x] += c
        for (h, xs) in posRes.items():
            for (x, c) in xs.items():
                print '%s\tpos\t%s\t%d' % (h, render(K, x), c)
        for (h, xs) in negRes.items():
            for (x, c) in xs.items():
                print '%s\tneg\t%s\t%d' % (h, render(K, x), c)
        for (h, xs) in ctxRes.items():
            for (x, c) in xs.items():
                print '%s\tctx\t%s\t%d' % (h, render(K, x), c)
        return

    for (fn1, fn2) in pairs(opts['<input>']):
        with openFile(fn1) as f1, openFile(fn2) as f2:
            for fq1, fq2 in both(readFastq(f1), readFastq(f2)):
                xs = kmersList(K, fq1[1]) + kmersList(K, fq2[1])
                hits = set([])
                for x in xs:
                    if x in posIdx:
                        for h in posIdx[x]:
                            posRes[h][x] += 1
                            hits.add(h)
                    if x in negIdx:
                        for h in negIdx[x]:
                            negRes[h][x] += 1
                            hits.add(h)

    for (h, xs) in posRes.items():
        for (x, c) in xs.items():
            #if c == 0:
            #    continue
            print '%s\tpos\t%s\t%d' % (h, render(K, x), c)
    for (h, xs) in negRes.items():
        for (x, c) in xs.items():
            #if c == 0:
            #    continue
            print '%s\tneg\t%s\t%d' % (h, render(K, x), c)

if __name__ == '__main__':
    main(sys.argv[1:])
