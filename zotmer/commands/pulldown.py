"""
Usage:
    zot pulldown [-p] <baits> <output> <input>...

Pull down reads (or pairs) that have homology to bait sequences.
The output is the name of a zipfile in to which the reads for each
of the baits will be placed.

Options:
    -p      treat reads as paired end. Requires an even number of inputs.
"""

import os
import sys
import zipfile

import docopt

from pykmer.basics import fasta, ham, kmersList, render
from pykmer.bits import popcnt
from pykmer.file import openFile, readFasta, readFastq, tmpfile
from pykmer.index import buildIndex, index
from pykmer.misc import unionfind, uniq
from pykmer.sparse import sparse

from zotmer.library.kmers import kmers
from zotmer.library.files import readKmers, writeKmers

def masks(K, tau):
    M = (1 << K) - 1
    ms = [M for i in xrange(tau+1)]
    for i in xrange(K):
        for j in xrange(tau+1):
            if (i % (tau+1)) == j:
                ms[j] ^= 1 << i
    return ms

def project(K, msk, x):
    y = 0
    i = 0
    j = 0
    while i < K:
        if msk & 1:
            y |= (x & 3) << j
            j += 2
        x >>= 2
        msk >>= 1
        i +=1
    return y

def cluster(K, xs, tau):
    msks = masks(K, tau)
    J = max([popcnt(msk) for msk in msks])
    uf = unionfind()
    idx = {}
    for x in xs:
        for y in [project(K, msk, x) for msk in msks]:
            if y not in idx:
                idx[y] = x
            else:
                uf.union(x, idx[y])

    idx = {}
    for x in xs:
        a = uf.find(x)
        if a not in idx:
            idx[a] = []
        idx[a].append(x)

    uf = unionfind()
    for ys in idx.values():
        for i in xrange(len(ys)):
            for j in xrange(i+1, len(ys)):
                if ham(ys[i], ys[j]) <= tau:
                    uf.union(ys[i], ys[j])
    
    idx = {}
    for x in xs:
        a = uf.find(x)
        if a not in idx:
            idx[a] = []
        idx[a].append(x)

    return (uf, idx)

def summarize(bxs):
    n = 0
    s = 0.0
    s2 = 0.0
    for acgt in bxs:
        n += 1
        v = sum(acgt)
        s += v
        s2 += v*v
    av = float(s)/float(n)
    vr = float(s2)/float(n) - av*av
    return (n, av, vr)

def fRender(bxs):
    r = []
    for acgt in bxs:
        b = 0
        for i in xrange(4):
            if acgt[i] > 0:
                b |= 1 << i
        r.append(fasta(b))
    return ''.join(r)

def thread(g, K, xs, uf, idx, cxs):
    M = (1 << (2*K)) - 1

    for (a, bs) in idx.iteritems():
        if len(bs) < 2:
            continue
        es = []
        for i in xrange(len(bs)):
            for j in xrange(i + 1, len(bs)):
                d = ham(bs[i], bs[j])
                exs = [(cxs[bs[i]], bs[i]), (cxs[bs[j]], bs[j])]
                exs.sort()
                es.append((d, -exs[0][0], -exs[1][0], exs[0][1], exs[1][1]))
        es.sort()
        seen = set([])
        wf = unionfind()
        print 'graph G {'
        for (d, c0, c1, a0, b0) in es:
            a = wf.find(a0)
            b = wf.find(b0)
            if a != b:
                ax = render(K, a0)
                bx = render(K, b0)
                if a0 not in seen:
                    print '\t%s [label="%s\n%d"];' % (ax, ax, -c0)
                    seen.add(a0)
                if b0 not in seen:
                    print '\t%s [label="%s\n%d"];' % (bx, bx, -c1)
                    seen.add(b0)
                print '\t%s -- %s [label="%d"];' % (ax, bx, d)
                wf.union(a, b)
        print '}'

    return

    eFwd = {}
    eRev = {}
    for i in xrange(xs.count()):
        x = xs.select(i)
        a = uf.find(x)
        eFwd[a] = set([])
        eRev[a] = set([])

    for i in xrange(xs.count()):
        x = xs.select(i)
        a = uf.find(x)
        y0 = (x << 2) & M
        y1 = y0 + 4
        (r0, r1) = xs.rank2(y0, y1)
        for j in xrange(r0, r1):
            y = xs.select(j)
            b = uf.find(y)
            eFwd[a].add(b)
            eRev[b].add(a)

    start = set([])
    for a in eFwd.keys():
        if len(eRev[a]) == 0:
            start.add(a)

    seq = {}
    mgd = {}
    for (a,ys) in idx.iteritems():
        seq[a] = [a]
        v = [[0, 0, 0, 0] for j in xrange(K)]
        for y in ys:
            c = cxs[y]
            for j in xrange(K):
                v[K - 1 - j][y & 3] += c
                y >>= 2
        mgd[a] = v

    vf = unionfind()
    for a in eFwd.keys():
        if len(eFwd[a]) != 1:
            continue
        b = list(eFwd[a])[0]
        if len(eRev[b]) != 1:
            continue
        p, q = vf.find(a), vf.find(b)
        s = seq[p] + seq[q]
        vf.union(p, q)
        r = vf.find(p)
        seq[r] = s
        if r == p:
            del seq[q]
        else:
            del seq[p]

#    seen = set([])
#    print 'digraph G {'
#
#    for (a0, b0s) in eFwd.iteritems():
#        a = vf.find(a0)
#        if a in seen:
#            continue
#        seen.add(a)
#        aa = seq[a]
#        s = mgd[aa[0]]
#        for i in xrange(1, len(aa)):
#            s.append(mgd[aa[i]][-1])
#        (n, av, vr) = summarize(s)
#        ss = fRender(s)
#        print '\tnode X%d [label="%d/%2.2f/%2.2f\n%s"];' % (a, n, av, vr, ss)
#
#    for (a0, b0s) in eFwd.iteritems():
#        a = vf.find(a0)
#        for b0 in b0s:
#            b = vf.find(b0)
#            k = (a,b)
#            if a != b and k not in seen:
#                print '\tX%d -> X%d;' % k
#                seen.add(k)
#    print '}'


    n = 0
    for (a, ys) in seq.items():
        s = [[0, 0, 0, 0] for i in xrange(len(ys) + K - 1)]
        for j in xrange(K):
            for l in xrange(4):
                s[j][l] += mgd[ys[0]][j][l]
        for i in xrange(1, len(ys)):
            for l in xrange(4):
                s[i+K-1][l] += mgd[ys[i]][K-1][l]
        for i in xrange(len(s)):
            print '%d\t%d\t%d\t%s' % (g, n, i, '\t'.join(map(str, s[i])))
        n += 1

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

def readKmers(K, paired, filenames):
    if paired:
        n = 0
        for (fn1, fn2) in pairs(filenames):
            with openFile(fn1) as f1, openFile(fn2) as f2:
                for fq1, fq2 in both(readFastq(f1), readFastq(f2)):
                    yield (kmersList(K, fq1[1]) + kmersList(K, fq2[2]), fq1, fq2)
                    n += 1
                    if (n & 65535) == 0:
                        print >> sys.stderr, 'reads processed:', n
    else:
        n = 0
        for fn in filenames:
            with openFile(fn) as f:
                for fq in readFastq(f):
                    yield kmersList(K, fq[1])
                    n += 1
                    if (n & 65535) == 0:
                        print >> sys.stderr, 'reads processed:', n

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = 27

    nms = []
    idx = {}
    for (nm, seq) in readFasta(openFile(opts['<baits>'])):
        n = len(nms)
        nms.append(nm)
        for x in kmersList(K, seq, True):
            if x not in idx:
                idx[x] = set([])
            idx[x].add(n)

    for x in idx.keys():
        idx[x] = list(idx[x])
        idx[x].sort()

    rn = 0
    if opts['-p']:

        hist = {}
        for (fn1, fn2) in pairs(opts['<input>']):
            tmps = [(tmpfile('_1.fastq'), tmpfile('_2.fastq')) for i in xrange(len(nms))]
            cache = [[[], []] for i in xrange(len(nms))]
            counts = [0 for i in xrange(len(nms))]
            with openFile(fn1) as f1, openFile(fn2) as f2:
                for fq1, fq2 in both(readFastq(f1), readFastq(f2)):
                    hits = set([])
                    for x in kmersList(K, fq1[1]):
                        for i in idx.get(x, []):
                            hits.add(i)
                    for x in kmersList(K, fq2[1]):
                        for i in idx.get(x, []):
                            hits.add(i)
                    n = len(hits)
                    hist[n] = 1 + hist.get(n, 0)
                    for i in hits:
                        counts[i] += 1
                        cache[i][0].append(fq1)
                        cache[i][1].append(fq2)
                        if len(cache[i][0]) >= 1024:
                            with open(tmps[i][0], 'a') as f:
                                for rd in cache[i][0]:
                                    print >> f, rd[0]
                                    print >> f, rd[1]
                                    print >> f, rd[2]
                                    print >> f, rd[3]
                            with open(tmps[i][1], 'a') as f:
                                for rd in cache[i][1]:
                                    print >> f, rd[0]
                                    print >> f, rd[1]
                                    print >> f, rd[2]
                                    print >> f, rd[3]
                            cache[i][0] = []
                            cache[i][1] = []
            for i in xrange(len(cache)):
                if len(cache[i][0]) > 0:
                    with open(tmps[i][0], 'a') as f:
                        for rd in cache[i][0]:
                            print >> f, rd[0]
                            print >> f, rd[1]
                            print >> f, rd[2]
                            print >> f, rd[3]
                    with open(tmps[i][1], 'a') as f:
                        for rd in cache[i][1]:
                            print >> f, rd[0]
                            print >> f, rd[1]
                            print >> f, rd[2]
                            print >> f, rd[3]
                    cache[i][0] = []
                    cache[i][1] = []
            with zipfile.ZipFile(opts['<output>'], 'w', zipfile.ZIP_DEFLATED) as z:
                for i in xrange(len(nms)):
                    if counts[i] > 0:
                        pth = '/'.join(nms[i].split())
                        z.write(tmps[i][0], pth + '/' + fn1)
                        os.remove(tmps[i][0])
                        z.write(tmps[i][1], pth + '/' + fn2)
                        os.remove(tmps[i][1])
        hist = hist.items()
        hist.sort()
        for (n,f) in hist:
            print '%d\t%d' % (n, f)
    else:
        raise "not implemented"

if __name__ == '__main__':

    main(sys.argv[1:])
