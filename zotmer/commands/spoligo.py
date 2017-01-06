"""
Usage:
    zot spoligo <input>...

Perform MTB spoligotyping on sets of k-mers.
"""

from pykmer.adaptors import kf2k
from pykmer.basics import kmer, render
from pykmer.misc import uniq
from pykmer.sparse import sparse
from pykmer.container import probe
import pykmer.kfset as kfset
import pykmer.kset as kset

import docopt
import sys

probesMTB = [
    kmer('ATAGAGGGTCGCCGGTTCTGGATCA'),
    kmer('CCTCATAATTGGGCGACAGCTTTTG'),
    kmer('CCGTGCTTCCAGTGATCGCCTTCTA'),
    kmer('ACGTCATACGCCGACCAATCATCAG'),
    kmer('TTTTCTGACCACTTGTGCGGGATTA'),
    kmer('CGTCGTCATTTCCGGCTTCAATTTC'),
    kmer('GAGGAGAGCGAGTACTCGGGGCTGC'),
    kmer('CGTGAAACCGCCCCCAGCCTCGCCG'),
    kmer('ACTCGGAATCCCATGTGCTGACAGC'),
    kmer('TCGACACCCGCTCTAGTTGACTTCC'),
    kmer('GTGAGCAACGGCGGCGGCAACCTGG'),
    kmer('ATATCTGCTGCCCGCCCGGGGAGAT'),
    kmer('GACCATCATTGCCATTCCCTCTCCC'),
    kmer('GGTGTGATGCGGATGGTCGGCTCGG'),
    kmer('CTTGAATAACGCGCAGTGAATTTCG'),
    kmer('CGAGTTCCCGTCAGCGTCGTAAATC'),
    kmer('GCGCCGGCCCGCGCGGATGACTCCG'),
    kmer('CATGGACCCGGGCGAGCTGCAGATG'),
    kmer('TAACTGGCTTGGCGCTGATCCTGGT'),
    kmer('TTGACCTCGCCAGGAGAGAAGATCA'),
    kmer('TCGATGTCGATGTCCCAATCGTCGA'),
    kmer('ACCGCAGACGGCACGATTGAGACAA'),
    kmer('AGCATCGCTGATGCGGTCCAGCTCG'),
    kmer('CCGCCTGCTGGGTGAGACGTGCTCG'),
    kmer('GATCAGCGACCACCGCACCCTGTCA'),
    kmer('CTTCAGCACCACCATCATCCGGCGC'),
    kmer('GGATTCGTGATCTCTTCCCGCGGAT'),
    kmer('TGCCCCGGCGTTTAGCGATCACAAC'),
    kmer('AAATACAGGCTCCACGACACGACCA'),
    kmer('GGTTGCCCCGCGCCCTTTTCCAGCC'),
    kmer('TCAGACAGGTTCGCGTCGATCAAGT'),
    kmer('GACCAAATAGGTATCGGCGTGTTCA'),
    kmer('GACATGACGGCGGTGCCGCACTTGA'),
    kmer('AAGTCACCTCGCCCACACCGTCGAA'),
    kmer('TCCGTACGCTCGAAACGCTTCCAAC'),
    kmer('CGAAATCCAGCACCACATCCGCAGC'),
    kmer('CGCGAACTCGTCCACAGTCCCCCTT'),
    kmer('CGTGGATGGCGGATGCGTTGTGCGC'),
    kmer('GACGATGGCCAGTAAATCGGCGTGG'),
    kmer('CGCCATCTGTGCCTCATACAGGTCC'),
    kmer('GGAGCTTTCCGGCTTCTATCAGGTA'),
    kmer('ATGGTGGGACATGGACGAGCGCGAC'),
    kmer('CGCAGAATCGCACCGGGTGCGGGAG')
]

def diff(xs, ys):
    zX = len(xs)
    zY = len(ys)
    i = 0
    j = 0
    k = 0
    while i < zX and j < zY:
        assert k <= i
        x = xs[i]
        y = ys[j]
        if x < y:
            if k != i:
                xs[k] = xs[i]
            i += 1
            k += 1
            continue
        if x > y:
            j += 1
            continue
        i += 1
        j += 1
    if k != i:
        while i < zX:
            xs[k] = xs[i]
            i += 1
            k += 1
        del xs[k:]

def neigh(K, x, d):
    if d == 0:
        return []
    xs = []
    for i in xrange(K):
        for j in xrange(3):
            xs.append(x ^ ((j + 1) << (2*i)))
    xs.sort()

    if d == 1:
        return xs

    zs = []
    for y in xs:
        zs += neigh(K, y, d - 1)
    zs.sort()
    uniq(zs)
    diff(zs, [x])
    diff(zs, xs)
    return zs

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    for inp in opts['<input>']:
        (m, _) = probe(inp)
        if m['type'] == 'k-mer set':
            (m, xs) = kset.read(inp)
            K = m['K']
        else:
            (m, xs) = kfset.read(inp)
            K = m['K']
            xs = kf2k(xs)

        # all MTB probes are 25bp
        J = 2*(K - 25)
        if J < 0:
            print >> sys.stderr, "K < 25. Cannot reliably spoligotype."
            sys.exit(1)

        idx = {}

        for i in xrange(len(probesMTB)):
            idx[probesMTB[i]] = i
            for y in neigh(25, probesMTB[i], 1):
                if y in idx:
                    print >> sys.stderr, 'ambiguity between probes %d and %d' % (i, idx[y])
                idx[y] = i
            for y in neigh(25, probesMTB[i], 2):
                if y in idx:
                    print >> sys.stderr, 'ambiguity between probes %d and %d' % (i, idx[y])
                idx[y] = i

        counts = [0 for i in xrange(len(probesMTB))]

        for x in xs:
            y = x >> J
            if y in idx:
                counts[idx[y]] += 1

        res = []
        for i in xrange(len(counts)):
            if counts[i] > 0:
                res.append('1')
            else:
                res.append('0')
        print inp + '\t' + ''.join(res)

if __name__ == '__main__':
    main(sys.argv[1:])
