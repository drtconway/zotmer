"""
Usage:
    zot spoligo [-lp PROBES] <input>...

Perform spoligotyping on sets of k-mers. By default, use the MTB
spoligotyping probes. If the -p option is used, then the argument
is a file name of a text file continaining an optional name and a
probe on each line. Lines starting with '#' are ignored.

Options:
    -l          use long output format with one line per probe
    -p PROBES   the name of a file containing probe sequences
"""

import array
import sys

import docopt

from pykmer.basics import kmer, render
from pykmer.container import container
from pykmer.container.std import readKmers
from pykmer.misc import uniq
from pykmer.sparse import sparse

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

def probe(s):
    return (len(s), kmer(s))

def findApprox(J, x, K, xs, D):
    assert J <= K

    s = 2*(K - J)
    y0 = x << s
    y1 = (x+1) << s
    (r0, r1) = xs.rank2(y0, y1)
    if r1 - r0 > 0:
        return True

    for d in xrange(1, D+1):
        for y in neigh(J, x, d):
            y0 = y << s
            y1 = (y+1) << s
            (r0, r1) = xs.rank2(y0, y1)
            if r1 - r0 > 0:
                return True
    return False

def findProbe(p, K, xs):
    Kp = p[0]
    x = p[1]
    if Kp <= K:
        return findApprox(Kp, x, K, xs, 2)

    # Kp > K.
    M = (1 << (2*K)) - 1
    y0 = x
    ys = []
    for i in xrange(1 + Kp - K):
        y = y0 & M
        if not findApprox(K, y, K, xs, 2):
            return False
        y0 >>= 2
    return True

probesMTB = [
    probe('ATAGAGGGTCGCCGGTTCTGGATCA'),
    probe('CCTCATAATTGGGCGACAGCTTTTG'),
    probe('CCGTGCTTCCAGTGATCGCCTTCTA'),
    probe('ACGTCATACGCCGACCAATCATCAG'),
    probe('TTTTCTGACCACTTGTGCGGGATTA'),
    probe('CGTCGTCATTTCCGGCTTCAATTTC'),
    probe('GAGGAGAGCGAGTACTCGGGGCTGC'),
    probe('CGTGAAACCGCCCCCAGCCTCGCCG'),
    probe('ACTCGGAATCCCATGTGCTGACAGC'),
    probe('TCGACACCCGCTCTAGTTGACTTCC'),
    probe('GTGAGCAACGGCGGCGGCAACCTGG'),
    probe('ATATCTGCTGCCCGCCCGGGGAGAT'),
    probe('GACCATCATTGCCATTCCCTCTCCC'),
    probe('GGTGTGATGCGGATGGTCGGCTCGG'),
    probe('CTTGAATAACGCGCAGTGAATTTCG'),
    probe('CGAGTTCCCGTCAGCGTCGTAAATC'),
    probe('GCGCCGGCCCGCGCGGATGACTCCG'),
    probe('CATGGACCCGGGCGAGCTGCAGATG'),
    probe('TAACTGGCTTGGCGCTGATCCTGGT'),
    probe('TTGACCTCGCCAGGAGAGAAGATCA'),
    probe('TCGATGTCGATGTCCCAATCGTCGA'),
    probe('ACCGCAGACGGCACGATTGAGACAA'),
    probe('AGCATCGCTGATGCGGTCCAGCTCG'),
    probe('CCGCCTGCTGGGTGAGACGTGCTCG'),
    probe('GATCAGCGACCACCGCACCCTGTCA'),
    probe('CTTCAGCACCACCATCATCCGGCGC'),
    probe('GGATTCGTGATCTCTTCCCGCGGAT'),
    probe('TGCCCCGGCGTTTAGCGATCACAAC'),
    probe('AAATACAGGCTCCACGACACGACCA'),
    probe('GGTTGCCCCGCGCCCTTTTCCAGCC'),
    probe('TCAGACAGGTTCGCGTCGATCAAGT'),
    probe('GACCAAATAGGTATCGGCGTGTTCA'),
    probe('GACATGACGGCGGTGCCGCACTTGA'),
    probe('AAGTCACCTCGCCCACACCGTCGAA'),
    probe('TCCGTACGCTCGAAACGCTTCCAAC'),
    probe('CGAAATCCAGCACCACATCCGCAGC'),
    probe('CGCGAACTCGTCCACAGTCCCCCTT'),
    probe('CGTGGATGGCGGATGCGTTGTGCGC'),
    probe('GACGATGGCCAGTAAATCGGCGTGG'),
    probe('CGCCATCTGTGCCTCATACAGGTCC'),
    probe('GGAGCTTTCCGGCTTCTATCAGGTA'),
    probe('ATGGTGGGACATGGACGAGCGCGAC'),
    probe('CGCAGAATCGCACCGGGTGCGGGAG')
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

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    nms = [str(i+1) for i in xrange(len(probesMTB))]
    probes = probesMTB

    if opts['-p'] is not None:
        nms = []
        probes = []
        bad = False
        with open(opts['-p']) as f:
            i = 0
            ln = 0
            for l in f:
                ln += 1
                if l[0] == '#':
                    continue
                i += 1
                t = l.split()
                if len(t) == 1:
                    nms.append(str(i))
                    probes.append(probe(t[0]))
                elif len(t) == 2:
                    nms.append(t[0])
                    probes.append(probe(t[1]))
                else:
                    bad = True
                    print >> sys.stderr, '%s line %d, badly formatted.' % (opts['-p'], i)
        if bad:
            sys.exit(1)

    for inp in opts['<input>']:
        with container(inp, 'r') as z:
            K = z.meta['K']
            xs = readKmers(z)
            xs = sparse(2*K, array.array('L', xs))

        res = []
        for i in xrange(len(probes)):
            if findProbe(probes[i], K, xs):
                res.append('1')
            else:
                res.append('0')
        if opts['-l']:
            for i in xrange(len(nms)):
                print '%s\t%s\t%s' % (inp, nms[i], res[i])
        else:
            print inp + '\t' + ''.join(res)

if __name__ == '__main__':
    main(sys.argv[1:])
