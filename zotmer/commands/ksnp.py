"""
Usage:
    zot ksnp [(-H DIST|-L DIST)] <output> <ref>

Use the kSNP algorithm (and extensions) to find candidate SNPs in a set
of reference k-mers.

Options:
    -H MAX          Perform Hamming based matching allowing a maximum
                    of MAX SNPs in the right hand side of k-mers.
    -L MAX          Perform Levenshtein based matching allowing a
                    maximum of MAX SNPs in the right hand side of
                    k-mers.
"""

import array
import sys

import docopt

from pykmer.basics import ham, lev
from pykmer.bits import popcnt
from pykmer.container import container
from pykmer.misc import unionfind
from pykmer.sparse import sparse

def ksnp(K, xs):
    J = (K - 1) // 2
    S = 2*(J + 1)
    M = (1 << (2*J)) - 1

    itms =  {}
    pfx = 0
    for x in xs:
        y = x >> S
        z = x & M
        if y != pfx:
            for ys in itms.values():
                if len(ys) > 1:
                    yield ys
            itms = {}
            pfx = y
        if z not in itms:
            itms[z] = []
        itms[z].append(x)
    for ys in itms.values():
        if len(ys) > 1:
            yield ys

def hamming(K, d, xs):
    J = (K - 1) // 2
    S = 2*(J + 1)
    T = 2*J
    M = (1 << (2*J)) - 1

    def grp(itms):
        z = len(itms)
        u = unionfind()
        for i in xrange(z):
            for j in xrange(i + 1, z):
                d0 = ham(itms[i], itms[j])
                if d0 <= d:
                    u.union(i, j)
        idx = {}
        for i in xrange(z):
            j = u.find(i)
            if j not in idx:
                idx[j] = []
            idx[j].append(i)
        for ys in idx.itervalues():
            if len(ys) == 1:
                continue
            zs = [itms[y] for y in ys]
            m = 0
            for z in zs:
                m |= 1 << ((z >> T) & 3)
            if popcnt(m) == 1:
                continue
            yield zs

    itms =  []
    pfx = 0
    for x in xs:
        y = x >> S
        if y != pfx:
            for g in grp(itms):
                yield g
            itms = []
            pfx = y
        itms.append(x)
    for g in grp(itms):
        yield g

def levenshtein(K, d, xs):
    J = (K - 1) // 2
    S = 2*(J + 1)
    T = 2*J
    M = (1 << (2*J)) - 1

    def grp(itms):
        z = len(itms)
        u = unionfind()
        for i in xrange(z):
            for j in xrange(i + 1, z):
                d0 = lev(K, itms[i], itms[j])
                if d0 <= d:
                    u.union(i, j)
        idx = {}
        for i in xrange(z):
            j = u.find(i)
            if j not in idx:
                idx[j] = []
            idx[j].append(i)
        for ys in idx.itervalues():
            if len(ys) == 1:
                continue
            zs = [itms[y] for y in ys]
            m = 0
            for z in zs:
                m |= 1 << ((z >> T) & 3)
            if popcnt(m) == 1:
                continue
            yield zs

    itms =  []
    pfx = 0
    for x in xs:
        y = x >> S
        if y != pfx:
            for g in grp(itms):
                yield g
            itms = []
            pfx = y
        itms.append(x)
    for g in grp(itms):
        yield g

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    with container(opts['<ref>'], 'r') as z:
        K = z.meta['K']
        xs = readKmers(z)

        if opts['-H'] is not None:
            d = int(opts['-H'])
            ref = hamming(K, d, xs)
        elif opts['-L'] is not None:
            d = int(opts['-L'])
            ref = levenshtein(K, d, xs)
        else:
            ref = ksnp(K, xs)

        xs = []
        for ys in ref:
            xs += ys
    xs.sort()
    with container(opts['<output>'], 'w') as z:
        writeKmers(K, xs, z)

if __name__ == '__main__':
    main(sys.argv[1:])
