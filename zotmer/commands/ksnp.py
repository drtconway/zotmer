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

from pykmer.basics import ham, lev
from pykmer.bits import popcnt
from pykmer.misc import unionfind
from pykmer.sparse import sparse
import pykmer.kset as kset
import pykmer.kfset as kfset

import array
import docopt
import sys

def ksnp(K, xs):
    J = (K - 1) // 2
    S = 2*(J + 1)
    M = (1 << (2*J)) - 1

    gxs =  {}
    gx = 0
    for x in xs:
        y = x >> S
        z = x & M
        if y != gx:
            for ys in gxs.values():
                if len(ys) > 1:
                    yield ys
            gxs = {}
            gx = y
        if z not in gxs:
            gxs[z] = []
        gxs[z].append(x)
    for ys in gxs.values():
        if len(ys) > 1:
            yield ys

def hamming(K, d, xs):
    J = (K - 1) // 2
    S = 2*(J + 1)
    T = 2*J
    M = (1 << (2*J)) - 1

    def grp(gxs):
        z = len(gxs)
        u = unionfind()
        for i in xrange(z):
            for j in xrange(i + 1, z):
                d0 = ham(gxs[i], gxs[j])
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
            zs = [gxs[y] for y in ys]
            m = 0
            for z in zs:
                m |= 1 << ((z >> T) & 3)
            if popcnt(m) == 1:
                continue
            yield zs

    gxs =  []
    gx = 0
    for x in xs:
        y = x >> S
        if y != gx:
            for g in grp(gxs):
                yield g
            gxs = []
            gx = y
        gxs.append(x)
    for g in grp(gxs):
        yield g

def levenshtein(K, d, xs):
    J = (K - 1) // 2
    S = 2*(J + 1)
    T = 2*J
    M = (1 << (2*J)) - 1

    def grp(gxs):
        z = len(gxs)
        u = unionfind()
        for i in xrange(z):
            for j in xrange(i + 1, z):
                d0 = lev(K, gxs[i], gxs[j])
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
            zs = [gxs[y] for y in ys]
            m = 0
            for z in zs:
                m |= 1 << ((z >> T) & 3)
            if popcnt(m) == 1:
                continue
            yield zs

    gxs =  []
    gx = 0
    for x in xs:
        y = x >> S
        if y != gx:
            for g in grp(gxs):
                yield g
            gxs = []
            gx = y
        gxs.append(x)
    for g in grp(gxs):
        yield g

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    refFn = opts['<ref>']
    (meta, xs) = kset.read(refFn)

    K = meta['K']

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

    kset.write(K, xs, opts['<output>'])

if __name__ == '__main__':
    main(sys.argv[1:])
