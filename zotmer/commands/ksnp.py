"""
Usage:
    zot ksnp [(-H DIST|-L DIST)] <input>...

Use the kSNP algorithm to find SNPs that distinguish the k-mer sets
in the argument list.

Options:
    -H DIST         Perform Hamming based matching allowing a maximum
                    of DIST SNPs in the right hand side of k-mers.
    -L DIST         Perform Levenshtein based matching allowing a
                    maximum of DIST SNPs in the right hand side of
                    k-mers.
"""

from pykmer.basics import ham, lev, render
from pykmer.bits import popcnt
from pykmer.misc import unionfind
from pykmer.sparse import sparse
import pykmer.kset as kset
from merge import merge1

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

    K = None
    ks = []
    ss = []
    for inp in opts['<input>']:
        (meta, xs) = kset.read(inp)
        if K is None:
            K = meta['K']
        else:
            K0 = meta['K']
            if K != K0:
                print >> sys.stderr, 'kset %s has unexpected K (%d)' % (inp, K0)
                sys.exit(1)
        ks.append(array.array('L', xs))
        ss.append(sparse(2*K, ks[-1]))

    u = ks[0].__iter__()
    for i in xrange(1, len(ks)):
        u = merge1(K, u, ks[i].__iter__())

    if opts['-H'] is not None:
        d = int(opts['-H'])
        groups = hamming(K, d, u)
    elif opts['-L'] is not None:
        d = int(opts['-L'])
        groups = levenshtein(K, d, u)
    else:
        groups = ksnp(K, u)

    for grp in groups:
        ms = []
        for x in grp:
            vs = long(0)
            for i in xrange(len(ss)):
                if ss[i].access(x):
                    vs |= 1 << i
            ms.append(vs)
        ok = True
        for i in xrange(len(ms)):
            for j in xrange(i + 1, len(ms)):
                if (ms[i] & ms[j]) != 0:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            print '\t'.join([render(K, x) for x in grp])

if __name__ == '__main__':
    main(sys.argv[1:])
