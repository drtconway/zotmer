"""
Usage:
    zot ksnp <input>...
"""

from pykmer.basics import render
from pykmer.sparse import sparse
import pykmer.kset as kset
from merge import merge1

import array
import docopt
import sys

def groups(K, xs):
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

    for grp in groups(K, u):
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
