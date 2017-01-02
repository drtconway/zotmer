"""
Usage:
    zot jaccard [-ap P] <input>...

Options:
    -a          print all pairwise distances
    -p P        Jaccard distance thresshhold for p-value computation.
                Defaults to 0.95
"""

from pykmer.adaptors import kf2k
from pykmer.container import probe
import pykmer.kset as kset
import pykmer.kfset as kfset
from pykmer.stats import logAdd, logChoose

import array
import docopt
import math
import sys

def getKmers(fn):
    (m, _) = probe(fn)
    K = m['K']
    if m['type'] == 'k-mer set':
        (_, xs) = kset.read(fn)
        return (K, array.array('L', xs))
    else:
        (_, xs) = kfset.read(fn)
        return (K, array.array('L', kf2k(xs)))

def jaccard(xs, ys):
    xz = len(xs)
    yz = len(ys)
    i = 0
    j = 0
    d = 0
    b = 0
    while i < xz and j < yz:
        x = xs[i]
        y = ys[j]
        if x < y:
            d += 1
            i += 1
            continue
        if x > y:
            d += 1
            j += 1
            continue
        b += 1
        i += 1
        j += 1
    d += xz - i
    d += yz - j
    return (b, b+d, float(b) / float(b + d))

def logIx(x, m, n):
    lx = math.log(x)
    j = m
    v = logChoose(n + j - 1, j)
    s = v + j * lx
    while True:
        j += 1
        v += math.log((n + j - 1.0) / j)
        t = v + j * lx
        #print j, t, s
        u = logAdd(s, t)
        if u == s:
            break
        s = u
    return n*math.log1p(-x) + s

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    fns = opts['<input>']

    z = 1
    if opts['-a'] is not None:
        z = len(fns)

    p = 0.95
    if opts['-p'] is not None:
        p = float(opts['-p'])

    for i in xrange(z):
        (xK, xs) = getKmers(fns[i])
        for j in xrange(i + 1, len(fns)):
            (yK, ys) = getKmers(fns[j])
            if xK != yK:
                print >> sys.stderr, 'mismatched K:', fns[j]
                sys.exit(1)
            (isec, union, d) = jaccard(xs, ys)
            pv = logIx(p, isec+1, (union - isec) + 1) / math.log(10)
            print '%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f' % (fns[i], fns[j], len(xs), len(ys), isec, union, d, pv)

if __name__ == '__main__':
    main(sys.argv[1:])
