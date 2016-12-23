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
from pykmer.stats import gammaP, logGammaP, log1mexp

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

def contFrac(f, g, n):
    t = 0
    for i in xrange(n):
        j = n - i
        t = g(j) / float(f(j) + t)
    return f(0) + t

def logBeta(a, b, x, n):
    inv = False
    if x >=  (a + 1.0)/(a + b + 1.0):
        x = 1 - x
        t = b
        b = a
        a = t
        inv = True
    def d(i):
        if i == 0:
            return 1.0
        m = i // 2
        if i & 1:
            return - ((a + m)*(a + b + m)*x)/((a + 2*m)*(a + 2*m + 1))
        else:
            return (m * (b - m) * x)/((a + 2*m - 1)*(a + 2*m))

    s = contFrac(d, lambda i: 1.0, n)
    ls = math.log(s)
    r = a*math.log(x) + b*math.log1p(x) + ls
    if inv:
        r = log1mexp(r)
    return r

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
            pv = logBeta(isec+1, (union - isec) + 1, p, 10) / math.log(10)
            print '%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f' % (fns[i], fns[j], len(xs), len(ys), isec, union, d, pv)

if __name__ == '__main__':
    main(sys.argv[1:])
