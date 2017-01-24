"""
Usage:
    zot jaccard [-abp P] <input>...

Compute Jaccard indexes between k-mer sets. By default, indexes are
computed only between the first k-mer set and all the remaining
k-mer sets. If the -a option is given, all pairwise indexes are
computed.  If the -p P option is given, a Null hypothesis test is
performed for the hypothesis that the underlying Jaccard Index is
less than P. This is particularly useful if subsets of k-mers are
being used (NB, if the k-mer sets are large, the statistics can be
very expensive to compute).

Options:
    -a          print all pairwise distances
    -p P        Jaccard distance thresshhold for p-value computation
"""

import array
import math
import sys

import docopt

from pykmer.container import container
from pykmer.container.std import readKmers
from pykmer.stats import logAdd, logChoose

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

def quantBeta(q, m, n):
    lq = math.log(q)
    l = 1e-10
    h = 1 - 1e-10
    while (h - l) > 1e-7:
        x = (h + l) / 2.0
        lp = logIx(x, m, n)
        if lp < lq:
            l = x
        else:
            h = x
    return l

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    fns = opts['<input>']

    Z = 1
    if opts['-a']:
        Z = len(fns)

    p = None
    if opts['-p'] is not None:
        p = float(opts['-p'])

    for i in xrange(Z):
        with container(fns[i], 'r') as z0:
            xK = z0.meta['K']
            xs = array.array('L', readKmers(z0))
            for j in xrange(i + 1, len(fns)):
                with container(fns[j], 'r') as z1:
                    yK = z1.meta['K']
                    ys = array.array('L', readKmers(z1))
                    if xK != yK:
                        print >> sys.stderr, 'mismatched K:', fns[j]
                        sys.exit(1)
                    (isec, union, d) = jaccard(xs, ys)
                    if p is None:
                        print '%s\t%s\t%d\t%d\t%d\t%d\t%f' % (fns[i], fns[j], len(xs), len(ys), isec, union, d)
                    else:
                        pv = logIx(p, isec+1, (union - isec) + 1) / math.log(10)
                        q05 = quantBeta(0.05, isec+1, (union - isec) + 1)
                        q95 = quantBeta(0.95, isec+1, (union - isec) + 1)
                        print '%s\t%s\t%d\t%d\t%d\t%d\t%f\t-%f\t+%f\t%f' % (fns[i], fns[j], len(xs), len(ys), isec, union, d, d - q05, q95 - d, pv)
                    sys.stdout.flush()

if __name__ == '__main__':
    main(sys.argv[1:])
