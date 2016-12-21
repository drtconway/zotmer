"""
Usage:
    zot jaccard [-A] <input>...

Options:
    -A          print all pairwise distances
"""

import pykmer.exset as exset

import docopt

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

def approxBeta(a, b, x):
    print (a + b + 1) * (1 - x)
    w1 = math.pow(b*x, 1.0/3.0)
    w2 = math.pow(a*(1.0 - x), 1.0/3.0)
    z = 3.0*(w1* (1.0 - 1.0/(9.0*b)) - w2*(1.0 - 1.0/(9.0*a))) / math.sqrt(w1*w1/b + w2*w2/a)
    sqrt2 = math.sqrt(2)
    r1 = math.erf(z/sqrt2)
    r2 = math.erfc(z/sqrt2)
    return (r1, r2)

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    fns = opts['<input>']
    z = 1
    if opts['--all'] is not None:
        z = len(fns)
    for i in xrange(z):
        (xm, xs) = exset.read(fns[i])
        xK = xm['K']
        for j in xrange(i + 1, len(fns)):
            (ym, ys) = exset.read(fns[j])
            yK = ym['K']
            if xK != yK:
                print >> sys.stderr, 'mismatched K:', fns[j]
                return
            (isec, union, d) = jaccard(xs, ys)
            (r1, r2) = approxBeta(isec+1, (union - isec) + 1, 0.99)
            print '%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f' % (fns[i], fns[j], len(xs), len(ys), isec, union, d, r1, r2)
