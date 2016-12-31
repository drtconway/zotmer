"""
Usage:
    zot merge [(-s|-f)] <output> <input>...

Options:
    -s          force output to be a k-mer set
    -f          force output to be a k-mer frequency set
"""

from pykmer.adaptors import k2kf, kf2k
from pykmer.container import probe
import pykmer.kset as kset
import pykmer.kfset as kfset

import docopt
import sys

def merge(k, xs, ys, access, combine):
    moreXs = True
    moreYs = True

    try:
        x = xs.next()
        xk = access(x)
    except StopIteration:
        moreXs = False
    try:
        y = ys.next()
        yk = access(y)
    except StopIteration:
        moreYs = False

    while moreXs and moreYs:
        if xk < yk:
            yield x
            try:
                x = xs.next()
                xk = access(x)
            except StopIteration:
                moreXs = False
            continue

        if xk > yk:
            yield y
            try:
                y = ys.next()
                yk = access(y)
            except StopIteration:
                moreYs = False
            continue

        assert xk == yk
        yield combine(x, y)
        try:
            x = xs.next()
            xk = access(x)
        except StopIteration:
            moreXs = False
        try:
            y = ys.next()
            yk = access(y)
        except StopIteration:
            moreYs = False

    while moreXs:
        yield x
        try:
            x = xs.next()
        except StopIteration:
            moreXs = False

    while moreYs:
        yield y
        try:
            y = ys.next()
        except StopIteration:
            moreYs = False

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = None
    types = []
    itrs = []
    ksets = 0
    kfsets = 0
    for inp in opts['<input>']:
        (m, _) = probe(inp)

        K0 = m['K']
        if K is None:
            K = K0
        elif K0 != K:
                raise MismatchedK(K, K0)

        if m['type'] == 'k-mer set':
            (_, xs) = kfset.read(inp)
            types.append(False)
            itrs.append(xs)
            ksets += 1
        else:
            (_, xs) = kfset.read(inp)
            types.append(True)
            itrs.append(xs)
            kfsets += 1

    outFmt = None
    if opts['-s'] is not None and opts['-s']:
        outFmt = False
    elif opts['-f'] is not None and opts['-f']:
        outFmt = True
    elif ksets > 0:
        outFmt = False
    else:
        outFmt = True

    out = opts['<output>']

    if outFmt == False:
        zs = None
        for i in xrange(len(types)):
            xs = itrs[i]
            if types[i] != outFmt:
                xs = kf2k(xs)
            if zs is None:
                zs = xs
            else:
                zs = merge(K, zs, xs, lambda x: x, lambda x, y: x)
        kset.write(K, zs, out)
    else:
        zs = None
        for i in xrange(len(types)):
            xs = itrs[i]
            if types[i] != outFmt:
                xs = k2kf(xs)
            if zs is None:
                zs = xs
            else:
                zs = merge(K, zs, xs, lambda x: x[0], lambda x, y: (x[0], x[1] + y[1]))
        kfset.write(K, zs, out)

if __name__ == '__main__':
    main(sys.argv[1:])
