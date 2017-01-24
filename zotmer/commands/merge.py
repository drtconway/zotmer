"""
Usage:
    zot merge <output> <input>...
"""

import os
import sys

import docopt

from pykmer.basics import render
from pykmer.container import container
from pykmer.container.std import readKmersAndCounts, writeKmersAndCounts
from pykmer.file import tmpfile

def pairs(xs):
    i = 0
    while i + 1 < len(xs):
        yield (xs[i], xs[i + 1])
        i += 2
    if i < len(xs):
        yield (xs[i], )

def merge(xs, ys):
    moreXs = True
    moreYs = True

    try:
        x = xs.next()
        xk = x[0]
    except StopIteration:
        moreXs = False
    try:
        y = ys.next()
        yk = y[0]
    except StopIteration:
        moreYs = False

    while moreXs and moreYs:
        if xk < yk:
            yield x
            try:
                x = xs.next()
                xk = x[0]
            except StopIteration:
                moreXs = False
            continue

        if xk > yk:
            yield y
            try:
                y = ys.next()
                yk = y[0]
            except StopIteration:
                moreYs = False
            continue

        assert xk == yk
        yield (xk, x[1] + y[1])

        try:
            x = xs.next()
            xk = x[0]
        except StopIteration:
            moreXs = False
        try:
            y = ys.next()
            yk = y[0]
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

def hist(zs, h, acgt):
    for z in zs:
        h[z[1]] = 1 + h.get(z[1], 0)
        acgt[z[0]&3] += 1
        yield z

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = None

    out = opts['<output>']

    px = list(pairs(opts['<input>']))
    if len(px) == 1:
        with container(out, 'w') as z:
            h = {}
            acgt = [0, 0, 0, 0]
            ix = px[0]
            if len(ix) == 1:
                with container(ix[0], 'r') as z0:
                    K = z0.meta['K']
                    xs = readKmersAndCounts(z0)
                    zs = hist(xs, h, acgt)
                    writeKmersAndCounts(K, xs, z)
            else:
                with container(ix[0], 'r') as z0:
                    K = z0.meta['K']
                    xs = readKmersAndCounts(z0)
                    with container(ix[1], 'r') as z1:
                        K1 = z1.meta['K']
                        if K1 != K:
                            print >> sys.stderr, "mismatched K"
                            sys.exit(1)
                        ys = readKmersAndCounts(z1)
                        zs = hist(merge(xs, ys), h, acgt)
                        writeKmersAndCounts(K, zs, z)
            n = float(sum(acgt))
            acgt = [c/n for c in acgt]
            z.meta['hist'] = h
            z.meta['acgt'] = acgt
        return

    tmps = []
    tmpnm = tmpfile('.pmc')
    with container(tmpnm, 'w') as z:
        for ix in px:
            if len(ix) == 1:
                nm = 'tmp-' + str(len(tmps))
                tmps.append(nm)
                with container(ix[0], 'r') as z0:
                    if K is None:
                        K = z0.meta['K']
                    else:
                        K0 = z0.meta['K']
                        if K0 != K:
                            print >> sys.stderr, "mismatched K"
                            sys.exit(1)
                    xs = readKmersAndCounts(z0)
                    writeKmersAndCounts(K, xs, z, nm)
            else:
                nm = 'tmp-' + str(len(tmps))
                tmps.append(nm)
                with container(ix[0], 'r') as z0:
                    if K is None:
                        K = z0.meta['K']
                    else:
                        K0 = z0.meta['K']
                        if K0 != K:
                            print >> sys.stderr, "mismatched K"
                            sys.exit(1)
                    xs = readKmersAndCounts(z0)
                    with container(ix[1], 'r') as z1:
                        K1 = z1.meta['K']
                        if K1 != K:
                            print >> sys.stderr, "mismatched K"
                            sys.exit(1)
                        ys = readKmersAndCounts(z1)
                        writeKmersAndCounts(K, merge(xs, ys), z, nm)

    assert K is not None

    with container(out, 'w') as z:
        h = {}
        acgt = [0, 0, 0, 0]
        with container(tmpnm, 'r') as z0:
            zs = None
            for fn in tmps:
                xs = readKmersAndCounts(z0, fn)
                if zs is None:
                    zs = xs
                else:
                    zs = merge(zs, xs)
            zs = hist(zs, h, acgt)
            writeKmersAndCounts(K, zs, z)
        n = float(sum(acgt))
        acgt = [c/n for c in acgt]
        z.meta['hist'] = h
        z.meta['acgt'] = acgt

    os.remove(tmpnm)

if __name__ == '__main__':
    main(sys.argv[1:])
