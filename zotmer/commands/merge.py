"""
Usage:
    zot merge <output> <input>...
"""

from pykmer.basics import render
import pykmer.kfset as kfset

import docopt

def merge1(k, xs, ys):
    moreXs = True
    moreYs = True

    try:
        x = xs.next()
    except StopIteration:
        moreXs = False
    try:
        y = ys.next()
    except StopIteration:
        moreYs = False

    while moreXs and moreYs:
        if x < y:
            yield x
            try:
                x = xs.next()
            except StopIteration:
                moreXs = False
            continue

        if x > y:
            yield y
            try:
                y = ys.next()
            except StopIteration:
                moreYs = False
            continue

        assert x == y
        yield x
        try:
            x = xs.next()
        except StopIteration:
            moreXs = False
        try:
            y = ys.next()
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

def merge2(k, xs, ys):
    moreXs = True
    moreYs = True

    try:
        (x,xf) = xs.next()
    except StopIteration:
        moreXs = False
    try:
        (y,yf) = ys.next()
    except StopIteration:
        moreYs = False

    while moreXs and moreYs:
        if x < y:
            yield (x, xf)
            try:
                (x,xf) = xs.next()
            except StopIteration:
                moreXs = False
            continue

        if x > y:
            yield (y, yf)
            try:
                (y,yf) = ys.next()
            except StopIteration:
                moreYs = False
            continue

        assert x == y
        yield (x, xf + yf)
        try:
            (x,xf) = xs.next()
        except StopIteration:
            moreXs = False
        try:
            (y,yf) = ys.next()
        except StopIteration:
            moreYs = False

    while moreXs:
        yield (x, xf)
        try:
            (x,xf) = xs.next()
        except StopIteration:
            moreXs = False

    while moreYs:
        yield (y, yf)
        try:
            (y,yf) = ys.next()
        except StopIteration:
            moreYs = False

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = None
    zs = None
    out = opts['<output>']
    for inp in opts['<input>']:
        (m, xs) = kfset.read(inp)
        K0 = m['K']
        if K is None :
            K = K0
            zs = xs
        else:
            if K != K0:
                raise MismatchedK(K, K0)
            zs = merge2(K, zs, xs)
    if K is not None:
        kfset.write(K, zs, out)
