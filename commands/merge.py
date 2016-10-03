from base import Cmd

from pykmer.basics import render
import pykmer.kfset as kfset

def merge(k, xs, ys):
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

class Merge(Cmd):
    def run(self, opts):
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
                zs = merge(K, zs, xs)
        if K is not None:
            kfset.write(K, zs, out)

def add(cmds):
    cmds['merge'] = Merge()

