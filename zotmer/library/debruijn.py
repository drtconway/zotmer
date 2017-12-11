"""
A collection of de Buijn graph manipulation tools.
"""

import pytest

from pykmer.basics import kmersList

def is_in(i, n):
    if isinstance(n, int):
        return i == n
    assert isinstance(n, slice)
    lo = 0
    if n.start is not None:
        lo = n.start
    if i < lo:
        return False
    if n.stop is not None and i >= n.stop:
        return False
    if n.step is not None:
        j = (i - lo) % n.step
        if j != 0:
            return False
    return True

def test_is_in_0():
    assert is_in(0, 0)

def test_is_in_1():
    assert is_in(1, 1)

def test_is_in_2():
    n = slice(10)
    for i in range(10):
        assert is_in(i, n)

def less_than_max(i, n):
    if isinstance(n, int):
        return i <= n
    assert isinstance(n, slice)
    if n.stop is not None and i > n.stop:
        return False
    return True

def test_less_than_max_0():
    n = 0
    i = 0
    assert not less_than_max(i, n)

def test_less_than_max_1():
    n = 1
    i = 0
    assert less_than_max(i, n)

def test_less_than_max_2():
    n = 10
    for i in range(n):
        assert less_than_max(i, n)

def test_less_than_max_3():
    n = slice(10)
    for i in range(10):
        assert less_than_max(i, n)

def follows(s, t):
    if s[1:] != t[:-1]:
        print '%s\t%s' % (s, t)
    assert s[1:] == t[:-1]

def cartesian(xs, ys):
    for x in xs:
        for y in ys:
            yield (x,y)

class Interpolator(object):
    def __init__(self, K, xs):
        self.K = K
        self.M = ((1 << (2*K)) - 1)
        self.S = 2*(K-1)
        self.xs = xs

    def succ(self, x):
        y0 = (x << 2) & self.M
        for i in range(4):
            y = y0 + i
            if y in self.xs:
                yield y

    def pred(self, x):
        y0 = x >> 2
        for i in range(4):
            y = y0 + (i << self.S)
            if y in self.xs:
                yield y

    def path(self, xb, xe, n):
        fwd = {}
        fwd[xb] = [[xb]]
        rev = {}
        rev[xe] = [[xe]]

        i = 1
        while less_than_max(i, n):
            if is_in(i, n):
                for x in set(fwd.keys()) & set(rev.keys()):
                    for (fst, snd) in cartesian(fwd[x], rev[x]):
                        yield fst + snd[1:]
            if (i & 1) == 0:
                nextFwd = {}
                for (x,ps) in fwd.iteritems():
                    for y in self.succ(x):
                        if y not in nextFwd:
                            nextFwd[y] = []
                        for p in ps:
                            nextFwd[y].append(p + [y])
                fwd = nextFwd
            else:
                nextRev = {}
                for (x,ps) in rev.iteritems():
                    for y in self.pred(x):
                        if y not in nextRev:
                            nextRev[y] = []
                        for p in ps:
                            nextRev[y].append([y] + p)
                rev = nextRev
            i += 1

    def interpolate(self, xb, xe, n):
        sMax = 0
        pMax = None
        for p in self.path(xb, xe, n):
            s = 0
            cc = []
            for x in p:
                s += self.xs[x]
                cc.append(self.xs[x])
            if s > sMax:
                sMax = s
                pMax = p
        return pMax

def paths(K, xs, xb, xe, n):
    """
    Find paths through the de bruijn graph implied by xs, beginning
    at xb and ending at xe. n specifies the length of the path. If
    it is an integer, then it gives an exact length of the path,
    otherwise it must be a slice which defines a half open interval
    within which the length of the path must lie.  If such a path
    exists, return the sequence of k-mers that form the path,
    including the endpoints. NB if xb == xe, then the path will be
    length 1.
    """
    I = Interpolator(K, xs)
    for p in I.path(xb, xe, n):
        yield p

def interpolate(K, xs, xb, xe, n):
    """
    Find a path through the de bruijn graph implied by xs, beginning
    at xb and ending at xe. n specifies the length of the path. If
    it is an integer, then it gives an exact length of the path,
    otherwise it must be a slice which defines a half open interval
    within which the length of the path must lie.  If such a path
    exists, return the sequence of k-mers that form the path,
    including the endpoints. NB if xb == xe, then the path will be
    length 1.

    If there are multiple paths, the one with the highest aggregate
    coverage is returned.
    """
    I = Interpolator(K, xs)
    p = I.interpolate(xb, xe, n)
    return p

def test_interpolate0Exacxt() :
    K = 25
    seq = "GGAGTTTCCAAGAGAAAATTTAGAGTTTGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAAC"
    ks = kmersList(K, seq, False)
    x = ks[0]
    xs = dict([(x, 1)])
    p = interpolate(K, xs, x, x, 1)
    assert p is not None
    assert len(p) == 1
    assert p[0] == x

def test_interpolate1Exacxt() :
    K = 25
    seq = "GGAGTTTCCAAGAGAAAATTTAGAGTTTGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAAC"
    ks = kmersList(K, seq, False)[:2]
    x = ks[0]
    y = ks[1]
    xs = dict([(z, 1) for z in ks])
    p = interpolate(K, xs, x, y, 2)
    assert p is not None
    assert len(p) == 2
    assert p[0] == x
    assert p[1] == y

def test_interpolate2Exacxt() :
    K = 25
    seq = "GGAGTTTCCAAGAGAAAATTTAGAGTTTGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAAC"
    ks = kmersList(K, seq, False)
    x = ks[0]
    y = ks[-1]
    xs = dict([(z, 1) for z in ks])
    p = interpolate(K, xs, x, y, len(ks))
    assert p is not None
    assert len(p) == len(ks)
    for i in range(len(ks)):
        assert ks[i] == p[i]

def test_interpolate3Exacxt() :
    K = 25
    seq = "GGAGTTTCCAAGAGAAAATTTAGAGTTTGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAAC"
    ks = kmersList(K, seq, False)
    x = ks[0]
    y = ks[-1]
    xs = dict([(z, 1) for z in ks])
    p = interpolate(K, xs, x, y, slice(len(ks)+1))
    assert p is not None
    assert len(p) == len(ks)
    for i in range(len(ks)):
        assert ks[i] == p[i]

