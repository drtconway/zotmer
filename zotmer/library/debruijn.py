"""
A collection of de Buijn graph manipulation tools.
"""

import random

import pytest

from pykmer.basics import ham, kmersList, render

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

def less_than_or_equal_max(i, n):
    if isinstance(n, int):
        return i <= n
    assert isinstance(n, slice)
    if n.stop is not None and i > n.stop:
        return False
    return True

def test_less_than_or_equal_max_0():
    n = 0
    i = 0
    assert less_than_or_equal_max(i, n)

def test_less_than_or_equal_max_1():
    n = 1
    i = 0
    assert less_than_or_equal_max(i, n)

def test_less_than_or_equal_max_2():
    n = 10
    for i in range(n):
        assert less_than_or_equal_max(i, n)

def test_less_than_or_equal_max_3():
    n = slice(10)
    for i in range(10):
        assert less_than_or_equal_max(i, n)

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
        while less_than_or_equal_max(i, n):
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

    def manyPaths(self, xbs, xes, n):
        fwd = {}
        for xb in xbs:
            fwd[xb] = [[xb]]
        rev = {}
        for xe in xes:
            rev[xe] = [[xe]]

        i = 1
        while less_than_or_equal_max(i, n):
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

def manyPaths(K, xs, xbs, xes, n):
    """
    Find paths through the de bruijn graph implied by xs, beginning
    at a member of xbs and ending at an element of xes. n specifies
    the length of the path. If it is an integer, then it gives an
    exact length of the path, otherwise it must be a slice which
    defines a half open interval within which the length of the
    path must lie.  If such a path exists, return the sequence of
    k-mers that form the path, including the endpoints. NB if xb
    == xe, then the path will be
    length 1.
    """
    I = Interpolator(K, xs)
    for p in I.manyPaths(xbs, xes, n):
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

def junction_kmers(K, st, en):
    """
    Generate a list of the k-mers that form a de Bruijn path corresponding
    to the sequence of st concatenated with the sequence of en.
    """
    M = (1 << (2*K)) - 1
    r = []
    for i in range(K+1):
        j = K - i
        x = st << (2*i)
        y = en >> (2*j)
        z = (x | y) & M
        r.append(z)
    return r

def test_junction_kmers() :
    K = 25
    seq = 'TACTTGCACTGGGAGGCACAGCGGCTTTTCAGTGTCACAGGTATTACGAG'
    xs = kmersList(K, seq)
    assert len(xs) == K+1
    ys = junction_kmers(K, xs[0], xs[-1])
    assert xs == ys

def debruijn_intersection(K, xs, ys):
    """
    For two lists of k-mers, xs and ys, return corresponding lists vs and ws
    such that vs contain only those elements of xs which have a de Bruijn
    neighbour in ys, and ws contains only those elements of ys that have a
    de Bruijn neighbour in xs.
    """
    K1 = K - 1
    M1 = (1 << (2*K1)) - 1

    # Project xs on to the trailing K-1 bases
    xIdx = {}
    for x in xs:
        x1 = x & M1
        if x1 not in xIdx:
            xIdx[x1] = []
        xIdx[x1].append(x)

    # Project ys on to the leading K-1 bases
    yIdx = {}
    for y in ys:
        y1 = y >> 2
        if y1 not in yIdx:
            yIdx[y1] = []
        yIdx[y1].append(y)

    vs = []
    ws = []
    for z1 in set(xIdx.keys()) & set(yIdx.keys()):
        vs += xIdx[z1]
        ws += yIdx[z1]

    return (vs, ws)

class fixed_path_finder(object):
    def __init__(self, K, D, X):
        self.K = K
        self.M = (1 << (2*K)) - 1
        self.D = D
        self.X = X

        self.lhsMasks = []
        self.rhsMasks = []
        for i in range(1,self.K+1):
            l = (1 << (2*i)) - 1
            r = l ^ self.M
            self.lhsMasks.append(l)
            self.rhsMasks.append(r)

    def lhs_mask(self, i, L):
        if i < self.K:
            return self.lhsMasks[i]
        return self.M

    def rhs_mask(self, i, L):
        j = i - (L - self.K + 1)
        if j < 0:
            return self.M
        assert 0 <= j and j < len(self.rhsMasks)
        return self.rhsMasks[j]

    def get_mask(self, i, L):
        return self.lhs_mask(i, L) & self.rhs_mask(i, L)

    def path_kmers(self, xs):
        L = len(xs)

        # Step 1: gather up all likely looking k-mers
        #
        Ys = []
        for i in range(L):
            x = xs[i]
            m = self.get_mask(i, L)
            x0 = x & m
            ys = []
            for y in self.X.iterkeys():
                if (y & m) != x0:
                    continue
                d = ham(x, y)
                if d < self.D:
                    ys.append(y)
            if len(ys) == 0:
                return None
            Ys.append(ys)

        # Step 2: eliminate k-mers that don't form part
        # of a complete path from xs[0]--xs[-1]
        #
        done = False
        while not done:
            done = True
            for j in range(1, L):
                i = j - 1
                vs0 = Ys[i]
                ws0 = Ys[j]
                (vs1, ws1) = debruijn_intersection(self.K, vs0, ws0)
                if len(vs1) == 0 or len(ws1) == 0:
                    return None
                if len(vs0) != len(vs1):
                    Ys[i] = vs1
                    done = False
                if len(ws0) != len(ws1):
                    Ys[j] = ws1
                    done = False

        return Ys

def fixed_path_kmers(K, D, X, xs):
    finder = fixed_path_finder(K, D, X)
    return finder.path_kmers(xs)

def test_fixed_path_kmers_0() :
    K = 25
    D = 3
    seq = 'TACTTGCACTGGGAGGCACAGCGGCTTTTCAGTGTCACAGGTATTACGAG'
    xs = kmersList(K, seq)
    X = dict([(x,1) for x in xs])
    Y = fixed_path_kmers(K, D, X, xs)
    assert Y is not None
    assert len(Y) == len(xs)
    for i in range(len(Y)):
        assert len(Y[i]) == 1
        assert Y[i][0] == xs[i]

def test_fixed_path_kmers_1() :
    random.seed(17)
    K = 25
    D = 3
    N = 100
    e = 0.01
    alts = {'A':['C','G','T'], 'C':['A','G','T'], 'G':['A','C','T'], 'T':['A','C','G']}
    seq = 'TACTTGCACTGGGAGGCACAGCGGCTTTTCAGTGTCACAGGTATTACGAG'
    L = len(seq)
    xs = kmersList(K, seq)
    X = {}
    for i in range(N):
        r = []
        for j in range(L):
            b = seq[j]
            if random.random() < e:
                b = random.choice(alts[b])
            r.append(b)
        s = ''.join(r)
        ys = kmersList(K, s)
        for y in ys:
            if y not in X:
                X[y] = 0
            X[y] += 1
    Y = fixed_path_kmers(K, D, X, xs)
    assert Y is not None
    assert len(Y) == len(xs)
    for i in range(len(Y)):
        assert xs[i] in Y[i]
        V = [(render(K, y), X[y]) for y in Y[i] if y != xs[i] and X[y] > X[xs[i]]]
        assert len(V) == 0

def test_fixed_path_kmers_2() :
    random.seed(17)
    K = 25
    D = 3
    N = 100
    e = 0.01
    seq0 = 'TACTTGCACTGGGAGGCACAGCGGCTTTTCAGTGTCACAGGTATTACGAG'
    seq1 = 'TACTTGCACTGGGAGGCCCAGCGGCTTTTCAGTGTCACAGGTATTAGGAG'
    xs = kmersList(K, seq0)
    ys = kmersList(K, seq1)
    X = dict([(y, 1) for y in ys])
    Y = fixed_path_kmers(K, D, X, xs)
    assert Y is not None
    assert len(Y) == len(xs)
    for i in range(len(Y)):
        assert len(Y[i]) == 1
        assert ys[i] in Y[i]
