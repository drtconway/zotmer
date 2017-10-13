"""
A collection of de Buijn graph manipulation tools.
"""

class Interpolator1(object):
    def __init__(self, K, xs):
        self.K = K
        self.M = ((1 << (2*K)) - 1)
        self.S = 2*(K-1)
        self.xs = xs

    def succ(self, x):
        y0 = (x << 2) & self.M
        for i in range(4):
            y = y0 + i
            if self.xs.get(y, 0) > 10:
                yield y

    def pred(self, x):
        y0 = x >> 2
        for i in range(4):
            y = y0 + (i << self.S)
            if self.xs.get(y, 0) > 10:
                yield y

    def path(self, xb, xe, n):
        stk = [(xb, n+self.K, [])]
        while len(stk) > 0:
            (x0, n0, p0) = stk.pop()
            if n0 == 0:
                if x0 == xe:
                    yield p0
                continue
            n1 = n0 - 1
            for x1 in self.succ(x0):
                p1 = p0 + [x0]
                stk.append((x1, n1, p1))

    def pathTo(self, xb, xe, N):
        stk = [(xb, 0, [xb])]
        while len(stk) > 0:
            (x0, n0, p0) = stk.pop()
            if x0 == xe:
                yield p0
            n1 = n0 + 1
            if n1 == N:
                continue
            for x1 in self.succ(x0):
                p1 = p0 + [x1]
                stk.append((x1, n1, p1))

    def interpolate(self, xb, xe, n):
        sMax = 0
        pMax = None
        for p in self.path(xb, xe, n):
            s = 0
            for x in p:
                s += self.xs[x]
            if s > sMax:
                sMax = s
                pMax = p
        return pMax

class Interpolator2(object):
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
        fwd[xb] = []
        rev = {}
        rev[xe] = []

        for i in range(n+self.K):
            if (i & 1) == 0:
                nextFwd = {}
                for (x,p) in fwd.iteritems():
                    pp = p + [x]
                    for y in self.succ(x):
                        nextFwd[y] = pp
                fwd = nextFwd
            else:
                nextRev = {}
                for (x,p) in rev.iteritems():
                    pp = [x] + p
                    for y in self.pred(x):
                        nextRev[y] = pp
                rev = nextRev

        for x in set(fwd.keys()) & set(rev.keys()):
            yield fwd[x] + [x] + rev[x]

    def interpolate(self, xb, xe, n):
        sMax = 0
        pMax = None
        for p in self.path(xb, xe, n):
            s = 0
            for x in p:
                s += self.xs[x]
            if s > sMax:
                sMax = s
                pMax = p
        return pMax

    def pathBetween(self, xb, xe, lo, hi):
        fwd = {}
        fwd[xb] = []
        rev = {}
        rev[xe] = []

        for i in range(hi):
            if lo <= i and i < hi:
                for x in set(fwd.keys()) & set(rev.keys()):
                    yield fwd[x] + rev[x]

            if (i & 1) == 0:
                nextFwd = {}
                for (x,p) in fwd.iteritems():
                    pp = p + [x]
                    for y in self.succ(x):
                        nextFwd[y] = pp
                fwd = nextFwd
            else:
                nextRev = {}
                for (x,p) in rev.iteritems():
                    pp = [x] + p
                    for y in self.pred(x):
                        nextRev[y] = pp
                rev = nextRev

    def interpolateBetween(self, xb, xe, lo, hi):
        sMax = 0
        pMax = None
        for p in self.pathBetween(xb, xe, lo, hi):
            s = 0
            for x in p:
                s += self.xs[x]
            if s > sMax:
                sMax = s
                pMax = p
        return pMax

def interpolate(K, xs, xb, xe, n):
    I2 = Interpolator2(K, xs)
    p2 = I2.interpolate(xb, xe, n)
    return p2

def interpolateBetween(K, xs, xb, xe, lo, hi):
    I2 = Interpolator2(K, xs)
    p2 = I2.interpolateBetween(xb, xe, lo, hi)
    return p2

def pathBetween(K, xs, xb, xe, N):
    I1 = Interpolator1(K, xs)
    return I1.pathTo(xb, xe, N)

def onlyPathBetween(K, xs, xb, xe, N):
    res = None
    for pth in pathBetween(K, xs, xb, xe, N):
        if res is None:
            res = pth
        else:
            # multiple paths
            return None
    return res
