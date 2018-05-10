import math

import pytest
import inspect
import random

def sgn(x):
    if x < 0.0:
        return -1
    return 1

def closeEnough(x, y):
    eps = 1.1920928955078125e-07

    a = math.fabs(x - y)
    if a < eps:
        return True

    if sgn(x) != sgn(y):
        return x == y

    (xm, xe) = math.frexp(x)
    (ym, ye) = math.frexp(y)

    if xe != ye:
        return False

    b = math.fabs(xm - ym)
    return b < eps

def setRndSeed():
    callerName = inspect.stack()[1][3]
    h = hash(callerName)
    random.seed(h)
    print callerName, h

def sqr(x):
    return x*x

class matrix(object):

    def __init__(self, kind, N, M):
        if M is None:
            M = N
        self.kind = kind
        self.N = N
        self.M = M

    def dim(self):
        return (self.N, self.M)

    def __iter__(self):
        for ij in self.iterKeys():
            yield (ij, self[ij])

    def __add__(self, other):
        if type(self) == type(other):
            r = self.copy()
        else:
            (N,M) = self.dim()
            r = dense(N, M, 0.0)
            for ij in self.iterKeys():
                r[ij] = self[ij]
        r += other
        return r

    def __iadd__(self, other):
        assert self.dim() == other.dim()
        for ij in other.iterKeys():
            self[ij] += other[ij]
        return self

    def __sub__(self, other):
        if type(self) == type(other):
            r = self.copy()
        else:
            (N,M) = self.dim()
            r = dense(N, M, 0.0)
            for ij in self.iterKeys():
                r[ij] = self[ij]
        r -= other
        return r

    def __isub__(self, other):
        print 'isub'
        assert self.dim() == other.dim()
        for ij in other.iterKeys():
            print ij
            self[ij] -= other[ij]
        return self

    def __mul__(self, other):
        (Na, Ma) = self.dim()
        (Nb, Mb) = other.dim()
        assert Ma == Nb

        r = dense(Na, Mb, 0.0)
        for i in xrange(Na):
            for j in xrange(Mb):
                for k in xrange(Ma):
                    r[(i,j)] += self[(i,k)] * other[(k,j)]

        return r

    def __str__(self):
        (N,M) = self.dim()
        return '\n'.join(['\t'.join([str(self[(i,j)]) for j in xrange(M)]) for i in xrange(N)])

def test_matrix0():
    N0 = 12
    M0 = 15
    m = matrix('v', N0, M0)
    d = m.dim()
    assert isinstance(d, tuple)
    assert len(d) == 2
    (N1, M1) = d
    assert N1 == N0
    assert M1 == M0
    assert m.kind == 'v'

class dense(matrix):
    def __init__(self, N, M = None, val = None):
        super(dense, self).__init__('d', N, M)

        if isinstance(val, matrix):
            self.data = [[val[(i,j)] for j in xrange(M)] for i in xrange(N)]
        else:
            if val is None:
                val = float('nan')
            elif isinstance(val, int):
                # coerce ints to floats to avoid all kinds of grief
                val = float(val)
            self.data = [[val for j in xrange(M)] for i in xrange(N)]

    def __len__(self):
        (N,M) = self.dim()
        return N*M

    def __getitem__(self, ij):
        (i,j) = ij
        #(N,M) = self.dim()
        #assert 0 <= i and i < N
        #assert 0 <= j and j < M
        return self.data[i][j]

    def __setitem__(self, ij, v):
        (i,j) = ij
        #(N,M) = self.dim()
        #assert 0 <= i and i < N
        #assert 0 <= j and j < M
        self.data[i][j] = v

    def copy(self):
        (N,M) = self.dim()
        return dense(N, M, self)

    def iterKeys(self):
        (N,M) = self.dim()
        for i in xrange(N):
            for j in xrange(M):
                yield (i,j)

    def iterRow(self, i):
        (N,M) = self.dim()
        assert 0 <= i and i < N
        for j in xrange(M):
            yield (i,j)

    def iterCol(self, j):
        (N,M) = self.dim()
        assert 0 <= j and j < N
        for i in xrange(N):
            yield (i,j)

    def transpose(self):
        (N, M) = self.dim()
        r = dense(M, N)
        for i in xrange(N):
            for j in xrange(M):
                r[(j,i)] = self[(i,j)]
        return r

def test_dense0():
    N = 3
    M = 2
    m = dense(N, M, 0.0)
    (N1, M1) = m.dim()
    assert N1 == N
    assert M1 == M
    for i in xrange(N):
        for j in xrange(M):
            assert m[(i,j)] == 0

def test_dense1():
    N = 3
    M = 2
    m = dense(N, M)
    (N1, M1) = m.dim()
    assert N1 == N
    assert M1 == M
    for i in xrange(N):
        for j in xrange(M):
            assert math.isnan(m[(i,j)])

def test_dense2():
    N = 3
    M = 2
    m = dense(N, M)
    (N1, M1) = m.dim()
    assert N1 == N
    assert M1 == M
    for i in xrange(N):
        for j in xrange(M):
            m[(i,j)] = 1.0
            assert m[(i,j)] == 1.0

def test_dense3():
    N = 3
    M = 2
    m = dense(N, M)
    (N1, M1) = m.dim()
    assert N1 == N
    assert M1 == M
    for i in xrange(N):
        for j in xrange(M):
            m[(i,j)] = 1.0
            assert m[(i,j)] == 1.0

def test_dense4():
    N = 3
    M = 2
    m = dense(N, M, 1.0)
    n = m.copy()
    assert n is not m
    assert n.dim() == m.dim()
    for i in xrange(N):
        for j in xrange(M):
            assert n[(i,j)] == m[(i,j)]

def test_dense5():
    N = 3
    M = 2
    m = dense(N, M, 1.0)
    ks = list(m.iterKeys())
    assert len(ks) == N*M
    ks = set(ks)
    assert len(ks) == N*M
    for i in xrange(N):
        for j in xrange(M):
            assert (i,j) in ks

def test_dense6():
    N = 3
    M = 2
    m = dense(N, M, 1.0)
    for i in xrange(N):
        rx = list(m.iterRow(i))
        assert len(rx) == M
        for j in xrange(M):
            assert (i,j) in rx
    for j in xrange(M):
        cx = list(m.iterCol(j))
        assert len(cx) == N
        for i in xrange(N):
            assert (i,j) in cx

def test_dense7():
    N = 3
    M = 2
    m = dense(N, M, 1.0)
    for i in xrange(N):
        for j in xrange(M):
            m[(i,j)] = i*M + j
    n = m.transpose()
    (Nt, Mt) = n.dim()
    assert Nt == M
    assert Mt == N
    for i in xrange(N):
        for j in xrange(M):
            assert n[(j,i)] == m[(i,j)]

def test_dense8():
    setRndSeed()
    N = 3
    M = 2
    m = dense(N, M)
    n = dense(N, M)
    for i in xrange(N):
        for j in xrange(M):
            m[(i,j)] = random.random()
            n[(i,j)] = random.random()

    r = m + n
    (N1, M1) = r.dim()
    assert N1 == N
    assert M1 == M
    for i in xrange(N):
        for j in xrange(M):
            ij = (i,j)
            assert r[ij] == m[ij] + n[ij]

def test_dense9():
    setRndSeed()
    N = 3
    M = 2
    m = dense(N, M)
    n = dense(N, M)
    for i in xrange(N):
        for j in xrange(M):
            m[(i,j)] = random.random()
            n[(i,j)] = random.random()

    r = m - n
    (N1, M1) = r.dim()
    assert N1 == N
    assert M1 == M
    for i in xrange(N):
        for j in xrange(M):
            ij = (i,j)
            assert r[ij] == m[ij] - n[ij]

def test_dense10():
    setRndSeed()
    N = 2
    M = 2
    m = dense(N, M)
    n = dense(N, M)
    m[(0,0)] = 1
    m[(0,1)] = 2
    m[(1,0)] = 3
    m[(1,1)] = 4
    n[(0,0)] = 5
    n[(0,1)] = 6
    n[(1,0)] = 7
    n[(1,1)] = 8
    r = m * n
    assert r[(0,0)] == 19
    assert r[(0,1)] == 22
    assert r[(1,0)] == 43
    assert r[(1,1)] == 50

class diagonal(matrix):
    def __init__(self, N, val = None):
        super(diagonal, self).__init__('i', N, N)
        if isinstance(val, diagonal):
            self.data = [val[(i,i)] for i in xrange(N)]
        else:
            if val is None:
                val = float('nan')
            elif not isinstance(val, float):
                val = float(val)
            self.data = [val for i in xrange(N)]

    def __len__(self):
        (N,M) = self.dim()
        assert N == M
        return N

    def __getitem__(self, ij):
        (N,M) = self.dim()
        assert N == M
        (i,j) = ij
        assert 0 <= i and i < N
        if i != j:
            return 0.0
        return self.data[i]

    def __setitem__(self, ij, v):
        (N,M) = self.dim()
        assert N == M
        (i,j) = ij
        assert i == j
        assert 0 <= i and i < N
        self.data[i] = v

    def copy(self):
        (N,M) = self.dim()
        assert N == M
        return diagonal(N, self)

    def iterKeys(self):
        (N,M) = self.dim()
        assert N == M
        for i in xrange(N):
            yield (i,i)

    def iterRow(self, i):
        (N,M) = self.dim()
        assert N == M
        assert 0 <= i and i < N
        yield (i,i)

    def iterCol(self, j):
        (N,M) = self.dim()
        assert N == M
        assert 0 <= j and j < N
        yield (j,j)

    def transpose(self):
        return self.copy()

def test_diagonal0():
    N = 5
    m = diagonal(N, 1.0)
    for i in xrange(N):
        for j in xrange(N):
            if i == j:
                assert m[(i,j)] == 1.0
            else:
                assert m[(i,j)] == 0.0

def test_diagonal1():
    N = 5
    m = diagonal(N, 1.0)
    n = dense(N, N, 0.0)
    for i in xrange(N):
        for j in xrange(N):
            if i == j:
                n[(i,j)] = 1.0

    for i in xrange(N):
        for j in xrange(N):
            ij = (i,j)
            assert m[ij] == n[ij]

    r1 = m + m
    assert isinstance(r1, diagonal)
    r2 = n + n
    assert isinstance(r2, dense)

    for i in xrange(N):
        for j in xrange(N):
            ij = (i,j)
            assert r1[ij] == r2[ij]

def test_diagonal2():
    setRndSeed()
    N = 5
    m = diagonal(N, 1.0)
    n = dense(N, N, 0.0)
    for i in xrange(N):
        for j in xrange(N):
            n[(i,j)] = random.random()
    r = m + n
    assert isinstance(r, dense)
    for i in xrange(N):
        for j in xrange(N):
            ij = (i,j)
            assert r[ij] == m[ij] + n[ij]

def test_diagonal3():
    setRndSeed()
    N = 5
    m = diagonal(N, 1.0)
    n = dense(N, N, 0.0)
    for i in xrange(N):
        for j in xrange(N):
            n[(i,j)] = random.random()
    r = m - n
    assert isinstance(r, dense)
    for i in xrange(N):
        for j in xrange(N):
            ij = (i,j)
            assert r[ij] == m[ij] - n[ij]

class sparse(matrix):
    def __init__(self, N, M, val = None):
        super(sparse, self).__init__('s', N, M)
        self.data = {}
        if val is not None:
            for (ij, v) in val:
                if v == 0.0:
                    continue
                self.data[ij] = v

    def __len__(self):
        return len(self.data)

    def __getitem__(self, ij):
        (N,M) = self.dim()
        (i,j) = ij
        assert 0 <= i and i < N
        assert 0 <= j and j < M
        if ij not in self.data:
            return 0.0
        return self.data[ij]

    def __setitem__(self, ij, v):
        (N,M) = self.dim()
        (i,j) = ij
        assert 0 <= i and i < N
        assert 0 <= j and j < M
        if v == 0.0 and ij in self.data:
            del self.data[ij]
        else:
            self.data[ij] = v

    def __mul__(self, other):
        (Na, Ma) = self.dim()
        (Nb, Mb) = other.dim()
        assert Ma == Nb

        r = dense(Na, Mb, 0.0)
        for ik in self.iterkeys():
            (i,k) = ik
            for kj in other.iterRow(k):
                assert kj[0] == k
                j = kj[1]
                r[(i,j)] += self[(i,k)] * other[(k,j)]
        return r

    def copy(self):
        (N,M) = self.dim()
        return sparse(N, M, self)

    def iterKeys(self):
        return self.data.iterkeys()

    def iterRow(self, I):
        (N,M) = self.dim()
        assert 0 <= I and I < N
        for ij in self.data.iterkeys():
            (i,j) = ij
            if i == I:
                yield ij

    def iterCol(self, J):
        (N,M) = self.dim()
        assert 0 <= J and J < M
        for ij in self.data.iterkeys():
            (i,j) = ij
            if j == J:
                yield ij

    def transpose(self):
        (N, M) = self.dim()
        r = sparse(N, M)
        for (ij, v) in self:
            (i,j) = ij
            r[(j,i)] = v
        return r

    def density(self):
        (N, M) = self.dim()
        return float(len(self.data)) / float(N*M)

def test_sparse0():
    setRndSeed()
    N = 1000
    m = sparse(N, N)
    (N1, M1) = m.dim()
    assert N1 == N
    assert M1 == N
    J = 100
    seen = {}
    for n in xrange(J):
        i = random.randint(0, N-1)
        j = random.randint(0, N-1)
        ij = (i,j)
        while ij in seen:
            i = random.randint(0, N-1)
            j = random.randint(0, N-1)
            ij = (i,j)
        seen[ij] = float(n)
        m[ij] = float(n)

    ijs = list(m.iterKeys())
    assert len(ijs) == J
    ijs = set(ijs)
    assert len(ijs) == J
    for ij in ijs:
        assert ij in seen
        assert m[ij] == seen[ij]

class lowerTri(matrix):
    def __init__(self, N, M = None, val = None, trans = None):
        super(lowerTri, self).__init__('l', N, M)
        (n,m) = self.dim()
        assert val is None or trans is None
        if isinstance(val, lowerTri):
            self.data = [[val[(i,j)] for j in xrange(min(i+1, M))] for i in xrange(N)]
        elif isinstance(trans, upperTri):
            self.data = [[val[(j,i)] for j in xrange(min(i+1, M))] for i in xrange(N)]
        else:
            if val is None:
                val = float('nan')
            elif not isinstance(val, float):
                # coerce ints to floats to avoid all kinds of grief
                val = float(val)
            self.data = [[val for j in xrange(min(i+1, M))] for i in xrange(N)]

    def __len__(self):
        (N, M) = self.dim()
        l = min(N, M)
        return l * (l+1) // 2

    def __getitem__(self, ij):
        (i,j) = ij
        #(N, M) = self.dim()
        #assert 0 <= i and i < N
        #assert 0 <= j and j < M
        if j > i:
            return 0.0
        return self.data[i][j]

    def __setitem__(self, ij, v):
        (i,j) = ij
        #(N, M) = self.dim()
        #assert 0 <= i and i < N
        #assert 0 <= j and j < M
        assert j <= i
        self.data[i][j] = v

    def __mul__(self, other):
        if not isinstance(other, lowerTri) and not isinstance(other, diagonal):
            return super(lowerTri, self).__mul__(other)

        (Na, Ma) = self.dim()
        (Nb, Mb) = other.dim()
        assert Ma == Nb

        r = lowerTri(Na, Mb, 0.0)
        for i in xrange(Na):
            for j in xrange(min(i+1, Mb)):
                for k in xrange(min(i+1, Ma)):
                    r[(i,j)] += self[(i,k)] * other[(k,j)]

        return r

    def iterKeys(self):
        (N, M) = self.dim()
        for i in xrange(N):
            for j in xrange(min(i+1, M)):
                yield (i,j)

    def iterRow(self, i):
        (N, M) = self.dim()
        assert 0 <= i and i < N
        for j in xrange(min(i+1, M)):
            yield (i,j)

    def iterCol(self, j):
        (N, M) = self.dim()
        assert 0 <= j and j < M
        for i in xrange(j, N):
            yield (i,j)

    def copy(self):
        (N,M) = self.dim()
        return lowerTri(N, M, self)

    def transpose(self):
        (N, M) = self.dim()
        return upperTri(M, N, val=None, trans=self)

    def fwdSub(self, b):
        (N, M) = self.dim()
        assert N == M
        assert len(b) == N
        x = [0.0 for i in xrange(N)]
        for i in xrange(N):
            t = b[i]
            for j in xrange(i):
                t -= self[(i,j)] * x[j]
            x[i] = t / self[(i,i)]
        return x

def test_lowerTri0():
    N = 11
    m = lowerTri(N, N, 1.0)
    for i in xrange(N):
        for j in xrange(N):
            if j <= i:
                assert m[(i,j)] == 1.0
            else:
                assert m[(i,j)] == 0.0

def test_lowerTri1():
    N = 11
    m = lowerTri(N, N, 1.0)
    n = m + m
    assert isinstance(n, lowerTri)
    for i in xrange(N):
        for j in xrange(N):
            if j <= i:
                assert n[(i,j)] == 2.0
            else:
                assert n[(i,j)] == 0.0

def test_lowerTri2():
    N = 11
    m = lowerTri(N, N, 1.0)
    n = m * m
    assert isinstance(n, lowerTri)
    print n
    for i in xrange(N):
        for j in xrange(N):
            if j <= i:
                assert n[(i,j)] == 1.0 + i - j
            else:
                assert n[(i,j)] == 0.0

class upperTri(matrix):
    def __init__(self, N, M = None, val = None, trans = None):
        super(upperTri, self).__init__('u', N, M)
        (n,m) = self.dim()
        assert val is None or trans is None
        if isinstance(val, upperTri):
            self.data = [[val[(i,j)] for j in xrange(i, M)] for i in xrange(N)]
        elif isinstance(trans, lowerTri):
            self.data = [[trans[(j,i)] for j in xrange(i, M)] for i in xrange(N)]
        else:
            if val is None:
                val = float('nan')
            elif not isinstance(val, float):
                # coerce ints to floats to avoid all kinds of grief
                val = float(val)
            self.data = [[val for j in xrange(i, M)] for i in xrange(N)]

    def __len__(self):
        (N, M) = self.dim()
        l = min(N, M)
        return l * (l+1) // 2

    def __getitem__(self, ij):
        (i,j) = ij
        #(N, M) = self.dim()
        #assert 0 <= i and i < N
        #assert 0 <= j and j < M
        if j < i:
            return 0.0
        return self.data[i][j-i]

    def __setitem__(self, ij, v):
        (i,j) = ij
        #(N, M) = self.dim()
        #assert 0 <= i and i < N
        #assert 0 <= j and j < M
        assert j >= i
        self.data[i][j-i] = v

    def __mul__(self, other):
        if not isinstance(other, upperTri) and not isinstance(other, diagonal):
            return super(upperTri, self).__mul__(other)

        (Na, Ma) = self.dim()
        (Nb, Mb) = other.dim()
        assert Ma == Nb

        r = upperTri(Na, Mb, 0.0)
        for i in xrange(Na):
            for j in xrange(i, Mb):
                for k in xrange(i, Ma):
                    r[(i,j)] += self[(i,k)] * other[(k,j)]
        return r

    def iterKeys(self):
        (N, M) = self.dim()
        for j in xrange(M):
            for i in xrange(min(j+1, N)):
                yield (i,j)

    def iterRow(self, i):
        (N, M) = self.dim()
        assert 0 <= i and i < N
        for j in xrange(i, M):
            yield (i,j)

    def iterCol(self, j):
        (N, M) = self.dim()
        assert 0 <= j and j < M
        for i in xrange(min(j+1, N)):
            yield (i,j)

    def copy(self):
        (N,M) = self.dim()
        return upperTri(N, M, self)

    def transpose(self):
        (N, M) = self.dim()
        return lowerTri(M, N, val=None, trans=self)

    def backSub(self, b):
        (N, M) = self.dim()
        assert N == M
        assert len(b) == N
        x = [0.0 for i in xrange(N)]
        Nm1 = N - 1
        for i0 in xrange(N):
            i = Nm1 - i0
            t = b[i]
            for j in xrange(i, N):
                t -= self[(i,j)] * x[j]
            x[i] = t / self[(i,i)]
        return x

def cholesky(A):
    (N,M) = A.dim()
    assert N == M

    L = dense(N, N, 0.0)
    for j in xrange(N):
        t = A[(j,j)]
        for k in xrange(j):
            t -= sqr(L[(j,k)])
        print t
        L[(j,j)] = math.sqrt(t)

        for i in xrange(j+1, N):
            t = A[(i,j)]
            for k in xrange(j):
                t -= L[(i,k)]*L[(j,k)]
            L[(i,j)] = t / L[(j,j)]
    return L

def choleskyLDL(A):
    (N,M) = A.dim()
    assert N == M

    D = diagonal(N, 0.0)
    L = lowerTri(N, N, 0.0)

    for j in xrange(N):
        t = A[(j,j)]
        for k in xrange(j):
            t -= sqr(L[(j,k)])*D[(k,k)]
        D[(j,j)] = t

        L[(j,j)] = 1.0
        for i in xrange(j+1, N):
            t = A[(i,j)]
            for k in xrange(j):
                t -= L[(i,k)]*L[(j,k)]*D[(k,k)]
            L[(i,j)] = t / D[(j,j)]
    return (L,D)

