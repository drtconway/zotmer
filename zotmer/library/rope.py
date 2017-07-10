import pytest

class rope(object):

    @staticmethod
    def atom(txt):
        return rope_atom(txt)

    @staticmethod
    def substr(parent, start, stop):
        return rope_substr(parent, start, stop).simplify()

    @staticmethod
    def concat(lhs, rhs):
        return rope_concat(lhs, rhs).simplify()

    @staticmethod
    def split(parent, idx):
        lhs = substr(parent, 0, idx)
        rhs = substr(parent, idx, len(parent))
        return (lhs, rhs)

    @staticmethod
    def join(itms):
        assert len(itms) > 0
        if len(itms) == 1:
            return itms[0]
        else:
            m = len(itms) // 2
            l = rope.join(itms[:m])
            r = rope.join(itms[m:])
            return rope.concat(l, r)

def test_atom0a():
    s = ''
    a = rope.atom(s)
    assert len(a) == 0
    assert a[:] == ''

def test_atom0b():
    s = ''
    a = rope.atom(s)
    with pytest.raises(AssertionError):
        c = a[0]

def test_atom0c():
    s = ''
    a = rope.atom(s)
    with pytest.raises(AssertionError):
        c = a[1:5]

def test_atom1a():
    s = 'GAAGAGGTGAATGTAATTCCTCCACACACTCCAGTTAGGTATGAATTTTCCTACTTTTAATTATATTATAATTTTG'
    a = rope.atom(s)
    assert len(s) == len(a)
    assert a[:] == s

def test_atom1b():
    s = 'GAAGAGGTGAATGTAATTCCTCCACACACTCCAGTTAGGTATGAATTTTCCTACTTTTAATTATATTATAATTTTG'
    a = rope.atom(s)
    assert a[0] == 'G'
    assert a[0:4] == 'GAAG'
    assert a[:4] == 'GAAG'
    assert a[58:] == 'AATTATATTATAATTTTG'

class rope_atom(rope):
    def __init__(self, text):
        self.text = text

    def __len__(self):
        return len(self.text)

    def __getitem__(self, idx):
        assert not isinstance(idx, slice) or idx.start is None or idx.start >= 0
        assert not isinstance(idx, slice) or idx.stop is None or idx.stop <= len(self.text)
        assert not isinstance(idx, slice) or idx.step is None or idx.step == 1
        assert isinstance(idx, slice) or idx >= 0
        assert isinstance(idx, slice) or idx < len(self.text)
        return self.text[idx]

    def __repr__(self):
        if len(self.text) < 50:
            return self.text
        else:
            return '%s...%d...%s' % (self.text[:5], len(self.text) - 10, self.text[-5:])

    def simplify(self):
        return self


def test_substr0a():
    s = ''
    a = rope.atom(s)
    r = rope.substr(a, 0, 0)
    assert r[:] == s

def test_substr0b():
    s = 'GAAGAGGTGAATGTAATTCCTCCACACACTCCAGTTAGGTATGAATTTTCCTACTTTTAATTATATTATAATTTTG'
    a = rope.atom(s)
    r = rope.substr(a, 0, 0)
    assert r[:] == ''

def test_substr1a():
    s = 'GAAGAGGTGAATGTAATTCCTCCACACACTCCAGTTAGGTATGAATTTTCCTACTTTTAATTATATTATAATTTTG'
    a = rope.atom(s)
    r = rope.substr(a, 0, 1)
    assert r[:] == 'G'

def test_substr1b():
    s = 'GAAGAGGTGAATGTAATTCCTCCACACACTCCAGTTAGGTATGAATTTTCCTACTTTTAATTATATTATAATTTTG'
    a = rope.atom(s)
    r = rope.substr(a, 10, 20)
    assert r[:] == 'ATGTAATTCC'

def test_substr2a():
    s = 'GAAGAGGTGAATGTAATTCCTCCACACACTCCAGTTAGGTATGAATTTTCCTACTTTTAATTATATTATAATTTTG'
    a = rope.atom(s)
    r0 = rope.substr(a, 10, 30)
    assert r0[:] == 'ATGTAATTCCTCCACACACT'
    r1 = rope.substr(r0, 5, 15)
    assert r1[:] == 'ATTCCTCCAC'
    assert repr(r1) == 'GAAGA...66...TTTTG[15:25]'

class rope_substr(rope):
    def __init__(self, parent, start, stop):
        assert start >= 0
        assert stop <= len(parent)
        self.parent = parent
        self.start = start
        self.stop = stop

    def __len__(self):
        return self.stop - self.start

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            s = 0
            if idx.start is not None:
                s = idx.start
            e = (self.stop - self.start)
            if idx.stop is not None:
                e = idx.stop
            assert idx.step is None or idx.step == 1
            s += self.start
            e += self.start
            return self.parent[s:e]
        else:
            assert idx >= 0 and idx < (self.stop - self.start)
            i = self.start + idx
            return self.parent[i]

    def __repr__(self):
        return '%s[%d:%d]' % (repr(self.parent), self.start, self.stop)

    def simplify(self):
        self.parent = self.parent.simplify()
        if isinstance(self.parent, rope_substr):
            self.start += self.parent.start
            self.stop += self.parent.start
            self.parent = self.parent.parent
            return self
        if isinstance(self.parent, rope_concat):
            ll = len(self.parent.lhs)
            if self.start >= ll:
                s = self.start - ll
                e = self.stop - ll
                n = rope_substr(self.parent.rhs, s, e)
                return n.simplify()
            if self.stop <= ll:
                n = rope_substr(self.parent.lhs, self.start, self.stop)
                return n.simplify()
        return self

class rope_concat(rope):
    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    def __len__(self):
        return len(self.lhs) + len(self.rhs)

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            s = 0
            if idx.start is not None:
                s = idx.start
            e = len(self)
            if idx.stop is not None:
                e = idx.stop
            assert idx.step is None or idx.step == 1
            if e < len(self.lhs):
                return self.lhs[s:e]
            elif s >= len(self.lhs):
                s -= len(self.lhs)
                e -= len(self.lhs)
                return self.rhs[s:e]
            else:
                e -= len(self.lhs)
                return self.lhs[s:] + self.rhs[:e]
        elif idx < len(self.lhs):
            return self.lhs[idx]
        else:
            idx -= len(self.lhs)
            assert idx < len(self.rhs)
            return self.rhs[idx]

    def __repr__(self):
        return '(%s)(%s)' % (repr(self.lhs), repr(self.rhs))

    def simplify(self):
        return self
