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

class rope_atom(rope):
    def __init__(self, text):
        self.text = text

    def __len__(self):
        return len(self.text)

    def __getitem__(self, idx):
        return self.text[idx]

    def __repr__(self):
        if len(self.text) < 50:
            return self.text
        else:
            return '%s...%d...%s' % (self.text[:5], len(self.text) - 10, self.text[-5:])

    def simplify(self):
        return self

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
