from zotmer.library.files import countsWriter, kmerWriter, readKmersAndCounts, writeKmersAndCounts, writeKmersAndCounts2

def merge2(z, K, xs, ys, h, nm = None):
    if nm is not None:
        knm = nm + '-kmers'
        cnm = nm + '-counts'
    else:
        knm = None
        cnm = None

    with kmerWriter(z, knm) as kx, countsWriter(z, cnm) as cx:
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
                h[x[1]] = 1 + h.get(x[1], 0)
                kx.append(xk), cx.append(x[1])
                try:
                    x = xs.next()
                    xk = x[0]
                except StopIteration:
                    moreXs = False
                continue

            if xk > yk:
                h[y[1]] = 1 + h.get(y[1], 0)
                kx.append(yk), cx.append(y[1])
                try:
                    y = ys.next()
                    yk = y[0]
                except StopIteration:
                    moreYs = False
                continue

            assert xk == yk
            c = x[1] + y[1]
            h[c] = 1 + h.get(c, 0)
            kx.append(xk), cx.append(c)

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
            h[x[1]] = 1 + h.get(x[1], 0)
            kx.append(x[0]), cx.append(x[1])
            try:
                x = xs.next()
            except StopIteration:
                moreXs = False

        while moreYs:
            h[y[1]] = 1 + h.get(y[1], 0)
            kx.append(y[0]), cx.append(y[1])
            try:
                y = ys.next()
            except StopIteration:
                moreYs = False

def merge2Block(z, K, xs, ys, h, nm = None):
    if nm is not None:
        knm = nm + '-kmers'
        cnm = nm + '-counts'
    else:
        knm = None
        cnm = None

    with std.kmerWriter(z, K, knm) as kx, std.countsWriter(z, K, cnm) as cx:
        moreXs = True
        moreYs = True

        try:
            xsBlk = xs.next()
            xsZ = len(xsBlk)
            i = 0
        except StopIteration:
            moreXs = False

        try:
            ysBlk = ys.next()
            ysZ = len(ysBlk)
            j = 0
        except StopIteration:
            moreYs = False

        B = 65536
        zBlk = []
        cBlk = []
        n = 0
        while moreXs and moreYs:
            while i < xsZ and j < ysZ and n < B:
                if xsBlk[i][0] < ysBlk[j][0]:
                    h[xsBlk[i][1]] = 1 + h.get(xsBlk[i][1], 0)
                    zBlk.append(xsBlk[i][0])
                    cBlk.append(xsBlk[i][1])
                    i += 1
                    n += 1
                    continue
                if xsBlk[i][0] > ysBlk[j][0]:
                    h[ysBlk[j][1]] = 1 + h.get(ysBlk[j][1], 0)
                    zBlk.append(ysBlk[j][0])
                    cBlk.append(ysBlk[j][1])
                    j += 1
                    n += 1
                    continue

                c = xsBlk[i][1] + ysBlk[j][1]
                h[c] = 1 + h.get(c, 0)
                zBlk.append(xsBlk[i][0])
                cBlk.append(c)
                i += 1
                j += 1
                n += 1

            if i == xsZ:
                try:
                    xsBlk = xs.next()
                    xsZ = len(xsBlk)
                    i = 0
                except StopIteration:
                    moreXs = False

            if j == ysZ:
                try:
                    ysBlk = ys.next()
                    ysZ = len(ysBlk)
                    j = 0
                except StopIteration:
                    moreYs = False

            if n == B:
                kx.append(zBlk)
                cx.append(cBlk)
                zBlk = []
                cBlk = []
                n = 0

        while moreXs:
            while i < xsZ and n < B:
                h[xsBlk[i][1]] = 1 + h.get(xsBlk[i][1], 0)
                zBlk.append(xsBlk[i][0])
                cBlk.append(xsBlk[i][1])
                i += 1
                n += 1

            if i == xsZ:
                try:
                    xsBlk = xs.next()
                    xsZ = len(xsBlk)
                    i = 0
                except StopIteration:
                    moreXs = False

            if n == B:
                kx.append(zBlk)
                cx.append(cBlk)
                zBlk = []
                cBlk = []
                n = 0

        while moreYs:
            while j < ysZ and n < B:
                h[ysBlk[j][1]] = 1 + h.get(ysBlk[j][1], 0)
                zBlk.append(ysBlk[j][0])
                cBlk.append(ysBlk[j][1])
                j += 1
                n += 1

            if j == ysZ:
                try:
                    ysBlk = ys.next()
                    ysZ = len(ysBlk)
                    j = 0
                except StopIteration:
                    moreYs = False

            if n == B:
                kx.append(zBlk)
                cx.append(cBlk)
                zBlk = []
                cBlk = []
                n = 0

        if n > 0:
            kx.append(zBlk)
            cx.append(cBlk)

class _kmerStream(object):
    def __init__(self, xs):
        self.more = True
        self.xs = xs
        self.x = None

    def next(self):
        try:
            self.x = self.xs.next()
        except StopIteration:
            self.more = False

    def __lt__(self, other):
        assert self.more
        assert other.more
        return self.x[0] < other.x[0]

def mergeN(xss, hist):
    ss = [_kmerStream(xs) for xs in xss]
    q = heap(ss)
    while len(q) > 0:
        x = q.front()[0]
        c = 0
        while q.front()[0] == x:
            c += q.front()[1]
            q.front().next()
            if q.front().more:
                q.modifyfront()
            else:
                q.pop()
        hist[c] = 1 + hist.get(c, 0)
        yield (x, c)
