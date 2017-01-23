import pykmer.container.std as std

def merge2(z, K, xs, ys, h, nm = None):
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
