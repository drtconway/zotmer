import sys
from zotmer.library.basics import ham, render
from zotmer.library.misc import unionfind

def make_masks(K, tau):
    tau1 = tau + 1
    M = (1 << (2*K)) - 1
    msks = [0 for i in range(tau1)]
    for i in range(K):
        v = i % tau1
        w = 3 << (2*i)
        for j in range(tau1):
            if v == j:
                msks[j] |= w
    msks = [m^M for m in msks]
    return msks

def trim(K, kx):
    tau = 3
    msks = make_masks(K, tau)
    frac = 0.05

    vx = {}
    for x in kx.iterkeys():
        for m in msks:
            y = (x & m)
            if y not in vx:
                vx[y] = []
            vx[y].append(x)
    
    wx = unionfind()
    for xs in vx.itervalues():
        for i in xrange(len(xs)):
            x = xs[i]
            for j in xrange(i+1, len(xs)):
                y = xs[j]
                d = ham(x, y)
                if d <= tau:
                    wx.union(x, y)

    vx = {}
    for x in kx.iterkeys():
        a = wx.find(x)
        if a not in vx:
            vx[a] = []
        vx[a].append((kx[x], x))

    wx = set([])
    for xs in vx.values():
        xs.sort()
        i = 0
        while i < len(xs) and xs[i][0] < frac*xs[-1][0]:
            wx.add(xs[i][1])
            i += 1
    for x in wx:
        del kx[x]
