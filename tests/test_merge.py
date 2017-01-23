from zotmer.helpers.merge import merge2
from pykmer.file import tmpfile
from pykmer.container import container
from pykmer.container.std import readKmersAndCounts, writeKmersAndCounts

import math
import os
import random

def pois(lam):
    lam = float(lam)
    x = 0
    p = math.exp(-lam)
    s = p
    u = random.random()
    while u > s:
        x += 1
        p *= lam/x
        s += p
    return x

def test_merg2_0():
    K = 27
    M = (1 << (2*K)) - 1
    N = 100000
    random.seed(17)
    xs = [(random.randint(0, M), pois(10)) for i in xrange(N)]
    xs.sort()
    ys = [(random.randint(0, M), pois(10)) for i in xrange(N)]
    ys.sort()

    nm0 = tmpfile()
    with container(nm0, 'w') as z:
        writeKmersAndCounts(K, xs, z, 'xs')
        writeKmersAndCounts(K, ys, z, 'ys')
    
    nm1 = tmpfile()
    h = {}
    with container(nm0, 'r') as z0, container(nm1, 'w') as z:
        merge2(z, K, readKmersAndCounts(z0, 'xs'), readKmersAndCounts(z0, 'ys'), h, 'zs')
    h = h.items()
    h.sort()

    ws = {}
    for (x,c) in xs:
        ws[x] = c + ws.get(x, 0)
    for (y,c) in ys:
        ws[y] = c + ws.get(y, 0)
    ws = ws.items()
    ws.sort()

    with container(nm1, 'r') as z:
        zs = list(readKmersAndCounts(z, 'zs'))

    assert len(ws) == len(zs)
    for i in xrange(len(ws)):
        assert ws[i] == zs[i]

    h1 = {}
    for (_, c) in ws:
        h1[c] = 1 + h1.get(c, 0)
    h1 = h1.items()
    h1.sort()

    assert len(h) == len(h1)
    for i in xrange(len(h)):
        assert h[i] == h1[i]

