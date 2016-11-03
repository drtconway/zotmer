from base import Cmd

import pykmer.kfset as kfset
from pykmer.container import probe

import math

def sqr(x):
    return x * x

def gauss(x0, x, n):
    return math.exp(-sqr(x0 - x)/(2.0*n*n))

def smooth(h):
    h1 = {}
    for x0 in h.keys():
        num = 0
        den = 0
        for (xi,yi) in h.items():
            z = gauss(x0, xi, 1.5)
            #print x0, xi, z, ly
            num += z*yi
            den += z
        y0 = num/den
        #print x0, h.get(x0, 0), y0
        h1[x0] = y0
    return h1

def infer(fn):
    (m, xs) = kfset.read(fn)
    K = m['K']
    h = {}
    for (x, f) in xs:
        c = h.get(f, 0)
        h[f] = c+1

    h1 = smooth(h)
    h1 = h1.items()
    h1.sort()
    v = h1[0]
    for v1 in h1[1:]:
        if v1[1] < v[1]:
            v = v1
        else:
            break
    return v[0]

def trim(xs, c):
    for (x,f) in xs:
        if f >= c:
            yield(x, f)

class Trim(Cmd):
    def run(self, opts):
        inp = opts['<kmers>']
        out = opts['<output>']
        c = int(opts['<cutoff>'])
        if c == 0:
            c = infer(inp)
            print 'inferred cutoff:', c
        (m, xs) = kfset.read(inp)
        K = m['K']
        kfset.write(K, trim(xs, c), out, m)

def add(cmds):
    cmds['trim'] = Trim()

