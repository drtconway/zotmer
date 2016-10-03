from base import Cmd

from pykmer.basics import kmersWithPos, lcp, render
from pykmer.file import readFasta
import pykmer.kset as kset
import pykmer.kfset as kfset
from pykmer.stats import logChi2CDF, counts2cdf, ksDistance2
from pykmer.exceptions import MismatchedK

import math
import sys

def getK(ins):
    k = None
    for fn in ins:
        k0 = kfset.probeK(fn)
        if k is None:
            k = k0
        elif k != k0:
            raise MismatchedK(k, k0)
    return k

def findLcp(idx, x):
    l = 0
    h = len(idx)
    while h > l:
        m = (h + l) // 2
        y = idx[m][0]
        if y == x:
            return m
        if y < x:
            l = m + 1
        else:
            h = m - 1
    return l

def resolveLcp(k, idx, lcps, x):
    i = findLcp(idx, x)
    z = min(i+1, len(idx))
    for h in range(z):
        j = z - h - 1
        l = lcp(k, x, idx[j][0])
        if l > lcps[j][0]:
            lcps[j] = (l, x)
        else:
            break
    for j in range(i + 1, len(idx)):
        l = lcp(k, x, idx[j][0])
        if l > lcps[j][0]:
            lcps[j] = (l, x)
        else:
            break

def gcat(k, x):
    cs = [0, 0, 0, 0]
    for i in range(k):
        cs[x&3] += 1
        x >>= 2
    return (float(cs[0] + cs[3])/k, float(cs[1] + cs[2])/k)

def null(g, s, j):
    if j == 0:
        return 0.0
    return math.exp(s*math.log1p(-math.pow(g,j)))

class Vars(Cmd):
    def run(self, opts):
        refFn = opts['<kmers>']
        k = getK(opts['<input>'])

def add(cmds):
    cmds['vars'] = Vars()
