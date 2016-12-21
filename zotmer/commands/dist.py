"""
Usage:
    zot dist [-M measure]... <k> <input>...

Options:
    -M measure  use "measure" for the distance between k-mer frequency sets.
                Use "-M list" to get a list of available measures.
"""

import pykmer.kfset as kfset
import pykmer.dist as dist
from pykmer.exceptions import MismatchedK

import array
import docopt
import fnmatch
import sys

measures = {}

class Measure:
    def __init__(self, name, desc, vec, func):
        self.name = name
        self.desc = desc
        self.vec = vec
        self.func = func

    def prep(self, K, fn):
        (meta, xs) = kfset.read(fn)
        fK = meta['K']
        if fK < K:
            raise MismatchedK(K, fK)
        if self.vec:
            S = 2*(fK - K)
            v = array.array('I', [0 for i in xrange(1 << (2*K))])
            for (x, c) in xs:
                y = x >> S
                v[y] += c
            return v
        
        S = 2*(fK - K)
        v = array.array('L', [])
        for (x,_) in xs:
            y = x >> S
            if len(v) == 0 or v[-1] != y:
                v.append(y)
        return v

    def measure(self, lhs, rhs):
        return self.func(lhs, rhs, self.vec)


def addMeasure(name, desc, vec, func):
    m = Measure(name, desc, vec, func)
    measures[name] = m

addMeasure('bray.curtis.quant',
           'Quantative Bray.Curtis distance', True, dist.brayCurtis)
addMeasure('bray.curtis.qual',
           'Qualitative Bray.Curtis distance', False, dist.brayCurtis)
addMeasure('chord.quant',
           'Quantative Chord distance', True, dist.chord)
addMeasure('chord.qual',
           'Qualitative Chord distance', False, dist.chord)
addMeasure('hellinger.quant',
           'Quantative Hellinger distance', True, dist.hellinger)
addMeasure('hellinger.qual',
           'Qualitative Hellinger distance', False, dist.hellinger)
addMeasure('jaccard.ab',
           'Abundance.based Jaccard distance', True, dist.jaccard)
addMeasure('jaccard.qual',
           'Qualitative Jaccard distance', False, dist.jaccard)
addMeasure('jensen.shannon',
           'Jensen.Shannon distance', True, dist.jensenShannon)
addMeasure('kulczynski.quant',
           'Quantative Kulczynski distance', True, dist.kulczynski)
addMeasure('kulczynski.qual',
           'Qualitative Kulczynski distance', False, dist.kulczynski)
addMeasure('ochiai.ab',
           'Abundance.based Ochiai distance', True, dist.ochiai)
addMeasure('ochiai.qual',
           'Qualitative Ochiai distance', False, dist.ochiai)
addMeasure('sorensen.ab',
           'Abundance.based Sorensen distance', True, dist.sorensen)
addMeasure('sorensen.qual',
           'Qualitative Sorensen distance', False, dist.sorensen)
addMeasure('whittaker.quant',
           'Quantative Whittaker distance', True, dist.whittaker)
addMeasure('whittaker.qual',
           'Qualitative Whittaker distance', False, dist.whittaker)

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    if "list" in opts['-M']:
        ms = measures.keys()
        ms.sort()
        msg = []
        for m in ms:
            msg.append(m + '\t' + measures[m].desc)
        print '\n'.join(msg)
        return

    allMs = measures.keys()
    allMs.sort()

    seen = set([])
    bad = False
    for mo in opts['-M']:
        found = False
        for m in allMs:
            if fnmatch.fnmatch(m, mo):
                seen.add(m)
                found = True
        if not found:
            print >> sys.stderr, 'warning: measure \'%s\' not found. Use -M list to see all measures.' % (mo,)
            bad = True
    ms = list(seen)
    ms.sort()

    if len(ms) == 0 or bad:
        return

    K = int(opts['<k>'])

    N = len(opts['<input>'])
    fns = opts['<input>']

    hdr = ['lhs.name', 'rhs.name']
    fmt = ['%s', '%s']
    vecNeeded = None
    setNeeded = None
    for m in ms:
        hdr.append(m)
        fmt.append('%g')
        if measures[m].vec:
            vecNeeded = m
        else:
            setNeeded = m
    fmt = '\t'.join(fmt)

    print '\t'.join(hdr)
    for i in xrange(N):
        lhsVec = None
        if vecNeeded is not None:
            lhsVec = measures[vecNeeded].prep(K, fns[i])
        lhsSet = None
        if setNeeded is not None:
            lhsSet = measures[setNeeded].prep(K, fns[i])

        for j in xrange(i + 1, N):
            rhsVec = None
            if vecNeeded is not None:
                rhsVec = measures[vecNeeded].prep(K, fns[j])
            rhsSet = None
            if setNeeded is not None:
                rhsSet = measures[setNeeded].prep(K, fns[j])

            vs = [fns[i], fns[j]]
            for m in ms:
                if measures[m].vec:
                    d = measures[m].measure(lhsVec, rhsVec)
                else:
                    d = measures[m].measure(lhsSet, rhsSet)
                vs.append(d)
            print fmt % tuple(vs)

if __name__ == '__main__':
    main(sys.argv[1:])
