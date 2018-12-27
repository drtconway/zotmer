import random

from zotmer.library.basics import render
from zotmer.library.file import autoremove, tmpfile
from zotmer.library.container.casket import casket
from zotmer.library.files import readVectorList, writeVector, readKmers, readKmersList, writeKmersList

def test_rwVector():
    N = 65536

    random.seed(17)
    xs = [int(0.5+10*random.expovariate(0.5)) for i in xrange(N)]

    with autoremove():
        t = tmpfile()

        with casket(t, 'w') as z:
            writeVector(z, xs, 'quux')

        with casket(t, 'r') as z:
            ys = readVectorList(z, 'quux')

        assert len(xs) == len(ys)
        for i in xrange(len(xs)):
            assert xs[i] == ys[i]

def test_kmersList():
    K = 25
    M = (1 << (2*K)) - 1
    N = 65536

    random.seed(17)
    xs = [random.randint(0, M) for i in xrange(N)]
    xs.sort()

    with autoremove():
        t = tmpfile()

        with casket(t, 'w') as z:
            zs = [x for x in xs]
            writeKmersList(z, zs)

        with casket(t, 'r') as z:
            ys = list(readKmers(z))

        assert len(xs) == len(ys)
        for i in xrange(len(xs)):
            assert xs[i] == ys[i], '%d\t%s\t%s' % (i, render(K, xs[i]), render(K, ys[i]))
