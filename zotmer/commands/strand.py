"""
zot strand - compute strand bias statistics

Usage:
    zot strand [options] <fastq>...

Options:
    -k K            k-mer size to use [default: 25]
    -r REF          use reference to anchor strands
    -p P            sampling fraction [default: 0.1]
    -s              single ended reads only
    -v              produce verbose output
"""

import random
import sys

import docopt

from zotmer.library.basics import kmersList, murmer, rc, render
from zotmer.library.file import openFile, readFasta, readFastq

def pairs(xs):
    assert len(xs) & 1 == 0
    i = 0
    while i + 1 < len(xs):
        yield (xs[i], xs[i + 1])
        i += 2

def both(xs, ys):
    while True:
        try:
            x = xs.next()
            y = ys.next()
            yield(x, y)
        except StopIteration:
            return

def parseFiles(K, paired, fns, verbose):
    M = (1 << 18) - 1
    rn = 0

    if not paired:
        for fn in fns:
            with openFile(fn) as f:
                rn += 1
                if verbose and (rn & M) == 0:
                    print >> sys.stderr, 'reads processed: %d' % (rn,)
                xs = kmersList(K, fq1[1], False)
                yield xs
        return

    for (fn1, fn2) in pairs(fns):
        with openFile(fn1) as f1, openFile(fn2) as f2:
            for fq1, fq2 in both(readFastq(f1), readFastq(f2)):
                rn += 1
                if verbose and (rn & M) == 0:
                    print >> sys.stderr, 'read pairs processed: %d' % (rn,)
                xs = kmersList(K, fq1[1], False) + [rc(K, x) for x in kmersList(K, fq2[1], False)]
                yield xs

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = int(opts['-k'])
    M = (1 << (2*K)) - 1

    paired = True
    if opts['-s']:
        paired = False

    p = float(opts['-p'])
    T = int(M * p)

    if opts['-r']:
        refs = []
        with openFile(opts['-r']) as f:
            for (nm, seq) in readFasta(f):
                refs += kmersList(K, seq, False)
        refs = set(refs)

        kill = set([])
        for x in refs:
            y = rc(K, x)
            if y in refs:
                kill.add(x)
                kill.add(y)
        print >> sys.stderr, 'removing %d/%d' % (len(kill), len(refs))

        refs -= set(kill)

        fwd = {}
        rev = {}
        for xs in parseFiles(K, paired, opts['<fastq>'], opts['-v']):
            fn = 0
            for x in xs:
                if x in refs:
                    fn += 1

            ys = [rc(K, x) for x in xs]
            rn = 0
            for y in ys:
                if y in refs:
                    rn += 1
            
            if fn + rn == 0:
                continue

            q = float(fn) / float(fn + rn)
            if random.random() < q:
                for x in xs:
                    fwd[x] = 1 + fwd.get(x, 0)
            else:
                for y in ys:
                    rev[y] = 1 + rev.get(y, 0)

        for (x,xc) in fwd.iteritems():
            y = rc(K, x)
            yc = 0
            if y in rev:
                yc = rev[y]
                del rev[y]
            print '%d\t%d' % (xc, yc)

        for (y,yc) in rev.iteritems():
            print '%d\t%d' % (0, yc)

        return

    kx = {}
    for xs in parseFiles(K, paired, opts['<fastq>'], opts['-v']):
        for x in xs:
            if x in kx:
                kx[x] += 1
                continue
            y = rc(K, x)
            z = murmer(min(x, y), 17)
            if (z & M) > T:
                continue
            kx[x] = 1

    for x in kx.keys():
        y = rc(K, x)
        if x > y:
            continue
        xc = kx[x]
        yc = kx.get(y, 0)
        if murmer(x, 17) >= murmer(y, 17):
            (a, b) = (x, y)
            (ac, bc) = (xc, yc)
        else:
            (a, b) = (y, x)
            (ac, bc) = (yc, xc)
        #print '%s\t%d\t%s\t%d' % (render(K, a), ac, render(K, b), bc)
        print '%d\t%d' % (ac, bc)

if __name__ == '__main__':
    main(sys.argv[1:])
