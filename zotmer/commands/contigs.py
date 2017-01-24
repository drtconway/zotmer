"""
Usage:
    zot contigs [-l LEN] <input>...

Perform a simple de Bruijn graph assemly of the input k-mers.  Only
non-branching paths in the de Bruijn graph are reported: no attempt
is made to resolve branching at all.

Options:
    -l LEN          minimum length of contig to report (default 2K)
"""

from pykmer.basics import rc, render
from pykmer.container import container
from pykmer.container.std import readKmers
from pykmer.sparse import sparse
from pykmer.bitvec import bitvec

import docopt
import array
import sys

def succ(K, xs, x):
    m = (1 << (2*K)) - 1
    y0 = (x << 2) & m
    y1 = y0 + 4
    (r0, r1) = xs.rank2(y0, y1)
    return range(r0, r1)

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    L0 = None
    if opts['-l']:
        L0 = int(opts['-l'])

    for inp in opts['<input>']:
        with container(inp, 'r') as z:
            K = z.meta['K']
            L = L0
            if L is None:
                L = 2*K

            xs = array.array('L', readKmers(z))
            S = sparse(2*K, xs)

            seen = bitvec(S.count())
            for i in xrange(S.count()):
                if seen[i]:
                    continue

                x = S.select(i)
                xb = rc(K, x)
                xp = succ(K, S, xb)
                if xp == 1:
                    # x isn't the start of a contig
                    continue

                pth = [x]
                seen[i] = 1
                xn = succ(K, S, x)
                while len(xn) == 1:
                    if seen[xn[0]] == 1:
                        break
                    x = S.select(xn[0])
                    pth.append(x)
                    seen[xn[0]] = 1
                    xb = rc(K, x)
                    j = S.rank(xb)
                    seen[j] = 1
                    xn = succ(K, S, x)

                if len(pth)+K-1 < L:
                    continue

                s = [render(K, pth[0])]
                for j in xrange(1, len(pth)):
                    s.append("ACGT"[pth[j]&3])

                print '>contig_%d\n%s' % (i, ''.join(s))

if __name__ == '__main__':
    main(sys.argv[1:])
