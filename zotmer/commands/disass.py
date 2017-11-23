"""
zot disass - "disassemble" contigs to produce some summary statistics

Usage:
    zot disass [options] <input>...

Options:
    -k K            value of k to use [default: 25]
    -c C            cutoff for distinguishing low and high frequency k-mers [default: 5]
    -p Prob         subsample k-mers with probability Prob [default: 1.0]
    -q Q            number of quantiles to compute [default: 10]
    -s              perform a single stranded analysis
    -S Seed         seed for subsampling [default: 17]
    -v              produce verbose output

"""

import math
import sys

import docopt
import yaml

from pykmer.basics import kmersList, sub
from pykmer.file import openFile, readFasta

def summarize(xs, cut, Q):
    res = {}
    h = {}
    cs = xs.values()
    cs.sort()
    for c in cs:
        h[c] = 1 + h.get(c, 0)

    res['histogram'] = sorted(h.items())
    res['mean'] = float(sum(cs))/float(max(1, len(cs)))
    res['median'] = 0
    if len(cs) > 0:
        if (len(cs) & 1) == 0:
            res['median'] = cs[len(cs)//2]
        else:
            m = len(cs)//2
            res['median'] = (cs[m] + cs[m+1]) / 2.0

    res['low-count'] = 0
    res['high-count'] = 0

    t = float(sum(h.values()))
    q0 = t/Q
    q = q0
    quant = []
    cum = 0
    for (c,f) in h.items():
        if c < cut:
            res['low-count'] += f
        else:
            res['high-count'] += f
        while cum + f > q:
            quant.append(c)
            q += q0
        cum += f
    res['quantiles'] = quant

    return res

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    K = int(opts['-k'])
    C = int(opts['-c'])
    Q = int(opts['-q'])
    S = int(opts['-S'])
    P = float(opts['-p'])

    verbose = opts['-v']

    both = True
    if opts['-s']:
        both = False

    res = []
    for fn in opts['<input>']:
        fres = {}
        fres['file'] = fn
        fres['contigs'] = []
        glob = {}
        ncontig = 0
        with openFile(fn) as f:
            for (nm, seq) in readFasta(f):
                ncontig += 1
                scaff = {}
                for x in kmersList(K, seq, both):
                    if sub(S, P, x):
                        scaff[x] = 1 + scaff.get(x, 0)
                summary = summarize(scaff, C, Q)
                summary['name'] = nm
                fres['contigs'].append(summary)
                for (x,c) in scaff.items():
                    glob[x] = c + glob.get(x, 0)
        fres['global'] = summarize(glob, C, Q)
        res.append(fres)

    yaml.safe_dump(res, sys.stdout)

if __name__ == '__main__':
    main(sys.argv[1:])
