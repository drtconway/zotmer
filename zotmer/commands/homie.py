"""
zot homie - scan for homopolymer variants.

Usage:
    zot homie -X [options] <base> <min-length> <max-length> <input>...
    zot homie -M [options] <input>...
    zot homie -S [options] <model> <input>...
    zot homie [options] <index> <input>...

Options:
    -X              scan reference fasta files for homopolymers
    -M              merge counts and construct a reference model for homopolymers
    -S              score samples against a model.
    -B BEDFILE      [indexing] only include homopolymers in the regions from the given BED file.
    -D D            maximum hamming distance for anchor k-mers [default: 1]
    -I FILENAME     [merging] read input files from the named file, 1 per line.
    -k K            k-mer size to use [default: 25]
    -m MODELFILE    score samples against the model produced by merging
    -o FILENMANE    write the output to the named file rather than stdout
    -s SAMPLE       give a sample name to include in the output
    -T              [merging] print output in a tabular form
    -v              produce verbose output
    -W WIDTH        [indexing] width of context to be used for read pulldown [default: 100]

<base> is a 16-letter FASTA description of which nucleotides to scan for.
<min-length> and <max-length> specify the minimum and maximum homopolymer lengths to scan for.
<index> is a yaml file with the output of a 'zot homie -X' invocation.
<model> is a yaml file with the output of a 'zot homie -M' invocation.

"""

import math
import re
import sys

import docopt
import yaml

from zotmer.library.basics import ham, kmer, kmers, kmersList, kmersWithPosList, murmer, rc, render, unFasta
from zotmer.library.capture import capture
from zotmer.library.file import openFile, readFasta
from zotmer.library.misc import intervalMap
from zotmer.library.reads import reads

verbose = False

def fa16ToList(B):
    x = unFasta(B)
    r = []
    for i in range(4):
        if x & (1 << i):
            r.append('ACGT'[i])
    return r

def readBed(fn):
    r = {}
    with open(fn) as f:
        for l in f:
            t = l.split()
            ch = t[0]
            st = int(t[1])
            en = int(t[2])
            nm = t[3]
            if ch not in r:
                r[ch] = intervalMap()
            r[ch].add((st, en), nm)
    return r

def scan(bs, lo, hi, W, bed, fns):
    global verbose

    pats = [(b, '(.{%d}[^%s])(%s{%d,%d})([^%s].{%d})' % (W-1, b, b, lo, hi, b, W-1)) for b in bs]
    ps = [(b,re.compile(pat, re.I)) for (b,pat) in pats]

    for fn in fns:
        if verbose:
            print >> sys.stderr, 'scanning', fn
        with openFile(fn) as f:
            for (nm, seq) in readFasta(f):
                for (b,p) in ps:
                    if verbose:
                        print >> sys.stderr, '    scanning for', b
                    j = 0
                    for m in p.finditer(seq):
                        (lhs, hom, rhs) =  m.groups()
                        n = len(hom)
                        s = m.start(2)
                        e = m.end(2)
                        if nm not in bed:
                            continue
                        g = None
                        if bed is not None:
                            g = bed[nm][s]
                            if g is None:
                                g = bed[nm][e]
                            if g is None:
                                continue
                        j += 1
                        r = {}
                        r['lhsFlank'] = lhs.upper()
                        r['length'] = n
                        r['rhsFlank'] = rhs.upper()
                        r['chr'] = nm
                        r['start'] = s
                        r['end'] = e
                        r['base'] = b
                        if g is not None:
                            r['gene'] = g
                        yield r
                    if verbose:
                        print >> sys.stderr, '    %s loci found.' % (j, )

def findAnchors(X, x0, D):
    for x in X.iterkeys():
        d = ham(x, x0)
        if d <= D:
            yield (x, d)

def spliceKmers(K, x, y, b, j):
    seq = render(K, x) + (render(1, b) * j) + render(K, y)
    return kmersList(K, seq)

def dist(cs):
    if type(cs) == list:
        t = sum(cs)
        v = 1.0/float(t)
        return [c*v for c in cs]
    else:
        t = sum(cs.values())
        v = 1.0/float(t)
        return dict([(j,c*v) for (j,c) in cs.items()])

def kld(ps, qs, N):
    assert type(ps) == type(qs)
    if type(ps) == list:
        d = 0.0
        for (p,q) in zip(ps, qs):
            if p == 0.0:
                continue
            assert q > 0.0
            d += p * math.log(p/q)
        return d
    else:
        q0 = 1.0/(N+1.0)
        d = 0.0
        for (j,p) in ps.items():
            assert p > 0.0
            q = qs.get(j, q0)
            d += p * math.log(p/q)
        return d

def main(argv):
    opts = docopt.docopt(__doc__, argv)

    global verbose
    verbose = opts['-v']

    if opts['-X']:
        W = int(opts['-W'])

        B = opts['<base>']
        lo = int(opts['<min-length>'])
        hi = int(opts['<max-length>'])
        bed = None
        if opts['-B']:
            bed = readBed(opts['-B'])

        out = sys.stdout
        if opts['-o']:
            out = openFile(opts['-o'], 'w')

        bs = fa16ToList(B)
        if len(bs) == 0:
            return

        yaml.safe_dump_all(scan(bs, lo, hi, W, bed, opts['<input>']), out, default_flow_style=False)

        return

    if opts['-M']:
        inputs = opts['<input>']
        if opts['-I']:
            with openFile(opts['-I']) as f:
                for l in f:
                    nm = l.strip()
                    inputs.append(nm)

        I = {}
        R = {}
        J = 0

        for fn in inputs:
            with open(fn) as f:
                if verbose:
                    print >> sys.stderr, 'processing', fn
                for itm in yaml.load_all(f):
                    loc = itm['locus']

                    I[loc] = dict([(k,itm[k]) for k in itm.keys() if k != 'counts' and k != 'sample'])

                    #j = sum(itm['counts'])
                    #if j < 10:
                    #    continue 

                    if loc not in R:
                        R[loc] = itm['counts']
                        J = max(J, max(itm['counts'].keys()))
                    else:
                        for (j,c) in itm['counts'].items():
                            if j in R[loc]:
                                R[loc][j] += c
                            else:
                                R[loc][j] = c
                            J = max(J, j)

        for loc in R.keys():
            cs = R[loc]
            ps = dist(cs)
            I[loc]['counts'] = cs
            I[loc]['dist'] = ps

        for loc in I.keys():
            if 'dist' not in I[loc]:
                del I[loc]

        out = sys.stdout
        if opts['-o']:
            out = openFile(opts['-o'], 'w')

        if opts['-T']:
            for loc in sorted(R.keys()):
                ps = [R[loc].get(j, 0.0) for j in range(J+1)]
                ss = ['%g' % (p, ) for p in ps]
                print '\t'.join([loc, I[loc]['base'], str(I[loc]['length'])] + ss)
            return

        yaml.safe_dump(I, out, default_flow_style=False)

        return

    if opts['-S']:
        with openFile(opts['<model>']) as f:
            model = yaml.load(f)

        out = sys.stdout
        if opts['-o']:
            out = openFile(opts['-o'], 'w')

        for fn in opts['<input>']:
            with openFile(fn) as f:
                for itm in yaml.load_all(f):
                    l = itm['locus']
                    if l not in model:
                        print >> sys.stderr, 'warning: file %s, locus %s was not in the model.' % (fn, l)
                        continue
                    cs = itm['counts']
                    ps = dist(cs)
                    if 'dist' not in model[l]:
                        print >> sys.stderr, model[l]
                    d = kld(ps, model[l]['dist'], sum(model[l]['counts'].values()))
                    nm = fn
                    if 'sample' in itm:
                        nm = itm['sample']
                    g = 'unknown'
                    if 'gene' in itm:
                        g = itm['gene']

                    fmts = ['%s', '%s', '%s', '%g']
                    vals = [nm, l, g, d]

                    print >> out, '\t'.join(fmts) % tuple(vals)
        return

    K = int(opts['-k'])
    D = int(opts['-D'])

    items = []
    locs = []
    cap = capture(K, reads=False, kmers=True)
    with openFile(opts['<index>']) as f:
        for itm in yaml.load_all(f):
            items.append(itm)
            lhs = itm['lhsFlank']
            rhs = itm['rhsFlank']
            loc = '%s:%d-%d' % (itm['chr'], itm['start'], itm['end'])
            locs.append(loc)
            cap.addBait(loc, lhs + 'N' + rhs)
    N = len(items)

    sample = None
    if opts['-s']:
        sample = opts['-s']

    for itm in reads(opts['<input>'], K=K, paired=True, reads=True, kmers=False, both=True, verbose=verbose):
        cap.addReadPairAndKmers(itm.reads[0], itm.reads[1])

    out = sys.stdout
    if opts['-o']:
        out = openFile(opts['-o'], 'w')

    def doit():
        for n in range(N):
            itm = items[n]
            X = cap.capKmers[n]

            lhs = itm['lhsFlank']
            rhs = itm['rhsFlank']
            m = itm['length']
            b = kmer(itm['base'])
            g = itm['gene']

            W = len(lhs)
            lhsAnc = kmer(lhs[-K:])
            lhsAncs = list(findAnchors(X, lhsAnc, D))
            rhsAnc = kmer(rhs[:K])
            rhsAncs = list(findAnchors(X, rhsAnc, D))

            if len(lhsAncs)*len(rhsAncs) > 500:
                print >> sys.stderr, 'warning: for locus %s there wre a large number of anchors (%d, %d)' % (locs[n], len(lhsAncs), len(rhsAncs))

            counts = [0 for j in range(K+1)]
            for (x, xd) in lhsAncs:
                for (y, yd) in rhsAncs:
                    for j in range(K+1):
                        zs = spliceKmers(K, x, y, b, j)
                        cs = [X.get(z, 0) for z in zs]
                        cs.sort()
                        counts[j] = max(counts[j], cs[0])

            if sum(counts) == 0:
                continue

            res = {}
            if sample is not None:
                res['sample'] = sample
            res['locus'] = locs[n]
            res['gene'] = g
            res['base'] = itm['base']
            res['length'] = m
            res['counts'] = dict([(j,c) for (j,c) in zip(range(K+1), counts) if c > 0])
            res['lhsFlank'] = lhs
            res['rhsFlank'] = rhs
            yield res

    yaml.safe_dump_all(doit(), out, default_flow_style=False)

if __name__ == '__main__':
    main(sys.argv[1:])
