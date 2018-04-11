import sys
from pykmer.file import openFile, readFasta, readFastq

root = sys.argv[1]

chs = set([])
bed = {}
with openFile(sys.argv[2]) as f:
    for l in f:
        t = l.split()
        ch = t[0]
        st = int(t[1])
        en = int(t[2])
        nm = t[3]
        bed[nm] = (ch, st, en)
        chs.add(ch)

chs = dict([(ch, list(readFasta(open(root + '/' + ch + '.fa')))[0][1]) for ch in sorted(chs)])

for fn in sys.argv[3:]:
    fst = fn.endswith('_1.fastq')
    with openFile(fn) as f:
        for rd in readFastq(f):
            t = rd[0].split()
            nm = t[1]
            pos = int(t[2])
            fragLen = int(t[3])
            strand = t[4]
            rdSeq = rd[1]
            rdLen = len(rdSeq)
            ch = bed[nm][0]
            frag = chs[ch][pos-1:pos+fragLen]
            print ch, strand, rdSeq, frag
