import os
import struct

from pykmer.container.casket import casket
from pykmer.file import tmpfile
import pykmer.codec64 as codec64

class fileEncoder(codec64.encoder):
    def __init__(self, f):
        codec64.encoder.__init__(self)
        self.f = f

    def write(self, w):
        self.f.write(struct.pack('Q', w))

    def __enter__(self):
        return self

    def __exit__(self, t, f, v):
        self.end()
        if t is not None:
            return False
        return True

class kmerWriter:
    def __init__(self, f):
        self.f = fileEncoder(f)
        self.p = 0

    def append(self, x):
        d = x - self.p
        self.f.append(d)
        self.p = x

    def __enter__(self):
        return self

    def __exit__(self, t, f, v):
        return self.f.__exit__(t, f, v)

class countsWriter:
    def __init__(self, f):
        self.f = fileEncoder(f)

    def append(self, x):
        self.f.append(x)

    def __enter__(self):
        return self

    def __exit__(self, t, f, v):
        return self.f.__exit__(t, f, v)

def readWords(f):
    B = 65536
    s = f.read(8 * B)
    while len(s) > 0:
        assert (len(s) & 7) == 0
        fmt = '%dQ' % (len(s) / 8,)
        ws = struct.unpack(fmt, s)
        for w in ws:
            yield w
        s = f.read(8 * B)

def writeWords(f, ws):
    B = 65536
    vs = []
    n = 0
    for w in ws:
        vs.append(w)
        n += 1

        if len(vs) == B:
            fmt = '%dQ' % (len(vs),)
            s = struct.pack(fmt, *vs)
            f.write(s)
            vs = []
    if len(vs) > 0:
        fmt = '%dQ' % (len(vs),)
        s = struct.pack(fmt, *vs)
        f.write(s)
        vs = []
    return n

def delta(xs):
    p = 0
    for x in xs:
        d = x - p
        yield d
        p = x

def undelta(ds):
    x = 0
    for d in ds:
        x += d
        yield x

def writeVector(z, xs, nm):
    with z.add_stream(nm) as f:
        ws = codec64.encode(xs)
        writeWords(f, ws)

def readVector(z, nm):
    print nm
    f = z.open(nm)
    ws = readWords(f)
    for x in codec64.decode(ws):
        yield x

def writeDeltas(z, xs, nm):
    return writeVector(z, delta(xs), nm)

def readDeltas(z, nm):
    return undelta(readVector(z, nm))

def writeKmers(z, xs, nm = 'kmers'):
    return writeDeltas(z, xs, nm)

def readKmers(z, nm = 'counts'):
    return readDeltas(z, nm)

def writeCounts(z, xs, nm = 'counts'):
    return writeVector(z, xs, nm)

def readCounts(z, nm = 'counts'):
    return readVector(z, nm)

def demuxKmersAndCounts(xs, f):
    cs = fileEncoder(f)
    for (x,c) in xs:
        yield x
        cs.append(c)
    cs.end()

def muxKmersAndCounts(xs, cs):
    moreXs = True
    try:
        x = xs.next()
    except StopIteration:
        moreXs = False
    moreCs = True
    try:
        c = cs.next()
    except StopIteration:
        moreCs = False
    assert moreXs == moreCs
    while moreXs:
        yield (x,c)
        try:
            x = xs.next()
        except StopIteration:
            moreXs = False
        try:
            c = cs.next()
        except StopIteration:
            moreCs = False
        assert moreXs == moreCs

def writeKmersAndCounts(z, xs, nm = None):
    if nm is None:
        xNm = 'kmers'
        cNm = 'counts'
    else:
        xNm = nm + '-kmers'
        cNm = nm + '-counts'
    t = tmpfile()
    with open(t, 'w') as f:
        ys = demuxKmersAndCounts(xs, f)
        writeKmers(z, ys, xNm)
    z.add_file(cNm, t)
    os.remove(t)

def writeKmersAndCounts2(z, xs, cs, nm = None):
    if nm is None:
        xNm = 'kmers'
        cNm = 'counts'
    else:
        xNm = nm + '-kmers'
        cNm = nm + '-counts'
    writeKmers(z, xs, xNm)
    writeCounts(z, cs, cNm)

def readKmersAndCounts(z, nm = None):
    if nm is None:
        xNm = 'kmers'
        cNm = 'counts'
    else:
        xNm = nm + '-kmers'
        cNm = nm + '-counts'

    return muxKmersAndCounts(readKmers(z, xNm), readCounts(z, cNm))
