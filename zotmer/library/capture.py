from zipfile import ZipFile, ZIP_DEFLATED
from tqdm import tqdm
from pykmer.basics import kmersList

class capture(object):
    def __init__(self, K, **kwargs):
        self.K = K
        self.names = []
        self.nameIdx = {}
        self.baits = {}
        self.reads = kwargs.get('reads', False)
        self.kmers = kwargs.get('kmers', True)
        self.verbose = kwargs.get('verbose', False)
        self.capReads = []
        self.capKmers = []

    def addBait(self, nm, seq, bothStrands=True):
        n = len(self.names)
        self.names.append(nm)
        self.nameIdx[nm] = n
        for x in kmersList(self.K, seq, bothStrands):
            if x not in self.baits:
                self.baits[x] = [n]
            elif self.baits[x][-1]  != n:
                self.baits[x].append(n)
        self.capReads.append([[], []])
        self.capKmers.append({})
        return n

    def addRead(self, rd):
        assert self.reads

        xs = kmersList(self.K, rd[1], True)
        ns = []
        for x in xs:
            if x not in self.baits:
                continue
            ns += self.baits[x]

        if len(ns) == 0:
            return

        ns = set(ns)

        for n in ns:
            self.capReads[n][0] += rd

    def addKmers(self, xs):
        assert self.kmers

        ns = []
        for x in xs:
            if x not in self.baits:
                continue
            ns += self.baits[x]

        if len(ns) == 0:
            return

        for n in ns:
            ys = self.capKmers[n]
            for x in xs:
                if x not in ys:
                    ys[x] = 0
                ys[x] += 1

    def addReadAndKmers(self, rd):
        xs = kmersList(self.K, rd[1], True)
        ns = []
        for x in xs:
            if x not in self.baits:
                continue
            ns += self.baits[x]

        if len(ns) == 0:
            return

        ns = set(ns)

        if self.reads:
            for n in ns:
                self.capReads[n][0] += rd

        if self.kmers:
            for n in ns:
                ys = self.capKmers[n]
                for x in xs:
                    if x not in ys:
                        ys[x] = 0
                    ys[x] += 1

    def addReadPairAndKmers(self, lhs, rhs):
        xs = kmersList(self.K, lhs[1], True)
        ys = kmersList(self.K, rhs[1], True)
        ns = []
        for x in xs:
            if x not in self.baits:
                continue
            ns += self.baits[x]
        for y in ys:
            if y not in self.baits:
                continue
            ns += self.baits[y]

        if len(ns) == 0:
            return

        ns = set(ns)

        if self.reads:
            for n in ns:
                self.capReads[n][0] += lhs
                self.capReads[n][1] += rhs

        if self.kmers:
            for n in ns:
                zs = self.capKmers[n]
                for x in xs:
                    if x not in zs:
                        zs[x] = 0
                    zs[x] += 1
                for y in ys:
                    if y not in zs:
                        zs[y] = 0
                    zs[y] += 1

    def saveReads(self, zipName):
        prog = None
        if self.verbose:
            prog = tqdm(total=len(self.names))
        with ZipFile(zipName, 'w', ZIP_DEFLATED) as z:
            for n in range(len(self.names)):
                if prog is not None:
                    prog.update()
                if len(self.capReads[n][0]) == 0:
                    continue
                if len(self.capReads[n][0]) != len(self.capReads[n][1]):
                    assert len(self.capReads[n][1]) == 0
                    nm = self.names[n]
                    z.writestr(nm + '/reads.fastq', '\n'.join(self.capReads[n][0]))
                else:
                    nm = self.names[n]
                    z.writestr(nm + '/reads_1.fastq', '\n'.join(self.capReads[n][0]) + '\n')
                    z.writestr(nm + '/reads_2.fastq', '\n'.join(self.capReads[n][1]) + '\n')
        if prog is not None:
            prog.close()
