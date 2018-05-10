from pykmer.basics import kmersList, kmersWithPosList, kmersWithPosLists

class kmerIndex(object):
    def __init__(self, K):
        self.K = K
        self.seqs = {}
        self.idx = {}

    def addSeq(self, nm, seq):
        assert nm not in self.seqs
        self.seqs[nm] = seq
        for x in kmersList(self.K, seq, False):
            if x not in self.idx:
                self.idx[x] = [nm]
            elif self.idx[x][-1] != nm:
                self.idx[x].append(nm)

    def getSeq(self, nm):
        return self.seqs[nm]

    def __getitem__(self, x):
        if x not in self.idx:
            return []
        return self.idx[x]

class kmerPosIndex(object):
    def __init__(self, K):
        self.K = K
        self.seqs = {}
        self.idx = {}
        self.tbl = {}

    def addSeq(self, nm, seq):
        assert nm not in self.seqs
        self.seqs[nm] = seq
        for (x,p) in kmersWithPosList(self.K, seq, False):
            p -= 1
            if x not in self.idx:
                self.idx[x] = []
            self.idx[x].append((nm,p))
            if nm not in self.tbl:
                self.tbl[nm] = {}
            self.tbl[nm][p] = x

    def getSeq(self, nm):
        return self.seqs[nm]

    def __getitem__(self, x):
        if x not in self.idx:
            return []
        return self.idx[x]

    def atPos(self, nm, p):
        if nm not in self.tbl:
            return None
        if p not in self.tbl[nm]:
            return None
        return self.tbl[nm][p]

    def simpleHits(self, xps):
        loc = {}
        for (x,p) in xps:
            if x not in self.idx:
                continue
            for (z,q) in self.idx[x]:
                if z not in loc:
                    loc[z] = {}
                if q not in loc[z]:
                    loc[z][q] = 0
                loc[z][q] += 1
        return loc

    def hits(self, xps):
        loc = {}
        for (x,p) in xps:
            p -= 1
            if x not in self.idx:
                continue
            for (z,q) in self.idx[x]:
                r = q - p
                if z not in loc:
                    loc[z] = {}
                if r not in loc[z]:
                    loc[z][r] = 0
                loc[z][r] += 1
        return loc

    def accumAll(self, loc, xps, acc):
        for z in loc.keys():
            if z not in acc:
                acc[z] = {}
            for r in loc[z].keys():
                for (x,p) in xps:
                    p -= 1
                    q = r + p
                    if q not in acc[z]:
                        acc[z][q] = {}
                    if x not in acc[z][q]:
                        acc[z][q][x] = 0
                    acc[z][q][x] += 1

    def readHits(self, seq):
        (fwd, rev) = kmersWithPosLists(self.K, seq)
        fwdLoc = self.hits(fwd)
        revLoc = self.hits(rev)
        return (fwd, fwdLoc, rev, revLoc)

    def pairHits(self, lhs, rhs, acc):
        (lhsFwd, lhsFwdHits, lhsRev, lhsRevHits) = self.readHits(lhs)
        self.accumAll(lhsFwdHits, lhsFwd, acc)
        self.accumAll(lhsRevHits, lhsRev, acc)

        (rhsFwd, rhsFwdHits, rhsRev, rhsRevHits) = self.readHits(rhs)
        self.accumAll(rhsFwdHits, rhsFwd, acc)
        self.accumAll(rhsRevHits, rhsRev, acc)

