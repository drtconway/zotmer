"""
"""

from os.path import basename

from tqdm import tqdm

from zotmer.library.basics import kmersList, kmersLists
from zotmer.library.file import openFile, readFasta, readFastq

compressionSuffixes = ['.gz', '.bz2']

class Reads(object):
    def __init__(self):
        self.reads = []
        self.kmers = []

def stripCompressionSuffix(nm):
    for suff in compressionSuffixes:
        if nm.endswith(suff):
            return nm[:-len(suff)]
    return nm

def isFasta(nm):
    bnm = stripCompressionSuffix(nm)
    if bnm.endswith(".fa"):
        return True
    if bnm.endswith(".fasta"):
        return True
    if bnm.endswith(".fas"):
        return True
    if bnm.endswith(".fna"):
        return True

class reads(object):
    def __init__(self, files, **kwargs):
        self.files = files
        self.K = kwargs.get('K', 25)
        self.paired = kwargs.get('paired', False)
        self.reads = kwargs.get('reads', False)
        self.kmers = kwargs.get('kmers', True)
        self.fwdOnly = kwargs.get('fwdOnly', False)
        self.separate = kwargs.get('separate', False)
        self.both = kwargs.get('both', False)
        self.verbose =  kwargs.get('verbose', False)

        assert not self.kmers or sum([self.fwdOnly, self.separate, self.both]) == 1

        self.N = 1 + (self.paired)
        self.currFilesInd = None
        self.currParsers = None
        self.currReads = None
        self.currKmers = None
        self.readNum = 0

        self.progress = None
        self.M = (1 << 16) - 1

    def __iter__(self):
        return self

    def next(self):
        self.readNum += 1
        if (self.readNum & self.M) == 0 and self.progress is not None:
            self.progress.update(self.M)

        while True:
            if self.currParsers is None:

                if self.currFilesInd is None:
                    self.currFilesInd = 0
                else:
                    self.currFilesInd += self.N

                if self.progress is not None:
                    self.progress.update(self.readNum & self.M)

                if self.currFilesInd + (self.N-1) >= len(self.files):
                    raise StopIteration

                if self.verbose:
                    pfx = ' & '.join([basename(self.files[i]) for i in range(self.currFilesInd, self.currFilesInd+self.N)])
                    self.progress = tqdm(unit=' reads', unit_scale=True)
                    self.progress.set_postfix(reading=pfx, refresh=True)

                self.currParsers = []
                for i in range(self.currFilesInd, self.currFilesInd+self.N):
                    fn = self.files[i]
                    f = openFile(fn)
                    if isFasta(fn):
                        self.currParsers.append(readFasta(f))
                    else:
                        self.currParsers.append(readFastq(f))
           
            self.currReads = []
            try:
                for p in self.currParsers:
                    self.currReads.append(p.next())
            except StopIteration:
                if len(self.currReads) != 0:
                    print >> sys.stderr, 'warning: files had unequal length'
                self.currParsers = None
                if self.progress is not None:
                    self.progress.close()
                    self.progress = None
                continue

            if self.kmers:
                self.currKmers = []
                for rd in self.currReads:
                    if self.fwdOnly:
                        self.currKmers.append(kmersList(self.K, rd[1], False))
                    elif self.both:
                        self.currKmers.append(kmersList(self.K, rd[1], True))
                    else:
                        assert self.separate
                        self.currKmers.append(kmersLists(self.K, rd[1]))

            res = Reads()

            if self.reads:
                res.reads = self.currReads
            if self.kmers:
                res.kmers = self.currKmers
            return res
