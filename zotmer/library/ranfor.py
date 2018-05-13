import math
import random

def entropy(d):
    if len(d) == 0:
        return 0.0
    t = 1.0/sum(d.values())
    e = 0.0
    for (x,n) in d.iteritems():
        p = n*t
        e += -p * math.log(p)
    return e

class decisionTree(object):
    def __init__(self, tree):
        self.root = tree

    @staticmethod
    def make(ns, ms, X, Y):
        root = decisionTree.makeRec(ns, ms, X, Y)
        return decisionTree(root)

    @staticmethod
    def makeRec(ns, ms, X, Y):
        j = random.choice(ms)
        r = decisionTree.split(ns, j, X, Y)
        if len(r) == 1:
            return r
        (x, e, lNs, rNs) = r
        lhs = decisionTree.makeRec(lNs, ms, X, Y)
        rhs = decisionTree.makeRec(rNs, ms, X, Y)
        return (j, x, lhs, rhs)

    @staticmethod
    def split(ns, j, X, Y):
        D = {}
        T = {}
        S = {}
        for i in ns:
            y = Y[i]
            if y not in S:
                S[y] = 0
                T[y] = 0
            T[y] += 1
            x = X[i][j]
            if x not in D:
                D[x] = {}
            if y not in D[x]:
                D[x][y] = 0
            D[x][y] += 1
        
        if len(T) == 1:
            return (T.keys()[0],)

        t = sum(T.values())

        e = 0.0
        for (y,n) in T.items():
            p = float(n) / float(t)
            if p > 0.0:
                e += n * -p*math.log(p)

        E = []
        for x in sorted(D.keys()):
            for (y,n) in D[x].items():
                S[y] += n
            t0 = sum(S.values())
            if t0 == t:
                break
            eL = 0.0
            eR = 0.0
            for y in S.keys():
                p = float(S[y]) / float(t0)
                if p > 0.0:
                    eL += S[y] * -p*math.log(p)
                q = float(T[y] - S[y]) / float(t - t0)
                if q > 0.0:
                    eR += (T[y] - S[y]) * -q*math.log(q)
            eGain = e - (eL + eR)
            E.append((eGain, x, e, eL, eR))
        E.sort()
        (eGain, x, e, eL, eR) = E[-1]
        lNs = []
        rNs = []
        for i in ns:
            if X[i][j] <= x:
                lNs.append(i)
            else:
                rNs.append(i)
        return (x, e, lNs, rNs)

    def decide(self, xs):
        n = self.root
        assert n is not None
        while len(n) > 1:
            (i, v, lhs, rhs) = n
            if xs[i] <= v:
                n = lhs
            else:
                n = rhs
        return n[0]

class ranfor(object):
    def __init__(self):
        self.T = []

    def make(self, N, M, X, Y):
        n = int(math.ceil(math.sqrt(N)))
        m = int(math.ceil(math.sqrt(M)))
        B = 1000
        for b in xrange(B):
            ns = [random.randint(0, N-1) for i in xrange(n)]
            ms = [random.randint(0, M-1) for i in xrange(m)]
            t = decisionTree.make(ns, ms, X, Y)
            self.T.append(t)

    def classify(self, v):
        votes = {}
        for t in self.T:
            j = t.decide(v)
            if j not in votes:
                votes[j] = 0
            votes[j] += 1
        maxJ = []
        maxN = 0
        for (j,n) in votes.iteritems():
            if n < maxN:
                continue
            if n > maxN:
                maxN = n
                maxJ = []
            maxJ.append(j)
        return random.choice(maxJ)

