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

class InformationGainSplitter(object):
    def __init__(self, X, Y, T):
        self.X = X
        self.Y = Y
        self.T = T

    def score(self, S):
        return self.entropy(S)

    def labelProbabilities(self, S):
        P = {}
        for i in S:
            y = self.Y[i]
            if y not in P:
                P[y] = 0
            P[y] += 1
        t = len(S)
        return dict([(y,float(n)/float(t)) for (y,n) in P.items()])

    def scoreSplits(self, S, j, E0):
        if self.T[j] == 'cat':
            return self.scoreCategorical(S, j, E0)
        if self.T[j] == 'ord':
            return self.scoreOrdinal(S, j, E0)
        return self.scoreContinuous(S, j, E0)

    def scoreCategorical(self, S, j, E0):
        if E0 is None:
            E0 = self.entropy(S)

        W = set([self.X[i][j] for i in S])

        if len(W) == 1:
            return

        for x in W:
            lhsS = []
            rhsS = []
            for i in S:
                if self.X[i][j] == x:
                    lhsS.append(i)
                else:
                    rhsS.append(i)
            assert len(lhsS) > 0
            assert len(rhsS) > 0
            lhsP = float(len(lhsS)) / float(len(S))
            lhsE = self.entropy(lhsS)
            rhsP = float(len(rhsS)) / float(len(S))
            rhsE = self.entropy(rhsS)
            ig = E0 - (lhsP*lhsE + rhsP*rhsE)
            yield (ig, j, 'eq', x, lhsE, rhsE, lhsS, rhsS)

    def scoreOrdinal(self, S, j):
        return self.scoreContinuous(S, j, E0)

    def scoreContinuous(self, S, j, E0):

        W = sorted(set([self.X[i][j] for i in S]))

        for x in W[:-1]:
            lhsS = []
            rhsS = []
            for i in S:
                if self.X[i][j] <= x:
                    lhsS.append(i)
                else:
                    rhsS.append(i)
            assert len(lhsS) > 0
            assert len(rhsS) > 0
            lhsP = float(len(lhsS)) / float(len(S))
            lhsE = self.entropy(lhsS)
            rhsP = float(len(rhsS)) / float(len(S))
            rhsE = self.entropy(rhsS)
            ig = E0 - (lhsP*lhsE + rhsP*rhsE)
            yield (ig, j, 'le', x, lhsE, rhsE, lhsS, rhsS)

    def entropy(self, S):
        P = {}
        for i in S:
            y = self.Y[i]
            if y not in P:
                P[y] = 0
            P[y] += 1
        return entropy(P)

class decisionTree(object):
    def __init__(self, tree):
        """
        Create a decision tree object

        A decision tree has internal nodes which are tuples
            (featureIndex, relation, value, lhs, rhs)
        where
            featureIndex    the feature to branch on: [0, M)
            relation        the compariston to make for a left branch: eq, le
            value           the value to compare against
            lhs,rhs         the left-hand and right-hand children

        A leaf node is a dictionary with labels as keys and
        the associated probability of that label as values.
        """
        self.root = tree

    @staticmethod
    def make(S, M, X, Y, T):
        """
        Make a new decision tree.

        S is the set of sample numbers from [0, N) to use from X
        M is the number of features
        X is the N*M matrix of training data
        Y is the N length vector of labels
        T the M length vector of feature types
        """
        splitter = InformationGainSplitter(X, Y, T)
        E0 = splitter.score(S)
        root = decisionTree.makeRec(S, M, E0, splitter)
        return decisionTree(root)

    @staticmethod
    def makeRec(S, M, E0, splitter):
        if E0 < 1e-3:
            # If the entropy is very low, then make a leaf node
            return splitter.labelProbabilities(S)

        # We're going to sample a set of root-M features to
        # choose from for creating this split.
        #
        n = int(1+math.sqrt(M))
        if M < 1000:
            # For smallish numbers of features just
            # shuffle a list and pull off the first n
            #
            J = range(M)
            random.shuffle(J)
            J = set(J[:n])
        else:
            # When the universe of features is large,
            # pick random ones.
            #
            # TODO: avoid degenerate features.
            #
            J = set([])
            while len(J) < n:
                j = random.randint(0, M-1)
                J.add(j)

        splits = []
        for j in J:
            for s in splitter.scoreSplits(S, j, E0):
                splits.append(s)

        if len(splits) == 0:
            # If there are no splits, then make a leaf
            return splitter.labelProbabilities(S)

        # Among the possible splits, make a stochastic choice
        # based on the relative score for each split.
        #
        chosenS = None
        totalScore = sum([s[0] for s in splits])
        u = random.random()
        for s in splits:
            p = s[0]/totalScore
            if u < p:
                chosenS = s
                break
            u -= p
        assert chosenS is not None

        (sc, j, rel, x, lhsE, rhsE, lhsS, rhsS) = chosenS

        lhs = decisionTree.makeRec(lhsS, M, lhsE, splitter)
        rhs = decisionTree.makeRec(rhsS, M, rhsE, splitter)
        return (j, rel, x, lhs, rhs)

    def decide(self, xs):
        n = self.root
        while isinstance(n, tuple):
            (i, rel, v, lhs, rhs) = n
            if rel == 'le':
                if xs[i] <= v:
                    n = lhs
                else:
                    n = rhs
            else:
                assert rel == 'eq'
                if xs[i] == v:
                    n = lhs
                else:
                    n = rhs
        return n

class ranfor(object):
    def __init__(self):
        """Initialize an empty random forest"""
        self.T = []

    def make(self, N, M, B, X, Y, T):
        """
        Train a random forest.

        N is the number of training samples (rows) in X
        M is the number of features (columns) in X
        B is the number of trees to create
        X is the training data matrix
        Y is the vector of (N) labels for the training data
        T is the vector of types of the features: cat[egorical], ord[inal], con[tinuous]
        """
        n = int(math.ceil(math.sqrt(N)))
        for b in xrange(B):
            S = [random.randint(0, N-1) for i in xrange(n)]
            J = range(M)
            t = decisionTree.make(S, J, X, Y, T)
            self.T.append(t)

    def classify(self, v):
        """
        Take a feature vector and return a dict with the labels as keys, and probability as values.
        """
        votes = {}
        for t in self.T:
            cps = t.decide(v)
            for (c,p) in cps.items():
                if c not in votes:
                    votes[c] = 0.0
                votes[c] += p
        pt = sum(votes.values())
        return dict([(c, p/pt) for (c,p) in votes.items()])
