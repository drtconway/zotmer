import array

def matrix(N, M):
    m = []
    for i in xrange(N):
        m.append(array.array('i', [0 for j in xrange(M)]))
    return m

def smithWaterman(a, b):
    W1 = 5
    n = len(a)
    m = len(b)
    H = matrix(n+1, m+1)
    maxH = 0
    maxP = (0,0)
    for i in xrange(1, n+1):
        for j in xrange(1, m+1):
            m0 = -1
            if a[i-1] == b[j-1]:
                m0 = 1
            h0 = H[i-1][j-1] + m0
            h1 = H[i-1][j] - W1
            h2 = H[i][j-1] - W1
            h3 = 0
            h = max([h0, h1, h2, h3])
            H[i][j] = h
            if h > maxH:
                maxH = h
                maxP = (i,j)

    (i, j) = maxP
    resA = []
    resB = []
    while i >= 0 and j >= 0 and H[i][j] > 0:
        h = H[i-1][j-1]
        i1 = i-1
        j1 = j-1
        rA = a[i-1]
        rB = b[j-1]
        if rA != rB:
            rA = rA.lower()
            rB = rB.lower()
        if H[i-1][j] > h:
            h = H[i-1][j]
            i1 = i-1
            j1 = j
            rA = a[i-1]
            rB = '-'
        if H[i][j-1] > h:
            h = H[i][j-1]
            i1 = i
            j1 = j-1
            rA = '-'
            rB = b[j-1]
        i = i1
        j = j1
        resA.append(rA)
        resB.append(rB)

    return (maxH, i, ''.join(resA[::-1]), j, ''.join(resB[::-1]))

def glocalAlignment(a, b):
    g = -5
    n = len(a)
    m = len(b)
    H = matrix(n+1, m+1)
    for j in xrange(1, m+1):
        H[0][j] = H[0][j-1] + g

    for i in xrange(1, n+1):
        for j in xrange(1, m+1):
            m0 = -1
            if a[i-1] == b[j-1]:
                m0 = 1
            h0 = ((H[i-1][j-1] >> 2) + m0, 3)
            h1 = ((H[i-1][j] >> 2) + g, 2)
            h2 = ((H[i][j-1] >> 2) + g, 1)
            (h, d) = max([h0, h1, h2])
            H[i][j] = (h << 2) + d

    (maxS, i) = max([(H[i][m] >> 2, i) for i in xrange(n)])
    j = m

    resA = []
    resB = []
    while i > 0 and j > 0:
        d = H[i][j] & 3
        if d == 3:
            i1 = i-1
            j1 = j-1
            rA = a[i-1]
            rB = b[j-1]
            if rA != rB:
                rA = rA.lower()
                rB = rB.lower()
        elif d == 2:
            i1 = i-1
            j1 = j
            rA = a[i-1]
            rB = '-'
        elif d == 1:
            i1 = i
            j1 = j-1
            rA = '-'
            rB = b[j-1]
        else:
            break
        i = i1
        j = j1
        resA.append(rA)
        resB.append(rB)

    return (maxS, i, ''.join(resA[::-1]), j, ''.join(resB[::-1]))

rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def revComp(s):
    r = []
    for c in s:
        r.append(rcDict.get(c, c))
    return ''.join(r[::-1])

