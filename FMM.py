#FMM.py
import math
import time
import random
import operator as op

def FMM(t, s, f, r, binom):

    ''''
    t = n-vector of targets
    s = n-vector of sources
    f = n-vector of points
    r = number of terms in Taylor expansion
    binom = Pascal's triangle for fast binomial coefficient lookup
    '''

    n = range(len(t))
    L = int(math.log(len(t), 2))
    p = [0 for _ in n]
    sDict = {}
    tDict = {}
    cellDict = {}
    r = range(r+1)
    
    'main loop'
    for l in range(2, L-1):
        power = int(pow(2,l))
        interval = 1.0/power
        cells = range(power)

        sDict = {key: 0 for key in s}
        tDict = {key: 0 for key in t}
        cellDict = {key: [[],[]] for key in cells}

        'sort s and t into cells and construct reverse indices'
        for i in t:
            cellNum = int(i/interval)
            tDict[i] = cellNum
            cellDict[cellNum][1].append(i)
        for j in s:
            cellNum = int(j/interval)
            sDict[j] = cellNum
            cellDict[cellNum][0].append(j)

        'mathemagic starts here'
        'construct moments'
        S = [[0 for i in r] for l in cells]
        for j in n:
            s_j = s[j]
            scell = sDict[s_j]
            sigma = interval*scell + (interval/2)
            for m in r:
                S[scell][m] += pow(s_j - sigma, m) * f[j]

        'calculate interactions'
        T = [[0 for i in r] for j in cells]
        for tcell in cells:
            possibles = [c for c in range(tcell-3, tcell+4) if c >= 0 and c < power]
            for scell in possibles:
                if abs(tcell-scell) > 1 and abs(tcell/2 - scell/2) <= 1:
                    sigma = interval*scell + interval/2
                    tau = interval*tcell + interval/2
                    for m in r:
                        for k in r:
                            T[tcell][m] += binom[m+k][k] * pow(tau-sigma, -m-k-1) *S[scell][k]
                            
        'increment return values'
        for i in n:
            t_i = t[i]
            tcell = tDict[t_i]
            tau = interval*tcell + interval/2
            for m in r:
                p[i] += T[tcell][m] * pow(tau - t_i, m)

    'local interactions'
    for i in n:
        t_i = t[i]
        tcell = tDict[t_i]
        for scell in range(tcell-1, tcell+2):
            if scell == -1 or scell == power:
                continue
            sources = cellDict[scell][0]
            if not sources:
                continue
            for s_j in sources:
                j = s.index(s_j)
                p[i] += f[j]/(t_i - s_j)
                
    return p

def cauchy(t, s):
    toReturn = [[0 for _ in range(len(t))] for _ in range(len(s))]
    for i in range(len(t)):
        for j in range(len(s)):
            toReturn[i][j] = 1.0/(t[i] - s[j])
    return toReturn

def matmult(A, x):
    toReturn = [0 for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(x)):
            toReturn[i] += A[i][j] * x[j]
    return toReturn

def twonorm(v):
    return math.sqrt(math.fsum(map(lambda x: x**2, v)))
    
def pascal(n):
    toReturn = []
    for _ in range(n):
        row = [1]
        if toReturn:
            last = toReturn[-1]
            row.extend([sum(pair) for pair in zip(last, last[1:])])
            row.append(1)
        toReturn.append(row)
    return toReturn

if __name__ == "__main__":
    nList = []
    t10List = []
    t20List = []
    t30List = []
    timeList = []
    error10List = []
    error20List = []
    error30List = []
    binom = pascal(61)
    
    for l in range(7, 14):
        n = int(pow(2,l))
        print n
        nList.append(n)
        t = [random.random() for _ in range(n)]
        s = [random.random() for _ in range(n)]
        f = [random.random() for _ in range(n)]
        
        tic = time.clock()
        p10 = FMM(t, s, f, 10, binom)
        t10 = time.clock() - tic
        t10List.append(t10)
        print t10

        tic = time.clock()
        p20 = FMM(t, s, f, 20, binom)
        t20 = time.clock() - tic
        t20List.append(t20)
        print t20

        tic = time.clock()
        p30 = FMM(t, s, f, 30, binom)
        t30 = time.clock() - tic
        t30List.append(t30)
        print t30

        cauchymat = cauchy(t,s)
        tic = time.clock()
        exact = matmult(cauchymat, f)
        etime = time.clock() - tic
        timeList.append(etime)
        print etime
        
        error10 = twonorm([abs(p10[i] - exact[i]) for i in range(n)])
        error20 = twonorm([abs(p20[i] - exact[i]) for i in range(n)])
        error30 = twonorm([abs(p30[i] - exact[i]) for i in range(n)])
        norm = twonorm(exact)
        error10List.append(error10/norm)
        error20List.append(error20/norm)
        error30List.append(error30/norm)

    print nList
    print t10List
    print t20List
    print t30List
    print timeList
    print error10List
    print error20List
    print error30List
