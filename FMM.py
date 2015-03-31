#FMM.py
import math
import time
import random
import operator as op

def FMM(t, s, f, r, binom):

    n = len(t)
    L = int(math.log(n, 2))
    p = [0 for _ in range(n)]
    sDict = {}
    tDict = {}
    cellDict = {}
    powers = range(r+1)
    
    'main loop'
    for l in range(2, L+1):
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

        'construct moments'
        S = [[0 for i in powers] for l in cells]
        for cell in cells:
            sigma = interval*cell + (interval/2)
            sources = cellDict[cell][0]
            for m in powers:
                for s_j in sources:
                    j = s.index(s_j)
                    S[cell][m] += pow(s_j - sigma, m) * f[j]
                    
        #tic = time.clock()
        'calculate interactions'
        T = [[0 for i in powers] for j in cells]
        for tcell in cells:
            for scell in cells:
                if abs(tcell-scell) > 1 and abs(tcell/2 - scell/2) <= 1:
                    sigma = interval*scell + interval/2
                    tau = interval*tcell + interval/2
                    for m in powers:
                        for k in powers:
                            combo = binom[m+k][k]
                            center = pow(tau-sigma, -m-k-1)
                            T[tcell][m] += combo * center *S[scell][k]
        #print time.clock() - tic
       
        'increment return values'
        for t_i in t:
            i = t.index(t_i)
            tcell = tDict[t_i]
            tau = interval*tcell + interval/2
            for m in powers:
                p[i] += T[tcell][m] * pow(t_i - tau, m)

    'local interactions'
    for t_i in t:
        tcell = tDict[t_i]
        i = t.index(t_i)
        for scell in range(tcell-1, tcell+2):
            if scell == -1 or scell == power:
                continue
            for s_j in cellDict[scell][0]:
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
    errorList = []
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
        maxError = max(error10/norm, error20/norm, error30/norm)
        errorList.append(maxError)

    print nList
    print t10List
    print t20List
    print t30List
    print timeList
    print errorList
