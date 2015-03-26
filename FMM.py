import math
import random
import numpy as np

def FMM(t, s, f, r):

    n = len(t)

    if not n is len(s):
        return None
    if not len(s) is len(f):
        return None
    
    L = int(math.log(n, 2))
    p = np.zeros(n)
    sDict = {}
    tDict = {}
    cellDict = {}
    
    'main loop'
    for l in range(2, L+1):
        sDict = {key: 0 for key in s}
        tDict = {key: 0 for key in t}
        cellDict = {key: [[],[]] for key in range(int(pow(2,l)))}

        power = int(pow(2,l))
        interval = 1.0/power

        'sort s and t into cells and construct reverse indices'
        for i in t:
            temp = interval
            cellNum = 0
            while i > temp:
                temp += temp
                cellNum += 1
            tDict[i] = cellNum
            cellDict[cellNum][1].append(i)
        for j in s:
            temp = interval
            cellNum = 0
            while j > temp:
                temp += temp
                cellNum += 1
            sDict[j] = cellNum
            cellDict[cellNum][0].append(j)

        'construct moments'
        S = [[0 for i in range(0, r)] for l in range(power)]
        for cell in range(power):
            sigma = interval*cell + (interval/2)
            for m in range(r+1):
                sInCell = cellDict[cell][0]
                for source in sInCell:
                    index = s.index(source)
                    S[cell][m] += pow(source - sigma, m) * f[index]

        'calculate interactions'
        T = [0 for k in range(p+1)]
        for k in range(p+1):
            for scell in range(power):
                sigma = interval*scell + (interval/2)
                for tcell in range(power):
                    tau = interval*tcell + (interval/2)
                    if abs(tcell - scell) > 1:
                        series = 0
                        for m in range(p+1):
                            combo = math.factorial(m+k)/math.factorial(k)/math.factorial(m)
                            center = pow(tau-sigma, -k-m-1)
                            series += center*combo* S[scell][m]
                        T[k] += series

        'increment return values'
        for tcell in range(power):
            tau = interval*tcell + (interval/2)
            for t_i in cellDict[tcell][1]:
                index = t.index(t_i)
                series = 0
                for k in range(r+1):
                    series += T[k] * pow(t_i - tau, k)
                p[index] += series

    'local interactions'
    for scell in range(1,int(pow(2,L))-1):
        for tcell in range(scell-1,scell+2):
            for t_i in cellDict[tcell][1]:
                index = t.index(t_i)
                series = 0
                for s_j in cellDict[scell][0]:
                    j = s.index(s_j)
                    series += f[j]/(t_i - s_j)
                p[index] += series

    return p

if __name__ == "__main__":
    for l in range(7, 14):
        n = int(pow(2,l))
        
