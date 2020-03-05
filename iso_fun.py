import numpy as np
from itertools import combinations

def deut_up(t,kint,lnP):
    '''Calculates deuterium uptake at time t given kint and lnP'''
    P = np.exp(lnP)
    tmp = 0
    namide = 0
    for i in range(0,len(kint)):
        if kint[i] >= 0:
            namide += 1
            tmp += np.exp(-kint[i]/P[i]*t*60)
    return (namide - tmp)/namide

def single_deut_up(t,kint,lnP):
    '''Calculates deuterium uptake at time t given kint and lnP'''
    P = np.exp(lnP)
    if kint > 0:
        return 1-np.exp(-kint/P*t*60)
    else:
        return 0

def set_A_group(a,k):
    return list(combinations(a,k))

def prob(k,t,kint,lnP):
    ''' Calculates the probability that k residues have exchanged at time t
    given intrinsic rates kint and protection factors lnP '''
    A = set_A_group([i for i in range(len(kint))],k)
    probability = 0
    for element in A:
        d1 = []; d2 = []
        for j in range(len(kint)):
            if j in element:
                d1.append(single_deut_up(t,kint[j],lnP[j]))
            else:
                d2.append(1-single_deut_up(t,kint[j],lnP[j]))
        p1 = np.prod(d1)
        p2 = np.prod(d2)
        probability += p1*p2
    return probability

def isotopic_envelope(t,kint,lnP):
    k_values = [i for i in range(len(kint))]
    freq = []
    for k in k_values:
        freq.append(prob(k,t,kint,lnP))
    return freq

def centered_isotopic_envelope(t,kint,lnP,fr0):   
    fr = isotopic_envelope(t, kint, lnP)
    f = np.zeros(len(kint))
    for i in range(len(f)):
        for j in range(i+1):
            f[i] += fr0[i-j]*fr[j]
    return f

# %%

