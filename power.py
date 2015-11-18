# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 09:58:16 2015

@author: Sara
"""

import numpy as np

def power_uo(GG,BB,N,NG,NL,vm,an,e,f,pgspec,plspec,qlspec):
    P = [0]*N
    Q = [0]*N
    for i in range(0, N):
       s1 = 0.0
       s2 = 0.0 
       for j in range(0, N):
           s1 = s1 + GG[i,j] * e[j] - BB[i,j] * f[j] 
           s2 = s2 + GG[i,j] * f[j] + BB[i,j] * e[j] 
       #print s1           
       P[i] = e[i] * s1 + f[i] * s2 
       Q[i] = f[i] * s1 - e[i] * s2
    P = np.transpose(P)
    Q = np.transpose(Q)
   
   # Calculate the deltas
    pb = [0]*(NG+NL)
    for i in range(0,NG-1):
        pb[i] = pgspec[i] - P[i+1] 
    for i in range(0,NL):
        j = NG + i - 1
        k = NL + j
        print j
        print k
        pb[j] = plspec[i] - P[NG+i]
        pb[k] = qlspec[i] - Q[NG+i]
        pb = np.transpose(pb)
    return P,Q,pb
   