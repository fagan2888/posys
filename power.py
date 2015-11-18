# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 09:58:16 2015

@author: Sara
"""
import numpy as np

def power_uo(GG,BB,N,NG,NL,vm,an,e,f,pgspec,plspec,qlspec):
    """
    Calculate the real and reactive power (P and Q) as well as Delta P and 
    Delta Q
    """

    # -- calculate P and Q
    a_sum = (GG*e - BB*f).sum(1)
    b_sum = (GG*f + BB*e).sum(1)
    P     = e*a_sum + f*b_sum
    Q     = f*a_sum - e*b_sum
    
    # -- calculate the Deltas
    pb = np.concatenate([pgspec[:NG-1] - P[1:NG],  
                         plspec[:NL] - P[NG:NG+NL],
                         qlspec[:NL] - Q[NG:NG+NL]])

    return P, Q, pb
   
