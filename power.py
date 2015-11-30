# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 09:58:16 2015

@author: Sara
"""
import numpy as np

def power_uo(GG,BB,N,NG,NL,e,f,pgspec,plspec,qlspec,ind_gen,ind_load):
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
    dtv   = np.concatenate([pgspec[:NG-1] - P[ind_gen],  
                            plspec[:NL] - P[ind_load],
                            qlspec[:NL] - Q[ind_load]])
                         
    return P, Q, dtv
   
