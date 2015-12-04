# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 08:40:43 2015

"""
from __future__ import division
import numpy as np

def jaco(GG,BB,P,Q,N,NG,NL,vm,e,f,ind_gen,ind_load):
    """
    Calculate the Jacobian
    """

    # -- utilities
    iden = np.eye(len(f))
    inds = iden>0
    a    = GG*e - BB*f
    b    = GG*f + BB*e
    ep   = iden*e
    fp   = iden*f
    ep_a = np.dot(ep,a)
    ep_b = np.dot(ep,b)
    fp_a = np.dot(fp,a)
    fp_b = np.dot(fp,b)

    # Forumate H, N, J and L matrices with standard sizes
    HH       =  fp_a - ep_b
    NN       =  ep_a + fp_b
    JJ       =  - ep_a - fp_b
    LL       =  fp_a - ep_b
    HH[inds] = -Q - BB[inds]*vm**2 ## CAREFUL
    NN[inds] =  P + GG[inds]*vm**2
    JJ[inds] =  P - GG[inds]*vm**2
    LL[inds] =  Q - BB[inds]*vm**2
    
    return np.vstack((np.hstack((HH[1:NG+NL,1:NG+NL],NN[1:NG+NL,ind_load])),
                      np.hstack((JJ[ind_load,1:NG+NL], LL[ind_load][:,ind_load] ))))
