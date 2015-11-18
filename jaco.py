# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 08:40:43 2015

@author: dr. Francisco de León, transcribed to Python by Sara Arango and dr. Gregory Dobler.
"""
from __future__ import division
import numpy as np
from input import GG, BB, N, NG, NL, vm, an, e, f

def JACO(GG,BB,P,Q,N,NG,NL,vm,an,e,f):

    a = GG * e - BB * f
    b = GG * f  + BB * e    
#    a = np.array([[0,0,0],[0,0,-3.1836],[0,-5.7068, 0]])
#    b = np.array([[0,0,0],[0,0,9.59],   [0,10.68,0]])    
    # Forumate H, N, J and L matrices with standard sizes
    di1 = np.diag_indices(N)
    fp = np.eye(N)
    ep = np.eye(N)
    fp[di1]=f
    ep[di1]=e
    
    HH = np.dot(fp,a) - np.dot(ep,b)
    NN = np.dot(ep,a) + np.dot(fp,b)
    di = np.diag_indices(N)
    HH[di] = -Q - BB[di]*vm**2
    NN[di] = P + GG[di]*vm**2 # diagonal
    
    JJ = -np.dot(ep,a) - np.dot(fp,b)
    LL = np.dot(fp,a) - np.dot(ep,b)
    JJ[di] = P - GG[di]*vm**2 # diagonal
    LL[di] = Q - BB[di]*vm**2 # diagonal
  
    H  = HH[1:N,1:N]
    L  = LL[NG:N,NG:N] 
    J  = JJ[NG:N,1:N]
    NN = NN[NG-1:N,NG:N]
    Jac = np.vstack((np.hstack((H,NN)),np.hstack((J,L))))
    return Jac
#Jac.flatten()#.reshape((N,N))