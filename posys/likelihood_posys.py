# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 09:26:19 2016

@author: saf537
"""
import numpy as np
import copy as cp
import pypower.api as pypo 
from get_ppc14 import *
import matplotlib.pyplot as plt
#import pdb

#def likelihood_ps(bus_vec,measur_vec):      
def likelihood_ps(measur_vec,bvec):

    # -- modify load buses
    ppc["bus"][ind,2] = bvec

    # -- estimate the transformer measurements
    ppopt   = pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0) 
    r       = pypo.runpf(ppc, ppopt)
    estim   = r[0]['gen'][:,2] 
    
    # -- calculate the likelihood
    sig     = 10.0

    return np.exp(-((estim - measur_vec)**2).sum()/(2*sig**2))
        
ppc0       = get_ppc14(op_change=1,dlt=0,busN=1) #trivial case: original solutin
ppc = cp.deepcopy(ppc0)
measur_vec = ppc0['gen'][:,2]
ind   = ppc0["bus"][:,1]==1 
np.random.seed(314)

# -- modify the load buses

loads = ppc0['bus'][:,0][ppc0['bus'][:,1]==1]
tmp = 0
LK_vec = []
for bs in loads:
    L_vec = []
    for mod in np.arange(0,2.1,0.1):
        bvec  = ppc0["bus"][ind,2].copy()
        bvec[tmp] *=  mod
        bvec *= bvec>0.0
        L     = likelihood_ps(measur_vec,bvec)
        L_vec.append(L)
    plt.plot(L_vec)
    LK_vec.append(L_vec)
    tmp += 1
    
#    print L  