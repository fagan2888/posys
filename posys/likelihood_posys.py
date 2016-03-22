# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 09:26:19 2016

@author: saf537
"""
import numpy as np
import pypower.api as pypo 
from get_ppc14 import *
import matplotlib.pyplot as plt

#def likelihood_ps(bus_vec,measur_vec):      
def likelihood_ps(measur_vec,bvec):

    # -- modify load buses
    ppc["bus"][ind,2] = bvec

    # -- estimate the transformer measurements
    #ppc     = get_ppc14(op_change=2,dlt=0,busN=0,invec=bus_vec)
    ppopt   = pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0) # Careful: have to use Newton Method!!!
    r       = pypo.runpf(ppc, ppopt)
    estim   = r[0]['gen'][:,2] 
    
    # -- calculate the likelihood
    sig     = np.std(estim - measur_vec)
    return np.exp(((estim - measur_vec)**2).sum()/sig**2)        
        
ppc0       = get_ppc14(op_change=1,dlt=0,busN=1) #trivial case: original solutin
measur_vec = ppc0['gen'][:,2]
ind   = ppc0["bus"][:,1]==1 # 1 is load bus
np.random.seed(314)

# -- modify the load buses

#loads = ppc0['bus'][:,0][ppc0['bus'][:,1]==1]
tmp = 0
LK_vec = []
for bs in loads:
#    L          = likelihood_ps(measur_vec,bvec)
    #bvec  = ppc0["bus"][ind,2] + 0.5*np.random.randn(ind.sum())
    L_vec = []
    for mod in np.arange(0,1,0.05):
        bvec  = ppc0["bus"][ind,2]
        bvec[tmp] =  mod*np.random.randn()
        bvec *= bvec>0.0
        L     = likelihood_ps(measur_vec,bvec)
        L_vec.append(L)
    plt.plot(L_vec)
    LK_vec.append(L_vec)
    tmp += 1
    
#    print L
    
    