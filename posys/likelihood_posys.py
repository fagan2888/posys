# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 09:26:19 2016

@author: saf537
"""
import numpy as np
import case14_mod
import pypower.api as pypo 

def likelihood_ps(bus_vec,measur_vec):      
    ppc     = case14_mod.case14_mod(busN = 0, dlt = 0, op_change=2, invec = measur_vec)
    ppopt   = pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0) # Careful: have to use Newton Method!!!
    r       = pypo.runpf(ppc, ppopt)
    estim   = r[0]['gen'][:,2]  
    sig     = np.std(estim - measur_vec)
    return np.exp(((estim - measur_vec)**2).sum()/sig**2)        
        
ppc0       = case14_mod.case14_mod(busN = 1,dlt = 0, op_change=1) # trivial case: original solutin
measur_vec = ppc0['gen'][:,2]
bus_vec    = measur_vec + np.random.normal(0,1,len(measur_vec))
L          = likelihood_ps(bus_vec,measur_vec)