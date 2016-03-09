# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 09:26:19 2016

@author: saf537
"""
import numpy as np

def likelihood_ps(bus_vec,measur_vec):      
    ppc     = case14_mod.case14_mod(busN = i,dlt = dlt_vec[j])
    ppopt   = pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0) # Careful: have to use Newton Method!!!
    r       = pypo.runpf(ppc, ppopt)
    estim   =
    sig     =  
    return np.exp(()**2/sig**2)        
    