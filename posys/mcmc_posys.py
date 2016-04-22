# -*- coding: utf-8 -*-
"""
Urban Observatory
"""

import numpy as np
import copy as cp
import pypower.api as pypo 
from get_ppc14 import *
import matplotlib.pyplot as plt
import emcee
import corner

# Refer to https://github.com/rwl/PYPOWER/blob/master/pypower/newtonpf.py

def likelihood_ps(theta,y):
    
    # -- modify load buses
    ppc["bus"][ind,2] = theta
    
    # -- estimate the transformer measurements
    ppopt = pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0) 
    sol   = pypo.runpf(ppc, ppopt)
    estim = sol[0]['gen'][:,2] 
    
    # -- calculate the likelihood
    sig     = 1.0
    if (theta>=0.0).all():
        return -((estim - y)**2).sum()/(2*sig**2)
    else:
        return -np.inf


# -- utilities
ndim     = 9
nwalkers = 30
nsteps   = 200
cut      = 50


# -- intialize the 14-bus system
ppc0  = get_ppc14(1,0,1) # 0 implies no change
ppc   = cp.deepcopy(ppc0)
y     = ppc0['gen'][:,2].copy() # default measured values of transformers
ind   = ppc0["bus"][:,1]==1 # building indices
binit = ppc0["bus"][ind,2].copy()

# -- Initialize sampler
print("initializing sampler...")
np.random.seed(314)
sampler = emcee.EnsembleSampler(nwalkers, ndim, likelihood_ps, args=[y])
pos     = [binit + 20.0*np.random.randn(ndim) for i in range(nwalkers)]
pos    *= np.transpose(np.transpose(pos) >0.00)


# -- run walkers
print("running walkers...")
sampler.run_mcmc(pos, nsteps)
print("walkers finished...")


# -- save chain
#oname = "../output/random_test3.npy"
#print("saving chain to {0}".format(oname))
#np.save(oname,sampler.chain)
