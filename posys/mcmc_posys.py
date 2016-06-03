# -*- coding: utf-8 -*-

# Urban Observatory

import numpy as np
import copy as cp 
import pypower.api as pypo 
from get_ppc14 import *
import matplotlib.pyplot as plt
import emcee
import corner

# Refer to https://github.com/rwl/PYPOWER/blob/master/pypower/newtonpf.py

def lnlike(theta,y):
    
    # -- modify load buses
    ppc["bus"][ind,2] = theta
    
    # -- estimate the transformer measurements 
    sol   = pypo.runpf(ppc, pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0))
    estim = sol[0]['gen'][:,2] 
    
    # -- calculate the likelihood
    sig     = 1.0
    if (theta>=0.0).all():
        return -((estim - y)**2).sum()/(2*sig**2)
    else:
        return -np.inf

# -- utilities
ndim     = 9
nwalkers = 300
nsteps   = 100
cut      = 50

# -- intialize the 14-bus system
ppc0  = get_ppc14(1,0,1) # 0 implies no change
ppc   = cp.deepcopy(ppc0)
#y     = ppc0['gen'][:,2].copy() # default measured values of transformers
ind   = ppc0["bus"][:,1]==1 # building indices
#binit = ppc0["bus"][ind,2].copy()

for val in [1.00, 2.00, 5.00, 10.00, 20.00]:
	# -- Initialize sample
	print("initializing sampler...")
	np.random.seed(314)
	ppc["bus"][ind,2] = val*np.ones(ndim)
	sol     = pypo.runpf(ppc, pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0))
	y       = sol[0]['gen'][:,2]
	binit   = ppc['bus'][ind,2].copy()	
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike, args=[y])
	pos     = np.array([binit*(1.0+0.2*np.random.randn(ndim)) for i in range(nwalkers)]).clip(min=0.0)
	# -- run walkers
	print("running walkers...")
	sampler.run_mcmc(pos, nsteps)
	print("walkers finished...")

	# -- save chain
	oname = "./output/random_test20percent_diffValues"+str(val)+".npy"
	print("saving chain to {0}".format(oname))
	np.save(oname,sampler.chain)
