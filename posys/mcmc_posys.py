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

#ndim, nwalkers = 9, 300
ndim, nwalkers = 9, 300
nsteps = 500
cut = 50

# Refer to https://github.com/rwl/PYPOWER/blob/master/pypower/newtonpf.py

def likelihood_ps(theta,y):
    
    # -- modify load buses
    ppc["bus"][ind,2] = theta
    
    # -- estimate the transformer measurements
    ppopt   = pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0) 
    r       = pypo.runpf(ppc, ppopt)
    estim   = r[0]['gen'][:,2] 
    
    # -- calculate the likelihood
    sig     = 10.0
    if (theta>=0.0).all():
        return np.exp(-((estim - y)**2).sum()/(2*sig**2))
    else:
        return -np.inf

# Choose the "true" parameters.
ppc0       = get_ppc14(op_change=1,dlt=0,busN=1) #trivial case: original solutin
ppc = cp.deepcopy(ppc0)
y = ppc0['gen'][:,2].copy() # Measured values in the transformers
ind   = ppc0["bus"][:,1]==1 
theta = []


for nit in range(0,5):
    for val in [20.00, 50.00]:
        # -- Inputting loads
        np.random.seed(314)
        del(theta)
        theta  = [val*np.ones(ndim)]#ppc0["bus"][ind,2].copy()#+ 1e-0*np.random.randn(ndim)
        theta *= theta>0.0
    
        # -- Initializing mcmc walkers
        #pos = [theta + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        #pos = [np.abs(val*np.ones(ndim)+ 2e+1*np.random.randn(ndim))  for i in range(nwalkers)]    
        pos = np.random.uniform(0,100,size=[nwalkers,ndim])       
        
        print("initializing sampler...")
        sampler = emcee.EnsembleSampler(nwalkers, ndim, likelihood_ps, args=[y])
        
        print("running walkers...")
        #sampler.run_mcmc(pos, 500)
        sampler.run_mcmc(pos, nsteps)
        print("walkers finished...")
        
        #samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
        #samples = sampler.chain[:, :, :].reshape((-1, ndim))
        nm = '%sWalk%sbuild%ssteps%svalRandomInit%s'%(nwalkers,ndim,nsteps,val,nit)
        np.save('random/'+nm, sampler.chain)
        
        np.save('random/median_'+nm,np.median(np.median(sampler.chain[:,cut:,:],axis=0),axis=0) - theta)
    

#fig = corner.corner(samples, labels=["$b1$", "$b2$", "$b3$","$b4$","$b5$",
#                                     "$b6$","$b7$","$b8$","$b9$"],
#                    truths=theta)
#                    


#for fl in Walk9build500steps1,Walk9build500steps2, Walk9build500steps50,Walk9build500steps5:
    
    


"""
# Choose the "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534
# Generate some synthetic data from the model.
N = 50
x = np.sort(10*np.random.rand(N))
yerr = 0.1+0.5*np.random.rand(N)
y = m_true*x+b_true
y += np.abs(f_true*y) * np.random.randn(N)
y += yerr * np.random.randn(N)
A = np.vstack((np.ones_like(x), x)).T
C = np.diag(yerr * yerr)
cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))
def lnlike(theta, x, y, yerr):
    m, b, lnf = theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))
import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
m_ml, b_ml, lnf_ml = result["x"]
def lnprior(theta):
    m, b, lnf = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
        return 0.0
    return -np.inf
    
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)
    
ndim, nwalkers = 3, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    
import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
import corner
fig = corner.corner(samples, labels=["$m$", "$b$", "$\ln\,f$"],
                      truths=[m_true, b_true, np.log(f_true)])
fig.savefig("triangle.png")
samples[:, 2] = np.exp(samples[:, 2])
m_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))
"""                                                