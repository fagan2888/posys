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

ndim, nwalkers = 9, 300

def likelihood_ps(theta,y):
    # -- modify load buses
    ppc["bus"][ind,2] = theta
    # -- estimate the transformer measurements
    ppopt   = pypo.ppoption(PF_ALG=2, VERBOSE=0, OUT_ALL=0) 
    r       = pypo.runpf(ppc, ppopt)
    estim   = r[0]['gen'][:,2] 
    # -- calculate the likelihood
    sig     = 10.0
    return np.exp(-((estim - y)**2).sum()/(2*sig**2))

# Choose the "true" parameters.
ppc0       = get_ppc14(op_change=1,dlt=0,busN=1) #trivial case: original solutin
ppc = cp.deepcopy(ppc0)
y = ppc0['gen'][:,2].copy()
ind   = ppc0["bus"][:,1]==1 
np.random.seed(314)
theta  = ppc0["bus"][ind,2].copy()+ 1e-4*np.random.randn(ndim)
theta *= theta>0.0

pos = [theta + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, likelihood_ps, args=y)

sampler.run_mcmc(pos, 500)

samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

"""

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