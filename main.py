# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 08:36:20 2015

@author: Francisco de León, transcribed by Sara Arango and Gregory Dobler

Inputs: from input file:

- GG and BB 
- P and Q
- N, NG and NL
- vm and an 
- e and f 
- pgspec, plspec and qlspec 

"""
#import numpy as np
#from power import power_uo
#from utilities import get_params
#from jaco import jaco
from matplotlib import pyplot
execfile('power.py')
execfile('utilities.py')
execfile('jaco.py')

# Be clear on initial conditions
setting = '3 bus'
params = get_params(setting)
GG, BB, N, NG, NL, vm, an, e, f, pgspec, plspec, qlspec, ind_gen, ind_load = params

# -- set tolerance params
min_tol  = 0.005
itr      = 0
iter_max = 100
tol      = 1
delta_magnitude = [0]*(iter_max+1)
dpq_magnitude = [0]*(iter_max+1)

P, Q, dpq = power_uo(GG,BB,N,NG,NL,e,f,pgspec,plspec,qlspec,ind_gen,ind_load)
jac       = jaco(GG,BB,P,Q,N,NG,NL,vm,e,f,ind_gen,ind_load)

while (itr <= iter_max and tol > min_tol): 
    # Compute powers and power mismatches
    P, Q, dpq = power_uo(GG,BB,N,NG,NL,e,f,pgspec,plspec,qlspec,ind_gen,
                         ind_load)
    #jac       = jaco(GG,BB,P,Q,N,NG,NL,vm,e,f,ind_gen,ind_load)
                        
    # Solve system of equations
    dtv = np.linalg.solve(jac, dpq)
#    dtv = np.linalg.lstsq(jac, dpq)[0]
#    dtv = np.dot(np.linalg.pinv(np.dot(jac,jac.T)),np.dot(dpq,jac.T))

    delta_magnitude[itr] = np.linalg.norm(dtv,1)
    dpq_magnitude[itr]   = np.linalg.norm(dpq,1)
    # Update values of voltages and angles
    vm[ind_load] *= 1 + dtv[NG:NG+NL] # load bus
    an[1:N]      += dtv[0:NL+NG-1]*180/np.pi # loads and generators
    
    vm[ind_load] = np.abs(vm[ind_load])
    an[1:N]      = an[1:N] % 360.

    # Update vectors e and f
    e = vm*np.cos(an*np.pi/180.)
    f = vm*np.sin(an*np.pi/180.)
    
    # Refresh
    itr += 1
    
    print np.max(np.abs(dpq))
    #print("iter, theta_2, voltage_3, theta_3 = {0}, {1}, {2}, {3}" \
    #          .format(itr,an[1],vm[2],an[2]))

pyplot.plot(dpq_magnitude)
pyplot.plot(delta_magnitude)
