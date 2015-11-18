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
import numpy as np
import jaco as jaco
from math import pi
from power import power_uo

# Be clear on initial conditions
from input import GG, BB, N, NG, NL, vm, an, e, f, plspec, qlspec, pgspec, Q_d,P_d

min_tol = 0.005
n_ite  = 1
ite_max = 1
tol = 1

#Compute powers and power mismatches
P,Q,Delta_known = power_uo(GG,BB,N,NG,NL,vm,an,e,f,pgspec,plspec,qlspec)

Jac = jaco.JACO(GG,BB,P,Q,N,NG,NL,vm,an,e,f)

while (n_ite <= ite_max and tol > min_tol): 
    # Define jacobian
    Q = np.transpose(np.array([0, 1.5365, -0.9397])) # Q that makes the jacobian work
    P = np.transpose(np.array([0, 0,-2.9313]))
    Jac = jaco.JACO(GG,BB,P,Q,N,NG,NL,vm,an,e,f)
    
    # Solve system of equations
    Delta_unknown = np.linalg.solve(Jac, Delta_known)
    
    # Update values of voltages and angles
    vm_p = Delta_unknown[NL+1:N+NL] #setup the subsets correctly
    an_p = Delta_unknown[0:N-1]*180/pi #setup the subsets correctly
    
    # Calculate new estimates for the unknowns
    vm[NG:NG+NL] = vm[N-1]*(1 + vm_p)
    an[1:N] = an[1:N] + an_p
    
    # Update vectors e and f
    e = np.transpose(vm*np.cos(np.radians(an)))
    f = np.transpose(vm*np.sin(np.radians(an)))
    
    # Compute powers and power mismatches
    P,Q,Delta_known = power_uo(GG,BB,N,NG,NL,vm,an,e,f,pgspec,plspec,qlspec)   
    
    # Refresh
    n_ite += 1
    print vm
    print an
    #@er_vol = np.sum(np.abs(np.abs(vm) - np.abs(vm_p))/np.abs(vm_p))
    #er_an = np.sum(np.abs(np.abs(an) - np.abs(an_p))/np.abs(an_p))
    tol = 1#np.min([er_vol, er_an])
