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
from power import power_uo
from utilities import get_params
from jaco import jaco
from network import get_network
from ybus import ybus

# Be clear on initial conditions

bus,line = get_network()
BB, GG = ybus(bus,line)

params = get_params()
GG, BB, Q_d, P_d, N, NG, NL, vm, an, e, f, pgspec, plspec, qlspec = params

# -- set tolerance params
min_tol  = 0.005
itr      = 1
iter_max = 5
tol      = 1

while (itr <= iter_max and tol > min_tol): 
    #Compute powers and power mismatches
    P, Q, dpq = power_uo(GG,BB,N,NG,NL,e,f,pgspec,plspec,qlspec)
    jac       = jaco(GG,BB,P,Q,N,NG,NL,vm,e,f)

    # Solve system of equations
    dtv = np.linalg.solve(jac, dpq)

    # Update values of voltages and angles
    vm[NG:NG+NL] *= 1 + dtv[NL+1:N+NL]
    an[1:N]      += dtv[0:N-1]*180/np.pi
    
    # Update vectors e and f
    e = vm*np.cos(an*np.pi/180.)
    f = vm*np.sin(an*np.pi/180.)

    # Refresh
    itr += 1

    print("iter, theta_2, voltage_3, theta_3 = {0}, {1}, {2}, {3}" \
              .format(itr,an[1],vm[2],an[2]))
