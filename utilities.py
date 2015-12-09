# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 17:47:04 2015

@author: Sara Arango and Gregory Dobler
"""
import numpy as np
from ybus import ybus
from network import get_network

# Sample input for a very simple network

def get_params(setting):
    """
    Return the parameters for the system.
    """
    if setting == '3 bus': 
            GG     = np.array([[4,0,-4],
                       [0,4,-4],
                       [-4,-4,8]])
            BB     = np.array([[-5,0,5], 
                       [0,-10,10],
                       [5,10,-15]])
            #Q_d    = np.array([0, 0, -0.9397])
            #P_d    = np.array([0, 2.6860,-2.9313])
            NG     = 2
            NL     = 1
            N      = NG - 1 + 2*NL
            #vm     = np.array([1.0,1.1249,0.93834])
            vm     = np.array([1.0, 1.0, 1.0])
            #vm     = np.array([1.0,1.1249,1.0])
            #an     = np.array([0.0, 6.3, -3.44])
            an     = np.array([0.0,0.0,0.0])  
            #an     = np.array([0.0, 0.0, 0.0])
            e      = vm*np.cos(an*np.pi/180.)
            f      = vm*np.sin(an*np.pi/180.)
#            pgspec = np.array([1.7])
#            plspec = np.array([-2])
#            qlspec = np.array([-1])
            pgspec = np.array([0.0])
            plspec = np.array([0.0])
            qlspec = np.array([0.0])
            ind_gen = np.array([ False, True, False])
            ind_load = np.array([False, False, True])
            
    elif setting == '14 bus':
            bus, line = get_network()
            GG, BB    = ybus(bus,line)   
            ind_gen   = np.array(bus[:,9]==2)
            ind_load  = np.array(bus[:,9]==3)
            NG        = sum(ind_gen)+1
            NL        = sum(ind_load)
            pgspec    = bus[ind_gen,3]
            plspec    = bus[ind_load,5]
            qlspec    = bus[ind_load,6]
            # vgspec
            N         = NG - 1 + 2*NL
            vm        = bus[:,1]
            an        = bus[:,2]
            e         = vm*np.cos(an*np.pi/180.)
            f         = vm*np.sin(an*np.pi/180.)
 

    return [GG, BB, N, NG, NL, vm, an, e, f, pgspec, plspec, qlspec, ind_gen, 
            ind_load]
    # GGD: what do these variables stand for?
    # GGD: also change the names to something more readable
