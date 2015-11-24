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
            N      = 3
            NG     = 2
            NL     = 1
            vm     = np.array([1.0,1.1249,0.93834])
            an     = np.array([0.0, 6.3, -3.44])
            e      = vm*np.cos(an*np.pi/180.)
            f      = vm*np.sin(an*np.pi/180.)
            pgspec = np.array([1.7])
            plspec = np.array([-2])
            qlspec = np.array([-1])
            
    elif setting == '14 bus':
            bus,line = get_network()
            BB, GG = ybus(bus,line)            
            NG = sum(bus[:,9]!=3)
            NL = sum(bus[:,9]==3)
            pgspec = bus[bus[:,9]==2,3]
            plspec = bus[bus[:,9]==3,5]
            qlspec = bus[bus[:,9]==3,6]
            # vgspec
            N = NG - 1 + 2*NL
            vm = bus[:,1]
            an = bus[:,2]
            e  = vm*np.cos(an*np.pi/180.)
            f  = vm*np.sin(an*np.pi/180.)

    return [GG, BB, N, NG, NL, vm, an, e, f, pgspec, plspec, qlspec]
    # GGD: what do these variables stand for?
    # GGD: also change the names to something more readable