# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 17:47:04 2015

@author: Sara Arango and Gregory Dobler
"""
import numpy as np
# Sample input for a very simple network

def get_params():
    """
    Return the parameters for the system.
    """
    # GGD: what do these variables stand for?
    # GGD: also change the names to something more readable
    GG     = np.array([[4,0,-4],
                       [0,4,-4],
                       [-4,-4,8]])
    BB     = np.array([[-5,0,5], 
                       [0,-10,10],
                       [5,10,-15]])

    Q_d    = np.array([0, 0, -0.9397])
    P_d    = np.array([0, 2.6860,-2.9313])
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

    return [GG, BB, Q_d, P_d, N, NG, NL, vm, an, e, f, pgspec, plspec, qlspec]
