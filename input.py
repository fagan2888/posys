# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 17:47:04 2015

@author: Sara Arango and Gregory Dobler
"""
import numpy as np
# Sample input for a very simple network

GG = np.array([[4,0,-4],
                 [0,4,-4],
                [-4,-4,8]])
BB = np.array([[-5,0,5], 
                 [0,-10,10],
                [5,10,-15]])

#P = np.transpose(np.array([0, 1.70,-2]))
#Q = np.transpose(np.array([0, 0, -1]))
Q_d = np.transpose(np.array([0, 0, -0.9397]))
P_d = np.transpose(np.array([0, 2.6860,-2.9313]))
#Q = np.transpose(np.array([0, 1.5365, -0.9397])) # Q that makes the jacobian work
#P = np.transpose(np.array([0, 0,-2.9313])) # P that makes the jacobian work
N = 3
NG = 2
NL = 1
vm = np.array([0,1.1249,0.93834])
an = np.array([0, 6.3, -3.44])
e = np.transpose(vm*np.cos(np.radians(an)))
f = np.transpose(vm*np.sin(np.radians(an)))
pgspec = np.array([1.7])
plspec = np.array([-2])
qlspec = np.array([-1])

# Single letter naming