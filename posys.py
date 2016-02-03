#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from pdb import set_trace as stop
from scipy import sparse
import matplotlib.pyplot as plt

"""

Functions:

- jaco: Assembles Jacobian .
- syspower: Calculates real and reactive power P and Q.
- ybus: 
- getnetwork
- get_params

"""

def jaco(GG,BB,P,Q,N,NG,NL,vm,e,f,ind_gen,ind_load):
    """
    Calculate the Jacobian
    """
    # -- utilities
    iden = np.eye(len(f))
    inds = iden>0
    a    = GG*e - BB*f
    b    = GG*f + BB*e
    ep   = iden*e
    fp   = iden*f
    ep_a = np.dot(ep,a)
    ep_b = np.dot(ep,b)
    fp_a = np.dot(fp,a)
    fp_b = np.dot(fp,b)

    # Form H, N, J and L matrices with standard sizes
    HH       =  fp_a - ep_b
    NN       =  ep_a + fp_b
    JJ       = -ep_a - fp_b
    LL       =  fp_a - ep_b
    HH[inds] = -Q - BB[inds]*vm**2 ## CAREFUL
    NN[inds] =  P + GG[inds]*vm**2
    JJ[inds] =  P - GG[inds]*vm**2
    LL[inds] =  Q - BB[inds]*vm**2

    return np.vstack((np.hstack((HH[1:NG+NL,1:NG+NL],NN[1:NG+NL,ind_load])),
                      np.hstack((JJ[ind_load,1:NG+NL], LL[ind_load][:,ind_load]))))



def syspower(GG,BB,N,NG,NL,e,f,pgspec,plspec,qlspec,ind_gen,ind_load):
    """
    Calculate the real and reactive power (P and Q) as well as Delta P and
    Delta Q
    """
    # -- calculate P and Q
    a_sum = (GG*e - BB*f).sum(1)
    b_sum = (GG*f + BB*e).sum(1)

    P     = e*a_sum + f*b_sum
    Q     = f*a_sum - e*b_sum

    # -- calculate the Deltas
    dpq   = np.concatenate([pgspec[:NG-1] - P[ind_gen],
                            plspec[:NL] - P[ind_load],
                            qlspec[:NL] - Q[ind_load]])

    return P, Q, dpq



def ybus(bus,line):
    """
    Return the GG and BB matrices from input bus and line info.

    This function computes the Ybus Matrix necessary
    for the Newtons load flow programs
    The imput comes from two matrices 'line' and 'bus'

     The input matrix 'line' is organized in 7 columns:

      (1)   (2)   (3)   (4)   (5)   (6)   (7)
       NS    NR    R     XL    BC    TR    TS

     NS = Sending node  (iterger number)
     NR = Receiving node  (integer number)
     R  = Resistance
     XL = Series reactance
     BC = Shunt susceptance
     TR = Transformers tap ratio
     TS = Transformers phase shift

    The input matrix 'bus' is organized in 10 columns:

      (1)   (2)   (3)   (4)   (5)   (6)   (7)   (8)   (9)   (10)
       Bn   Vol   Ang    Pg    Qg    Pl    Ql    Bg    Bb     Bt

     Bn  = Bus number
     Vol = Voltage
     Ang = Angle
     Pg  = Active Power Generated
     Qg  = Reactive Power Generated
     Pl  = Active Power Consumed
     Ql  = Reactive Power Consumed
     Bg  = Bus conductance
     Bb  = Bus susceptance
     BT  = Bus type ( 1 = slack , 2 = generator , 3 = load )

    The ouptputs are two Matices GG and BB components YBUS
    YBUS = GG + j BB
    """

    # -- extract line info
    ns, nr, r, xl, bc, tr, ts = (line.T)[:7]

    # -- utilities
    jimag = 1j
    nline = line.shape[0]
    nbus  = bus.shape[0]

    # -- series impedance & admittance
    z = r + jimag*xl
    y = 1. / z

    # -- initialize an empty sparse Ybus
    # GGD: Francisco used a sparse matrix here!!!
    Y = np.zeros([nbus,nbus],dtype=np.complex)

    # -- buid Y bus
    for ii in range(nline):
        tr[ii] = 1 if tr[ii]==0 else tr[ii]
        tps    = tr[ii]*np.exp(jimag*ts[ii]*np.pi/180)
        send   = ns[ii]-1
        rece   = nr[ii]-1
        Y[send,rece] -= y[ii]/tps.conjugate() # Francisco
        Y[rece,send] -= y[ii]/tps # Francisco
        #Y[send,rece] -= y[ii]/tps;
        #Y[rece,send] = Y[send,rece]
        Y[send,send] += (y[ii]+jimag*bc[ii]/2.)/(tps*tps.conjugate()) # Francisco
        #Y[send,send] += y[ii]/(tps**2) + bc[ii];
        Y[rece,rece] += y[ii]+jimag*bc[ii]/2. # Francisco
        #Y[send,send] += + y[ii] + bc[ii];

    # -- include the bus to ground conductance and susceptance
    # -- bus[:,7] - bus conductance
    # -- bus[:,8] - bus susceptance
    # GGD: not sure if this is translated correctly!!!
    if bus[:,7:9].any():
        Y[np.eye(Y.shape[0])>0] += bus[:,7]+jimag*bus[:,8]

    # -- separate real and imaginary parts since my power flow
    # -- programs compute the Jacobian and powers using them
    return Y.real, Y.imag



def get_network(systype):
    """
    Set the line and bus matrices for the system.
    """
        
    # Case: 14 bus. 3 bus does not need this because it comes assembled.
        
    #      Sending, receiving, R,       X,       B,    tap ratio, phase shift.

#     -- original line matrix from 
#    line = np.array([[1,   2, 0.01938, 0.05917, 0.0000, 0.0000, 0.0],
#                     [1,   5, 0.05403, 0.22304, 0.0000, 0.0000, 0.0],
#                     [2,   3, 0.04699, 0.19797, 0.0000, 0.0000, 0.0],
#                     [2,   4, 0.05811, 0.17632, 0.0000, 0.0000, 0.0],
#                     [2,   5, 0.05695, 0.17388, 0.0000, 0.0000, 0.0],
#                     [3,   4, 0.06701, 0.17103, 0.0000, 0.0000, 0.0],
#                     [4,   5, 0.01335, 0.04211, 0.0000, 0.0000, 0.0],
#                     [4,   7, 0.00000, 0.20912, 0.0000, 1.0000, 0.0],
#                     [4,   9, 0.00000, 0.55618, 0.0000, 1.0000, 0.0],
#                     [5,   6, 0.00000, 0.25202, 0.0000, 1.0000, 0.0],
#                     [6,  11, 0.09498, 0.19890, 0.0000, 0.0000, 0.0],
#                     [6,  12, 0.12291, 0.25581, 0.0000, 0.0000, 0.0],
#                     [6,  13, 0.06615, 0.13027, 0.0000, 0.0000, 0.0],
#                     [7,   8, 0.00000, 0.17615, 0.0000, 0.0000, 0.0],
#                     [7,   9, 0.00000, 0.11001, 0.0000, 0.0000, 0.0],
#                     [9,  10, 0.03181, 0.0845,  0.0000, 0.0000, 0.0],
#                     [9,  14, 0.12711, 0.27038, 0.0000, 0.0000, 0.0],
#                     [10, 11, 0.08205, 0.19207, 0.0000, 0.0000, 0.0],
#                     [12, 13, 0.22092, 0.19988, 0.0000, 0.0000, 0.0],
#                     [13, 14, 0.17093, 0.34802, 0.0000, 0.0000, 0.0]])

    # -- from http://publish.illinois.edu/smartergrid/ieee-14-bus-system/
#    line = np.array([[ 1,  2, 0.01938, 0.05917, 0.0000, 0.0000, 0.0],
#                     [ 1,  5, 0.05403, 0.22304, 0.0000, 0.0000, 0.0],
#                     [ 2,  3, 0.04699, 0.19797, 0.0000, 0.0000, 0.0],
#                     [ 2,  4, 0.05811, 0.17632, 0.0000, 0.0000, 0.0],
#                     [ 2,  5, 0.05695, 0.17388, 0.0000, 0.0000, 0.0],
#                     [ 3,  4, 0.06701, 0.17103, 0.0000, 0.0000, 0.0],
#                     [ 4,  5, 0.01335, 0.04211, 0.0000, 0.0000, 0.0],
#                     [ 6, 11, 0.09498, 0.19890, 0.0000, 0.9780, 0.0],
#                     [ 6, 12, 0.12291, 0.25581, 0.0000, 0.9690, 0.0],
#                     [ 6, 13, 0.06615, 0.13027, 0.0000, 0.9320, 0.0],
#                     [ 7,  8, 0.00000, 0.17615, 0.0000, 0.0000, 0.0],
#                     [ 7,  9, 0.00000, 0.11001, 0.0000, 0.0000, 0.0],
#                     [ 9, 10, 0.03181, 0.08450, 0.0000, 0.0000, 0.0],
#                     [ 9, 14, 0.12711, 0.27038, 0.0000, 0.0000, 0.0],
#                     [10, 11, 0.08205, 0.19207, 0.0000, 0.0000, 0.0],
#                     [12, 13, 0.22092, 0.19988, 0.0000, 0.0000, 0.0],
#                     [13, 14, 0.17093, 0.34802, 0.0000, 0.0000, 0.0]])

#    # -- from PYPOWER
    line = np.array([[ 1,  2, 0.01938, 0.05917, 0.0000, 0.0000, 0.0],
                     [ 1,  5, 0.05403, 0.22304, 0.0000, 0.0000, 0.0],
                     [ 2,  3, 0.04699, 0.19797, 0.0000, 0.0000, 0.0],
                     [ 2,  4, 0.05811, 0.17632, 0.0000, 0.0000, 0.0],
                     [ 2,  5, 0.05695, 0.17388, 0.0000, 0.0000, 0.0],
                     [ 3,  4, 0.06701, 0.17103, 0.0000, 0.0000, 0.0],
                     [ 4,  5, 0.01335, 0.04211, 0.0000, 0.0000, 0.0],
                     [ 4,  7, 0,       0.20912, 0.0000, 0.9780, 0.0],
                     [ 4,  9, 0,       0.55618, 0.0000, 0.9690, 0.0],
                     [ 5,  6, 0,       0.25202, 0.0000, 0.9320, 0.0],
                     [ 6, 11, 0.09498, 0.19890, 0.0000, 0.0000, 0.0],
                     [ 6, 12, 0.12291, 0.25581, 0.0000, 0.0000, 0.0],
                     [ 6, 13, 0.06615, 0.13027, 0.0000, 0.0000, 0.0],
                     [ 7,  8, 0,       0.17615, 0.0000, 0.0000, 0.0],
                     [ 7,  9, 0,       0.11001, 0.0000, 0.0000, 0.0],
                     [ 9, 10, 0.03181, 0.08450, 0.0000, 0.0000, 0.0],
                     [ 9, 14, 0.12711, 0.27038, 0.0000, 0.0000, 0.0],
                     [10, 11, 0.08205, 0.19207, 0.0000, 0.0000, 0.0],
                     [12, 13, 0.22092, 0.19988, 0.0000, 0.0000, 0.0],
                     [13, 14, 0.17093, 0.34802, 0.0000, 0.0000, 0.0]])
                     
                         # -- weird experiment n. 1
#    line = np.array([[ 1,  2, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 1,  5, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 2,  3, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 2,  4, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 2,  5, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 3,  4, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 4,  5, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 4,  7, 1.00000, 1.00000, 0.0000, 1.0000, 0.0],
#                     [ 4,  9, 1.00000, 1.00000, 0.0000, 1.0000, 0.0],
#                     [ 5,  6, 1.00000, 1.00000, 0.0000, 1.0000, 0.0],
#                     [ 6, 11, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 6, 12, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 6, 13, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 7,  8, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 7,  9, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 9, 10, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [ 9, 14, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [10, 11, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [12, 13, 1.00000, 1.00000, 0.0000, 0.0000, 0.0],
#                     [13, 14, 1.00000, 1.00000, 0.0000, 0.0000, 0.0]])
                     
                     # -- weird experiment n. 2
#    line = np.array([[ 1,  2, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 1,  5, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 2,  3, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 2,  4, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 2,  5, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 3,  4, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 4,  5, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
##                     [ 4,  7, 0.50000, 0.50000, 0.0000, 0.9780, 0.0],
##                     [ 4,  9, 0.50000, 0.50000, 0.0000, 0.9690, 0.0],
##                     [ 5,  6, 0.50000, 0.50000, 0.0000, 0.9320, 0.0],
#                     [ 4,  7, 0.50000, 0.50000, 0.0000, 0.0300, 0.0],
#                     [ 4,  9, 0.50000, 0.50000, 0.0000, 0.0300, 0.0],
#                     [ 5,  6, 0.50000, 0.50000, 0.0000, 0.0300, 0.0],
#                     [ 6, 11, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 6, 12, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 6, 13, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 7,  8, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 7,  9, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 9, 10, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [ 9, 14, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [10, 11, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [12, 13, 0.50000, 0.50000, 0.0000, 0.0000, 0.0],
#                     [13, 14, 0.50000, 0.50000, 0.0000, 0.0000, 0.0]])




    #               Bus, Vol,   An,   P,    Q,    PL,     QL, Cond, Susc, Type
     
    # -- test with no power
    bus = np.array([[1,  1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 1],
                    [2,  1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 2],
                    [3,  1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 2],
                    [4,  1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 3],
                    [5,  1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 3],
                    [6,  1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 2],
                    [7,  1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 3],
                    [8,  1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 2],
                    [9,  1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 3],
                    [10, 1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 3],
                    [11, 1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 3],
                    [12, 1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 3],
                    [13, 1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 3],
                    [14, 1.000, 0.0, 0.0, 0.000, 0.000, 0.000, 0.0, 0.0, 3]])


    # -- Francisco's inital guesses
#    bus = np.array([[1,  1.060, 0.0, 0.0,   0.000, 0.000,  0.000, 0.0, 0.0, 1],
#                    [2,  1.045, 0.0, 0.4,   0.424, 0.217,  0.127, 0.0, 0.0, 2],
#                    [3,  1.010, 0.0, 0.0,   0.234, 0.942,  0.190, 0.0, 0.0, 2],
#                    [4,  1.000, 0.0, 0.0,   0.000, 0.478, -0.039, 0.0, 0.0, 3],
#                    [5,  1.000, 0.0, 0.0,   0.000, 0.076,  0.016, 0.0, 0.0, 3],
#                    [6,  1.070, 0.0, 0.0,   0.122, 0.112,  0.075, 0.0, 0.0, 2],
#                    [7,  1.000, 0.0, 0.0,   0.000, 0.000,  0.000, 0.0, 0.0, 3],
#                    [8,  1.090, 0.0, 0.0,   0.174, 0.000,  0.000, 0.0, 0.0, 2],
#                    [9,  1.000, 0.0, 0.0,   0.000, 0.295,  0.166, 0.0, 0.0, 3],
#                    [10, 1.000, 0.0, 0.0,   0.000, 0.090,  0.058, 0.0, 0.0, 3],
#                    [11, 1.000, 0.0, 0.0,   0.000, 0.035,  0.018, 0.0, 0.0, 3],
#                    [12, 1.000, 0.0, 0.0,   0.000, 0.061,  0.016, 0.0, 0.0, 3],
#                    [13, 1.000, 0.0, 0.0,   0.000, 0.135,  0.058, 0.0, 0.0, 3],
#                    [14, 1.000, 0.0, 0.0,   0.000, 0.149,  0.050, 0.0, 0.0, 3]])

    # -- true solution
#    bus = np.array([[1,  1.060,   0.00, 2.324,-0.169, 0.000,  0.000, 0.0, 0.0, 1],
#                    [2,  1.045,  -4.98, 0.4,   0.424, 0.217,  0.127, 0.0, 0.0, 2],
#                    [3,  1.010, -12.72, 0.0,   0.234, 0.942,  0.190, 0.0, 0.0, 2],
#                    [4,  1.019, -10.33, 0.0,   0.000, 0.478, -0.039, 0.0, 0.0, 3],
#                    [5,  1.020,  -8.78, 0.0,   0.000, 0.076,  0.016, 0.0, 0.0, 3],
#                    [6,  1.070, -14.22, 0.0,   0.122, 0.112,  0.075, 0.0, 0.0, 2],
#                    [7,  1.062, -13.37, 0.0,   0.000, 0.000,  0.000, 0.0, 0.0, 3],
#                    [8,  1.090, -13.36, 0.0,   0.174, 0.000,  0.000, 0.0, 0.0, 2],
#                    [9,  1.056, -14.94, 0.0,   0.000, 0.295,  0.166, 0.0, 0.0, 3],
#                    [10, 1.051, -15.10, 0.0,   0.000, 0.090,  0.058, 0.0, 0.0, 3],
#                    [11, 1.057, -14.79, 0.0,   0.000, 0.035,  0.018, 0.0, 0.0, 3],
#                    [12, 1.055, -15.07, 0.0,   0.000, 0.061,  0.016, 0.0, 0.0, 3],
#                    [13, 1.050, -15.16, 0.0,   0.000, 0.135,  0.058, 0.0, 0.0, 3],
#                    [14, 1.036, -16.04, 0.0,   0.000, 0.149,  0.050, 0.0, 0.0, 3]])

    return bus,line


def get_params(systype):
    
    """
    Return the parameters for the system.
    """
    
    if systype=='3bus':
        GG       = np.array([[4.,   0,  -4],
                             [0 ,   4,  -4],
                             [-4,  -4,   8]])
        BB       = np.array([[-5.,   0,   5],
                             [ 0, -10,  10],
                             [ 5,  10, -15]])
        GG       += 1e-15*np.random.randn(*GG.shape)
        BB       += 1e-15*np.random.randn(*GG.shape)
        NG       = 2
        NL       = 1
        N        = NG - 1 + 2*NL
#        vm       = np.array([1.0, 1.1249, 0.93834])
#        an       = np.array([0.0,6.3,-3.44])
        vm       = np.array([1.0, 1.0, 1.0])
        an       = np.array([0.0, 0.0, 0.0])
        e        = vm*np.cos(an*np.pi/180.)
        f        = vm*np.sin(an*np.pi/180.)
        pgspec   = np.array([1.7])
        plspec   = np.array([-2])
        qlspec   = np.array([-1])
        ind_gen  = np.array([False, True, False])
        ind_load = np.array([False, False, True])
        
    elif systype == '14bus':
        bus, line = get_network(systype)
        GG, BB    = ybus(bus,line)
        ind_gen   = np.array(bus[:,9]==2)
        ind_load  = np.array(bus[:,9]==3)
        NG        = sum(ind_gen)+1
        NL        = sum(ind_load)
        pgspec    = bus[ind_gen,3]
        plspec    = bus[ind_load,5]
        qlspec    = bus[ind_load,6]
        N         = NG - 1 + 2*NL
        vm        = bus[:,1]
        an        = bus[:,2]
        e         = vm*np.cos(an*np.pi/180.)
        f         = vm*np.sin(an*np.pi/180.)

    return [GG, BB, N, NG, NL, vm, an, e, f, pgspec, plspec, qlspec, ind_gen,
            ind_load]



#####################
## main for now... ##
#####################

# -- set the system type and get the parameters
systype = '14bus'
params  = get_params(systype)
GG, BB, N, NG, NL, vm, an, e, f, pgspec, plspec, qlspec, ind_gen, ind_load = \
    params

# -- set tolerance params
min_tol  = 0.005
itr      = 0
iter_max = 100
tol      = 1 # Is this correct?! In pypower they calculate this at the input.
delta_magnitude = [0]*(iter_max+1)
dpq_magnitude = [0]*(iter_max+1)
vol_list = []

while (itr <= iter_max and tol > min_tol):
    # Compute powers and power mismatches
    P, Q, dpq = syspower(GG,BB,N,NG,NL,e,f,pgspec,plspec,qlspec,ind_gen,
                         ind_load)
    jac       = jaco(GG,BB,P,Q,N,NG,NL,vm,e,f,ind_gen,ind_load)

    # Solve system of equations
#    dtv = np.linalg.solve(jac, dpq)
    dtv       = np.linalg.lstsq(jac, dpq)[0]
#    dtv = sparse.linalg.spsolve(sparse.csr_matrix(jac), dpq)
#    dtv = np.dot(np.linalg.inv(np.dot(jac,jac.T)),np.dot(dpq,jac.T))
#    stop()
    # Track delta magnitude changes
#    delta_magnitude[itr] = np.linalg.norm(dtv,1)
#    dpq_magnitude[itr]   = np.linalg.norm(dpq,1)

    # Update values of voltages and angles
    vm[ind_load] *= 1 + dtv[NG:NG+NL] # load bus
    an[1:N]      += dtv[0:NL+NG-1]*180/np.pi # loads and generators
#    vm[ind_load]  = np.abs(vm[ind_load])
#    an[1:N]       = an[1:N] % 360.

    # Update vectors e and f
    e = vm*np.cos(an*np.pi/180.)
    f = vm*np.sin(an*np.pi/180.)

    # Refresh
    itr += 1

    print("\niteration: {0}".format(itr))
    print("  power : {0}".format(P))
    print("  reac  : {0}".format(Q))
    print("  dpq   : {0}".format(dpq))
    if systype=='3bus':
        print("  jac   : {0}".format(jac[0]))
        print("        : {0}".format(jac[1]))
        print("        : {0}".format(jac[2]))
    print("  dtv   : {0}".format(dtv))
    print("  vm    : {0}".format(vm))
    print("  an    : {0}".format(an))
    print("  ConNumb    : {0}".format(np.linalg.cond(jac)))
    print("  det   : {0}".format(np.linalg.det(jac)))
    
    
    
    #plt.plot(vm)
