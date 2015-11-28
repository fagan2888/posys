import numpy as np

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
    from network import get_network
    
    bus,line = get_network()

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
