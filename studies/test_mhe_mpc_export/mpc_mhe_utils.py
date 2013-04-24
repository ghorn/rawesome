# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 15:19:57 2013

@author: mzanon
"""
import numpy as np
import scipy
import copy
import casadi as C

def LinearizeSystem(mpcrt,dae):
    raise Exception('not implemented yet')
#    (ode, _) = dae.solveForXDotAndZ()
#    
#    odeValue = []
#    for name in dae.xNames():
#        odeValue.append(ode[name])
#    
#    Ode = C.SXFunction([dae.xVec(),dae.uVec()],[C.veccat(odeValue)])
#    Ode.init()
#    Afcn = Ode.jacobian(0)
#    Afcn.init()
#    Bfcn = Ode.jacobian(1)
#    Bfcn.init()
#    
#    dae.Afcn = Afcn
#    dae.Bfcn = Bfcn

def solve_d_are(A, B, Q, R, N):
    
    P = np.eye(A.shape[0])
    
    APBN = np.dot(np.dot(A.T, P), B) + N
    BPBRinv = scipy.linalg.inv( np.dot(np.dot(B.T, P), B) + R )
    
    tol = 1e-10#1e-6
    dPnorm = 1
    k=0
    while dPnorm >= tol:
        P1 = copy.deepcopy( np.dot(np.dot(A.T, P), A) - np.dot(np.dot(APBN, BPBRinv), APBN.T) + Q)
        dPnorm = scipy.linalg.norm(P-P1)
        P = P1
        if k>1000:
            raise Exception('fixed point iteration for DARE solution exceeded max iters')
    
    return P

def dlqr(A, B, Q, R, N=None):
    
    if N == None:
        N = np.zeros((Q.shape[0],R.shape[0]))
    
#    P = scipy.linalg.solve_discrete_are(A, B, Q, R)
    P = solve_d_are(A, B, Q, R, N)
    
    k1 = np.dot(np.dot(B.T, P), B) + R
    k2 = np.dot(np.dot(B.T, P), A)
    K = np.linalg.solve(k1,k2)
    print np.linalg.inv(k1)*k2 - K
    
    return K, P
    
def ComputeTerminalCost(dae, xlin, ulin, Q, R, N=None): 
    
#    dae.Afcn.init()
#    dae.Bfcn.init()
#    
#    dae.Afcn.setInput(xlin,0)
#    dae.Afcn.setInput(ulin,1)
#    dae.Bfcn.setInput(xlin,0)
#    dae.Bfcn.setInput(ulin,1)
#    
#    A = dae.Afcn.output()
#    B = dae.Bfcn.output()
    
    K, P = dlqr(A, B, Q, R, N=None)


def InitializeMPC(mpcrt,dae):
    
    mpcrt.x = np.zeros(mpcrt.x.shape)
    mpcrt.u = np.zeros(mpcrt.u.shape) 
    
    mpcrt.S[0,0] = 2.0#/(xRms)**2
    mpcrt.S[1,1] = 1.0#/(vRms)**2
    mpcrt.S[2,2] = 1.0#/(fRms)**2

    mpcrt.SN[0,0] = 1.0#/(xRms)**2
    mpcrt.SN[1,1] = 1.0#/(vRms)**2
    
#    LinearizeSystem(mpcrt,dae)
   