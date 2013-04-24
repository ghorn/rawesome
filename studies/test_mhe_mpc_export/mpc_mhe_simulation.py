# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:48:33 2013

@author: mzanon
"""

import rawe
import numpy as np
#import casadi as C

import matplotlib.pyplot as plt

from NMPC import makeNmpc
#from MHE import makeMhe
from mpc_mhe_utils import *

#from rawe.ocp.ocprtTools import OcpTools

#Q = np.eye(2)
#R = np.eye(1)
#N = np.zeros((2,1))
A = np.array([[0,1],[0,0]])
#B = np.array([[0,1]]).T
#A = np.array([[0,0],[0,0]])
#B = np.array([[0,0]]).T

dae = rawe.dae.Dae()
[x,v] = dae.addX(['x','v'])
[u] = dae.addU(['u'])


dae.setResidual([dae.ddt("x")-v,
                 dae.ddt("v")-u])

# Simulation parameters
N_mpc = 10  # Number of MPC control intervals
N_mhe = 10  # Number of MHE control intervals
Ts = 0.5    # Sampling time
Tf = 5.    # Simulation duration

# Create the MPC class
mpcRT, intOpts = makeNmpc(dae,N=N_mpc,dt=Ts)

InitializeMPC(mpcRT,dae)

mpcLog = rawe.ocp.ocprt.Logger(mpcRT,dae)

# Create the MHE class
#mheRT = makeMhe(dae,N=N_mhe,dt=Ts)

# Initialize the MPC-MHE scheme
#mpcRT.initialize()
#mheRT.initialize()

# Create a simulation class
sim = rawe.sim.Sim(dae,Ts)

from rawe.dae.rienIntegrator import RienIntegrator
integrator = RienIntegrator(dae,ts=Ts, numIntegratorSteps=400, integratorType='INT_IRK_GL2')
    
mpcRT.x0 = np.array([0,1])
    
time = 0
while time < Tf:
#    print time
    
#    mheRT.preparationStep()
#    fbret = mheRT.feedbackStep()
#    if fbret != 0:
#        raise Exception("MHE feedbackStep returned error code "+str(fbret))
#    
#    mheRT.log(sim)

    for k in range(5):
        mpcRT.preparationStep()
        fbret = mpcRT.feedbackStep()
        if fbret != 0:
            raise Exception("MPC feedbackStep returned error code "+str(fbret))
    
#    mpcRT.log(sim)
    mpcLog.log(mpcRT)
    
    new_x = sim.step(mpcRT.x[0,:],mpcRT.u[0,:],{})  
    new_xR = integrator.step(mpcRT.x[0,:],mpcRT.u[0,:],{})
    print mpcRT.x0 - mpcRT.x[0,:]
    print new_x - mpcRT.x[1,:]
    print new_xR - mpcRT.x[1,:]
    
    new_xE = np.dot(scipy.linalg.expm(A*Ts),mpcRT.x0)
    assert(1==0)
#    outs = sim.getOutputs(mpcRT.x[0,:],mpcRT.u[0,:],{})

    # Shift MPC
    mpcRT.x0 = mpcRT.x[1,:]
    mpcRT.shift()
    
    # Shift MHE
#    mheRT.shift(new_x=new_state,new_u=new_control,new_y=new_reference,new_yN=new_referenceN,new_S=new_matrices,new_SN=new_matricesN)

    time += Ts
    
plt.ion()
mpcLog.plot(['x','v'])
mpcLog.subplot([['x'],['v']])
mpcLog.plot('u')
plt.show()
