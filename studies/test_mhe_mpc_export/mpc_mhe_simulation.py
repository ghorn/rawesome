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

#from rawe.ocp.ocprtTools import OcpTools


dae = rawe.dae.Dae()
[x,v] = dae.addX(['x','v'])
[u] = dae.addU(['u'])


dae.setResidual([dae.ddt("x")-v,dae.ddt("v")-u])

# Simulation parameters
N_mpc = 10  # Number of MPC control intervals
N_mhe = 10  # Number of MHE control intervals
Ts = 0.5    # Sampling time
Tf = 5.    # Simulation duration
dt_sim = 0.01  # Simulation integrator step size

# Create the MPC class
mpcRT = makeNmpc(dae,N=N_mpc,dt=Ts)


mpcLog = rawe.ocp.ocprt.Logger(mpcRT,dae)

# Create the MHE class
#mheRT = makeMhe(dae,N=N_mhe,dt=Ts)

# Initialize the MPC-MHE scheme
#mpcRT.initialize()
#mheRT.initialize()

# Create a simulation class
sim = rawe.sim.Sim(dae,dt_sim)

mpcRT.x0 = np.array([1,0])
    
time = 0
while time < Tf:
#    print time
    
#    mheRT.preparationStep()
#    fbret = mheRT.feedbackStep()
#    if fbret != 0:
#        raise Exception("MHE feedbackStep returned error code "+str(fbret))
#    
#    mheRT.log(sim)

    if time > 0:    
        mpcRT.x0 = mpcRT.x[1,:]
    
    mpcRT.preparationStep()
    fbret = mpcRT.feedbackStep()
    if fbret != 0:
        raise Exception("MPC feedbackStep returned error code "+str(fbret))
    
#    mpcRT.log(sim)
    mpcLog.log(mpcRT)
    
    new_x = sim.step(mpcRT.x[0,:],mpcRT.u[0,:],{})    
#    outs = sim.getOutputs(mpcRT.x[0,:],mpcRT.u[0,:],{})

    # Shift MPC
    mpcRT.shift()
    
    # Shift MHE
#    mheRT.shift(new_x=new_state,new_u=new_control,new_y=new_reference,new_yN=new_referenceN,new_S=new_matrices,new_SN=new_matricesN)

    time += Ts
    
plt.ion()
mpcLog.plot(['x','v'])
mpcLog.subplot([['x'],['v']])
mpcLog.plot('u')
plt.show()