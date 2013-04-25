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
from MHE import makeMhe
from mpc_mhe_utils import *

#from rawe.ocp.ocprtTools import OcpTools

Q = np.eye(2)
Q[0,0] = 2
R = np.eye(1)
N = np.zeros((2,1))
A = np.array([[0,1],[0,0]])
B = np.array([[0,1]]).T
A = np.array([[1,0.1],[0,1]])
B = np.array([[0.005,0.1]]).T
#A = np.array([[0,0],[0,0]])
#B = np.array([[0,0]]).T

dae = rawe.dae.Dae()
[x,v] = dae.addX(['x','v'])
[u] = dae.addU(['u'])

dae['measurements'] = C.vertcat([x,u])
dae['measurementsN'] = C.vertcat([x])

dae.setResidual([dae.ddt("x")-v,
                 dae.ddt("v")-u])

# Simulation parameters
N_mpc = 10  # Number of MPC control intervals
N_mhe = 10  # Number of MHE control intervals
Ts = 0.1    # Sampling time
Tf = 5.    # Simulation duration

# Create the MPC class
mpcRT, intOpts = makeNmpc(dae,N=N_mpc,dt=Ts)
mheRT, _ = makeMhe(dae,N=N_mpc,dt=Ts)

mpcLog = InitializeMPC(mpcRT,dae)
mheLog = InitializeMHE(mheRT,dae)



# Create the MHE class
#mheRT = makeMhe(dae,N=N_mhe,dt=Ts)

# Initialize the MPC-MHE scheme
#mpcRT.initialize()
#mheRT.initialize()

# Create a simulation class
#intOptions = {'type':'Rintegrator', 'numIntegratorSteps':400, 'integratorType':'INT_IRK_GL2', 'ts':Ts}
intOptions = {'type':'Idas', 'ts':Ts}

sim, simLog = InitializeSim(dae,intOptions)

mpcRT.x0 = np.array([1,0])
outs = sim.getOutputs(mpcRT.x[0,:],mpcRT.u[0,:],{})
new_y  = np.squeeze(outs['measurements'])
simLog.log(mpcRT.x0,new_y,[])
    
time = 0
while time < Tf:
#    print time
    
    mheRT.preparationStep()
    fbret = mheRT.feedbackStep()
    if fbret != 0:
        raise Exception("MHE feedbackStep returned error code "+str(fbret))
    
    mheLog.log(mheRT)
    
    # Get the latest measurement
    mpcRT.x0 = mheRT.x[-1,:]

    mpcRT.preparationStep()
    fbret = mpcRT.feedbackStep()
    if fbret != 0:
        raise Exception("MPC feedbackStep returned error code "+str(fbret))
    
    mpcLog.log(mpcRT)
    
    # Get the measurement BEFORE simulating
    outs = sim.getOutputs(mpcRT.x[0,:],mpcRT.u[0,:],{})
    new_y  = np.squeeze(outs['measurements'])
#    new_y = mhe.makeMeasurements(x,u)
    # Simulate the system
    new_x = sim.step(mpcRT.x[0,:],mpcRT.u[0,:],{})  
    # Get the last measurement AFTER simulating
    outs = sim.getOutputs(new_x,mpcRT.u[0,:],{})
    new_yN = np.array([outs['measurementsN']])
    
    simLog.log(new_x,new_y,new_yN)
    
    # Assign the new initial value and shift
    mpcRT.x0 = np.array(new_x)[:,0]
    mpcRT.shift()
    mheRT.shift(new_y=new_y,new_yN=new_yN)
    
    time += Ts
    


plt.ion()

Fig_plot(['x','v'],what=['sim','mhe','mpc'],simLog=simLog,mheLog=mheLog,mpcLog=mpcLog)

mpcLog.plot(['x','v'])
mpcLog.subplot([['x'],['v']])
mpcLog.plot('u')

mpcLog.plot(['x','v'],when='all')

plt.figure()
for k in range(1,N_mpc):
    plt.plot(np.array(mpcLog._log['x'])[k,:,0],np.array(mpcLog._log['x'])[k,:,1])
plt.grid()

plt.show()
