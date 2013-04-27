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
dae['measurementsN'] = x

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

from rawe.dae.rienIntegrator import RienIntegrator
Rint = RienIntegrator(dae,ts=Ts, numIntegratorSteps=400, integratorType='INT_IRK_GL2')
#Rint.getOutputs()
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
    
err1 = []
err2 = []

PL = np.eye(2)
WL = np.eye(2)*1e10
VL = mheRT.S
xL = mpcRT.x0

time = 0
while time < Tf:
    mheRT.preparationStep()
    fbret = mheRT.feedbackStep()
    if fbret != 0:
        raise Exception("MHE feedbackStep returned error code "+str(fbret))
    
    mpcRT.x0 = np.squeeze(mheRT.x[-1,:])
    
    mpcRT.preparationStep()
    fbret = mpcRT.feedbackStep()
    if fbret != 0:
        raise Exception("MPC feedbackStep returned error code "+str(fbret))
    
    yL = mheRT.y[0,:]
    x = mheRT.x[0,:]
    u = mheRT.u[0,:]
    PL1, xL1 = UpdateArrivalCost(Rint, x, u, xL, yL, PL, VL, WL)
    print xL1
    PL = PL1
    xL = xL1
    
    SimulateAndShift(mpcRT,mheRT,sim,mpcLog,mheLog,simLog)
    
    time += Ts
    print time
    


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
