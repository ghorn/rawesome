# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 15:19:57 2013

@author: mzanon
"""
import numpy as np
import scipy
import copy
import matplotlib.pyplot as plt

import casadi as C
import rawe


def dlqr(A, B, Q, R, N=None):
    
    if N == None:
        N = np.zeros((Q.shape[0],R.shape[0]))
    
    P = scipy.linalg.solve_discrete_are(A, B, Q, R)
    
#    nx = A.shape[0]
#    nu = B.shape[1]
#    PL = np.eye(nx)
#    for k in range(10):
#        M = np.bmat([[             Q, np.zeros((nx,nu)), np.zeros((nx,nx)) ],
#                     [             0,                 R, np.zeros((nu,nx)) ],
#                     [ -np.dot(PL,A),     -np.dot(PL,B),                PL ]])
#        np.linalg.inv(np.dot(M.T,M))*M.T
        
    
    k1 = np.dot(np.dot(B.T, P), B) + R
    k2 = np.dot(np.dot(B.T, P), A)
    K = np.linalg.solve(k1,k2)
    
    return K, P
    
def ComputeTerminalCost(integrator, xlin, ulin, Q, R, N=None): 
    
    integrator.x = xlin
    integrator.u = ulin
    integrator.step()
    A = integrator.dx1_dx0
    B = integrator.dx1_du
#    integrator.getOutputs()
    
    K, P = dlqr(A, B, Q, R, N=None)
    
    return K, P, A, B
    
def UpdateArrivalCost(integrator, x, u, xL, yL, PL, VL, WL): 
    ''' Arrival cost implementation.
        Approximate the solution of:
        min_{xL_,uL_,xL1_} ||  PL ( xL_-xL )         ||
                           ||  VL ( yL-h(xL_,uL_) )  ||
                           ||  WL wx                 ||
                     s.t.  wx = xL1_ - f(xL_,uL_)
        
        Linearization (at the last MHE estimate x,u which is different from xL,uL):
        f(xL_,uL_) ~= f(x,u) + df(x,u)/dx (xL_-x) + df(x,u)/du (uL_-u)
                   ~= f(x,u) +         Xx (xL_-x) +         Xu (uL_-u)
                   ~= f(x,u) - Xx x - Xu u + Xx xL_ + Xu uL_
                   ~= x_tilde              + Xx xL_ + Xu uL_
        h(xL_,uL_) ~= h(x,u) + dh(x,u)/dx (xL_-x) + dh(x,u)/du (uL_-u)
                   ~= f(x,u) +         Hx (xL_-x) +         Hu (uL_-u)
                   ~= h(x,u) - Hx x - Hu u + Hx xL_ + Hu uL_
                   ~= h_tilde              + Hx xL_ + Hu uL_
                   
        Linearized problem:
        min_{xL_,uL_,xL1_} ||  PL ( xL_ - xL )                          ||
                           ||  VL ( yL - h_tilde - Hx xL_ - Hu uL_ )    ||
                           ||  WL ( xL1_ - x_tilde - Xx xL_ - Xu uL_ )  ||
        
        Rewrite as:
        min_{xL_,uL_,xL1_} ||  M ( xL_, uL_, xL1_ ) + res  ||
        
        After QR factorization of M:
        min_{xL_,uL_,xL1_} ||  R ( xL_, uL_, xL1_ ) + rho  ||
            
        '''
    nx = x.shape[0]
    nu = u.shape[0]
    nV = VL.shape[0]
    
    integrator.x = x
    integrator.u = u
    out = integrator.getOutputs()
    h = np.squeeze(out['measurements'])
    x1 = integrator.step()
    Xx = integrator.dx1_dx0
    Xu = integrator.dx1_du
    
    Hx = np.diag([1,0])
    Hu = np.array([[0,1]]).T
    
    x_tilde = x1 - np.dot(Xx,x) - np.dot(Xu,u)
    h_tilde =  h - np.dot(Hx,x) - np.dot(Hu,u)
    
    res = np.bmat([ -np.dot(PL, xL),
                     np.dot(VL, yL - h_tilde),
                    -np.dot(WL, x_tilde) ])
    res = np.squeeze(np.array(res))
    
    M = np.bmat([[             PL,  np.zeros((nx,nu)), np.zeros((nx,nx)) ],
                 [ -np.dot(VL,Hx),     -np.dot(VL,Hu), np.zeros((nV,nx)) ],
                 [ -np.dot(WL,Xx),     -np.dot(WL,Xu),                WL ]])
    
    Q, R = np.linalg.qr(M)
    
#    R1  = R[:nx+nu,:nx+nu]
#    R12 = R[:nx+nu,nx+nu:]
    R2  = R[nx+nu:,nx+nu:]
    
#    rho = np.linalg.solve(Q,res)
    rho = np.squeeze(np.array(np.dot(Q.T,res)))
    rho2 = rho[nx+nu:]
    
    PL1 = R2
    xL1 = -np.linalg.solve(R2,rho2)
    
    return np.array(PL1), np.array(xL1)
    
def GenerateReference(dae,conf,refP):
    
    z0 = refP['z0']
    r0 = refP['r0']
    ddelta0 = refP['ddelta0']
    steadyState, _ = getSteadyState(dae,conf,ddelta0,r0,z0)
    
    xref = {}
    uref = {}
    for name in dae.xNames():
        xref[name] = steadyState[name]
    for name in dae.uNames():
        uref[name] = steadyState[name]
    
    return xref, uref

def InitializeMPC(mpcrt,dae):
    
    mpcrt.x = np.zeros(mpcrt.x.shape)
    mpcrt.u = np.zeros(mpcrt.u.shape) 
    
    mpcrt.S[0,0] = 2.0#/(xRms)**2
    mpcrt.S[1,1] = 1.0#/(vRms)**2
    mpcrt.S[2,2] = 1.0#/(fRms)**2

    mpcrt.SN[0,0] = 1.0#/(xRms)**2
    mpcrt.SN[1,1] = 1.0#/(vRms)**2
    P = np.array([[2.870184272463204e+01,     1.415980225850630e+01], 
                  [1.415980225850630e+01,     2.011263075885158e+01]])# 28.7018-2,   14.1598], [14.1598,   20.1126-1]])
    mpcrt.SN = P
    
#    LinearizeSystem(mpcrt,dae)

#    mpcLog = rawe.ocp.ocprt.Logger(mpcrt,dae)
#    
#    return mpcLog
   
def InitializeMHE(mhert,dae):
    
    mhert.x = np.zeros(mhert.x.shape)
    mhert.u = np.zeros(mhert.u.shape) 
    mhert.y = np.zeros(mhert.y.shape)
    mhert.y[:,0] += 1
    mhert.yN = np.zeros(mhert.yN.shape)+1 
    
    mhert.S[0,0] = 1.0#/(xRms)**2
    mhert.S[1,1] = 1.0#/(vRms)**2

    mhert.SN[0,0] = 1.0#/(xRms)**2
    
#    mheLog = rawe.ocp.ocprt.Logger(mhert,dae)
    
#    return mheLog
    
def SimulateAndShift(mpcRT,mheRT,sim,simLog):

    mheRT.log()    
    mpcRT.log()
    
    # Get the measurement BEFORE simulating
    outs = sim.getOutputs(mpcRT.x[0,:],mpcRT.u[0,:],{})
    new_y  = np.squeeze(outs['measurements']) + np.random.randn(2)*0.01
    # Simulate the system
    new_x = sim.step(mpcRT.x[0,:],mpcRT.u[0,:],{})
    # Get the last measurement AFTER simulating
    outs = sim.getOutputs(new_x,mpcRT.u[0,:],{})
#    outsR = Rint.getOutputs(x=np.squeeze(new_x),u=mpcRT.u[0,:])
    new_yN = np.array([outs['measurementsN']]) + np.random.randn(1)*0.01
    
    simLog.log(new_x=new_x,new_y=new_y,new_yN=new_yN,new_out=outs)
    
    # shift
    mpcRT.shift()
    mheRT.shift(new_y=new_y,new_yN=new_yN)
    
class SimLog(object):
    def __init__(self,dae,sim):
        self.xNames = dae.xNames()
        self.outputNames = dae.outputNames()
#        self.uNames = dae.uNames()
        self.Ts = sim._ts
        l=[]
        for n in self.outputNames: l.append([])
        self._log = {'x':[],'y':[],'yN':[],'outputs':dict(zip(self.outputNames,l))}
        
#        self.log()
    
    def log(self,new_x=None,new_y=None,new_yN=None,new_out=None):
        if new_x != None:
            self._log['x'].append(np.array(new_x))
        if new_y != None:
            self._log['y'].append(np.array(new_y))
        if new_yN != None:
            self._log['yN'].append(np.array(new_yN))
        if new_out != None:
            for name in new_out.keys():
                self._log['outputs'][name].append(np.array(new_out[name]))
        
    def _plot(self,names,title,style,when=0,showLegend=True):
        if isinstance(names,str):
            names = [names]
        assert isinstance(names,list)

        legend = []
        for name in names:
            assert isinstance(name,str)
            legend.append(name)

            # if it's a differential state
            if name in self.xNames:
                index = self.xNames.index(name)
                ys = np.squeeze(self._log['x'])[:,index]
                ts = np.arange(len(ys))*self.Ts
                plt.plot(ts,ys,style)
                
            if name in self.outputNames:
                index = self.outputNames.index(name)
                ys = np.squeeze(self._log['outputs'][name])
                ts = np.arange(len(ys))*self.Ts
                plt.plot(ts,ys,style)

        if title is not None:
            assert isinstance(title,str), "title must be a string"
            plt.title(title)
        plt.xlabel('time [s]')
        if showLegend is True:
            plt.legend(legend)
        plt.grid()
    
def InitializeSim(dae,intOptions):
    
    Ts = intOptions['ts']
    
    if intOptions['type'] == 'Idas':
        sim = rawe.sim.Sim(dae,Ts)
    elif intOptions['type'] == 'Rintegrator':
        from rawe.dae.rienIntegrator import RienIntegrator
        nSteps = intOptions['numIntegratorSteps']
        Type = intOptions['integratorType']
        sim = RienIntegrator(dae,ts=Ts, numIntegratorSteps=nSteps, integratorType=Type)
    else:
        raise Exception('integrator not supported')
    
    simLog = SimLog(dae,sim)
    
    return sim, simLog
    
def Fig_plot(names,title=None,style='',when=0,showLegend=True,what=[],mpcLog=None,mheLog=None,simLog=None):
    assert isinstance(what,list)
    
    fig = plt.figure()
    
    if title is None:
        if isinstance(names,str):
            title = names
        else:
            assert isinstance(names,list)
            if len(names) == 1:
                title = names[0]
            else:
                title = str(names)
    fig.canvas.set_window_title(str(title))

    plt.clf()
    
    if 'mpc' in what:
        if mpcLog == None: raise Exception('you must provide a mpc log to plot its variables')
        mpcLog._plot(names,None,'k',when='all',showLegend=True)
    if 'sim' in what:
        if simLog == None: raise Exception('you must provide a sim log to plot its variables')
        simLog._plot(names,None,'',when=0,showLegend=True)
    if 'mhe' in what:
        if mheLog == None: raise Exception('you must provide a mhe log to plot its variables')
        N = mheLog._log['x'][0].shape[0]
        if not isinstance(names,list):
            names = [names]
        if names[0] in mheLog.xNames:
            mheLog._plot(names,None,'o',when=N-1,showLegend=True)
        elif names[0] in mheLog.uNames:
            mheLog._plot(names,None,'o',when=N-2,showLegend=True)
        
def Fig_subplot(names,title=None,style='',when=0,showLegend=True,what=[],mpcLog=None,mheLog=None,simLog=None):
    assert isinstance(what,list)
    assert isinstance(names,list)
    
    fig = plt.figure()
    
    if title is None:
        if isinstance(names,str):
            title = names
        else:
            assert isinstance(names,list)
            if len(names) == 1:
                title = names[0]
            else:
                title = str(names)
    fig.canvas.set_window_title(str(title))

    plt.clf()
    
    n = len(names)
    for k,name in enumerate(names):
        plt.subplot(n,1,k+1)
        if 'mpc' in what:
            if mpcLog == None: raise Exception('you must provide a mpc log to plot its variables')
            mpcLog._plot(name,None,'k',when='all',showLegend=True)
        if 'sim' in what:
            if simLog == None: raise Exception('you must provide a sim log to plot its variables')
            simLog._plot(name,None,'',when=0,showLegend=True)
        if 'mhe' in what:
            if mheLog == None: raise Exception('you must provide a mhe log to plot its variables')
            N = mheLog._log['x'][0].shape[0]
            if not isinstance(name,list):
                name = [name]
            if name[0] in mheLog.xNames:
                mheLog._plot(name,None,'o',when=N-1,showLegend=True)
            elif name[0] in mheLog.uNames:
                mheLog._plot(name,None,'o',when=N-2,showLegend=True)
        