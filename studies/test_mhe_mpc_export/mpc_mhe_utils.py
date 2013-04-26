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


def dlqr(A, B, Q, R, N=None):
    
    if N == None:
        N = np.zeros((Q.shape[0],R.shape[0]))
    
    P = scipy.linalg.solve_discrete_are(A, B, Q, R)
    
    k1 = np.dot(np.dot(B.T, P), B) + R
    k2 = np.dot(np.dot(B.T, P), A)
    K = np.linalg.solve(k1,k2)
    
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
    
    _, P = dlqr(A, B, Q, R, N=None)
    
    return P


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

    mpcLog = rawe.ocp.ocprt.Logger(mpcrt,dae)
    
    return mpcLog
   
def InitializeMHE(mhert,dae):
    
    mhert.x = np.zeros(mhert.x.shape)
    mhert.u = np.zeros(mhert.u.shape) 
    mhert.y = np.zeros(mhert.y.shape)
    mhert.y[:,0] += 1
    mhert.yN = np.zeros(mhert.yN.shape)+1 
    
    mhert.S[0,0] = 1.0#/(xRms)**2
    mhert.S[1,1] = 1.0#/(vRms)**2

    mhert.SN[0,0] = 1.0#/(xRms)**2
    
    mheLog = rawe.ocp.ocprt.Logger(mhert,dae)
    
    return mheLog
    
class SimLog(object):
    def __init__(self,dae,sim):
        self.xNames = dae.xNames()
#        self.uNames = dae.uNames()
        self.Ts = sim._ts
        self._log = {'x':[],'y':[],'yN':[]}
        
#        self.log()
    
    def log(self,new_x,new_y,new_yN):
        self._log['x'].append(np.array(new_x))
        self._log['y'].append(np.array(new_y))
        self._log['yN'].append(np.array(new_yN))
        
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
        sim = rawe.dae.rienIntegrator(dae,ts=Ts, numIntegratorSteps=400, integratorType='INT_IRK_GL2')
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
        mheLog._plot(names,None,'o',when=N-1,showLegend=True)
        
        
    