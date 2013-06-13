# Copyright 2012-2013 Greg Horn
#
# This file is part of rawesome.
#
# rawesome is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rawesome is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with rawesome.  If not, see <http://www.gnu.org/licenses/>.

import time
import casadi as C
import numpy
import matplotlib.pyplot as plt

def getitemMsg(d,name,msg):
    try:
        return d[name]
    except KeyError:
        raise Exception('you forgot to put "'+name+'" in '+msg)

def vectorizeXUP(x,u,p,dae):
    if type(x) == dict:
        xVec = C.veccat([getitemMsg(x,name,"the states") for name in dae.xNames()])
    else:
        xVec = x
    if type(u) == dict:
        uVec = C.veccat([getitemMsg(u,name,"the controls") for name in dae.uNames()])
    else:
        uVec = u
    if type(p) == dict:
        pVec = C.veccat([getitemMsg(p,name,"the parameters") for name in dae.pNames()])
    else:
        pVec = p
    return (xVec,uVec,pVec)

def maybeToScalar(x):
    if x.size() == 1:
        return x.at(0)
    else:
        return x

class Timer(object):
    def __init__(self,dt):
        self.dt = dt

    def start(self):
        self._t0 = time.time()
        self.nextTime = self._t0 + self.dt

    def sleep(self):
        time_now = time.time()
        tToWait = self.nextTime - time_now
        if tToWait > 0:
            time.sleep(tToWait)
        else:
            self.nextTime = time_now
        self.nextTime = self.nextTime + self.dt

    def get(self):
        return self.nextTime - self.dt - self._t0

class Sim(object):
    def __init__(self, dae, ts):
        print "creating integrator"
        self.dae = dae
        self._ts = ts
        self.integrator = C.IdasIntegrator(self.dae.casadiDae())
        self.integrator.setOption("reltol",1e-6)
        self.integrator.setOption("abstol",1e-8)
        self.integrator.setOption("t0",0)
        self.integrator.setOption("tf",ts)
        self.integrator.setOption('name','integrator')
        self.integrator.setOption("linear_solver",C.CSparse)
        self.integrator.init()

        print "creating outputs function"
        (fAll, (f0,outputs0names)) = self.dae.outputsFun()
        self.outputsFunAll = fAll
        self.outputsFun0 = f0
        self.outputs0names = outputs0names
        
        self.xNames = dae.xNames()
        self.uNames = dae.uNames()
        self.outputNames = dae.outputNames()
#        self.uNames = dae.uNames()
        listOut=[]
        for n in self.outputNames: listOut.append([])
        self._log = {'x':[],'u':[],'y':[],'yN':[],'outputs':dict(zip(self.outputNames,listOut))}
        
    def step(self, x, u, p):
        (xVec,uVec,pVec) = vectorizeXUP(x,u,p,self.dae)
        self.integrator.setInput(xVec,C.INTEGRATOR_X0)
        self.integrator.setInput(C.veccat([uVec,pVec]),C.INTEGRATOR_P)
        self.integrator.evaluate()
        xNext = C.DMatrix(self.integrator.output())
        if type(x) == dict:
            ret = {}
            for k,name in enumerate(self.dae.xNames()):
                ret[name] = xNext[k].at(0)
        else:
            ret = xNext
        return ret

    def getOutputs(self, x, u, p):
        if self.outputsFun0 == None:
            return {}
        (xVec,uVec,pVec) = vectorizeXUP(x,u,p,self.dae)
        self.outputsFun0.setInput(xVec, 0)
        self.outputsFun0.setInput(uVec, 1)
        self.outputsFun0.setInput(pVec, 2)
        self.outputsFun0.evaluate()
        ret = {}
        for k,name in enumerate(self.outputs0names):
            ret[name] = maybeToScalar(C.DMatrix(self.outputsFun0.output(k)))
        return ret
    
    def log(self,new_x=None,new_u=None,new_y=None,new_yN=None,new_out=None):
        if new_x != None:
            self._log['x'].append(numpy.array(new_x))
        if new_u != None:
            self._log['u'].append(numpy.array(new_u))
        if new_y != None:
            self._log['y'].append(numpy.array(new_y))
        if new_yN != None:
            self._log['yN'].append(numpy.array(new_yN))
        if new_out != None:
            for name in new_out.keys():
                self._log['outputs'][name].append(numpy.array(new_out[name]))
    
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
                ys = numpy.squeeze(self._log['x'])[:,index]
                ts = numpy.arange(len(ys))*self._ts
                plt.plot(ts,ys,style)
                
            # if it's a control
            if name in self.uNames:
                index = self.uNames.index(name)
                ys = numpy.squeeze(self._log['u'])[:,index]
                ts = numpy.arange(len(ys))*self._ts
                plt.step(ts,ys,style)
                
            if name in self.outputNames:
                index = self.outputNames.index(name)
                ys = numpy.squeeze(self._log['outputs'][name])
                ts = numpy.arange(len(ys))*self._ts
                plt.plot(ts,ys,style)

        if title is not None:
            assert isinstance(title,str), "title must be a string"
            plt.title(title)
        plt.xlabel('time [s]')
        if showLegend is True:
            plt.legend(legend)
        plt.grid()
