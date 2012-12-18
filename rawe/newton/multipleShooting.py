import numpy as np

import casadi as C

import msmaps
from ocputils import Constraints

class Nmpc(object):
    def __init__(self,dae,nk):
        self.dae = dae
        self.nk = nk

        self._gaussNewtonObjF = []

        X = C.ssym('x',dae.xVec().size(),self.nk+1)
        U = C.ssym('u',dae.zVec().size(),self.nk)
        p = C.ssym('p',dae.pVec().size())
        self._dvMap = msmaps.NmpcMap(self.dae,self.nk,X,U,p)

        X = np.resize(np.array([None]),(dae.xVec().size(),self.nk+1))
        U = np.resize(np.array([None]),(dae.zVec().size(),self.nk))
        p = np.resize(np.array([None]),dae.pVec().size())
        self._boundMap = msmaps.NmpcMap(self.dae,self.nk,X,U,p)

        X = np.resize(np.array([None]),(dae.xVec().size(),self.nk+1))
        U = np.resize(np.array([None]),(dae.zVec().size(),self.nk))
        p = np.resize(np.array([None]),dae.pVec().size())
        self._guessMap = msmaps.NmpcMap(self.dae,self.nk,X,U,p)

        self._outputMapGenerator = msmaps.OutputMapGenerator(self)
        self._outputMap = msmaps.OutputMap(self._outputMapGenerator, self._dvMap.vectorize())

        self._constraints = Constraints()
        
    def __call__(self,*args,**kwargs):
        return self.lookup(*args,**kwargs)

    def lookup(self,name,timestep=None):
        try:
            return self._dvMap.lookup(name,timestep=timestep)
        except NameError:
            pass
        try:
            return self._outputMap.lookup(name,timestep)
        except NameError:
            pass
        raise NameError("unrecognized name \""+name+"\"")
        

    def bound(self,name,(lb,ub),timestep=None):
        self._boundMap.setVal(name,(lb,ub),timestep=timestep)

    def constrain(self,lhs,comparison,rhs,tag='unnamed_constraint'):
        self._constraints.add(lhs,comparison,rhs,tag)
        
    def setObj(self,obj):
        if hasattr(self,'_obj'):
            raise ValueError("don't change the objective function")
        self._obj = obj

    def addGaussNewtonObjF(self,gnF):
        self._gaussNewtonObjF.append(gnF)

