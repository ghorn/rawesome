import numpy as np

import casadi as C

import msmaps

class Nmpc(object):
    def __init__(self,dae,nk):
        self.dae = dae
        self.nk = nk

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
        
#    def setObj(self,obj):
#        if hasattr(self,'_obj'):
#            raise ValueError("don't change the objective function")
#        self._obj = obj

    def setGaussNewtonObjF(self,gnF):
        if hasattr(self,'_gnF'):
            raise ValueError("don't change the objective function")
        self._gnF = gnF
