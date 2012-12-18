import casadi as C
import numpy as np

class Nmpc(object):
    def __init__(self,dae,nk):
        self.dae = dae
        self.nk = nk

        X = C.ssym('x',dae.xVec().size(),self.nk+1)
        U = C.ssym('u',dae.zVec().size(),self.nk)
        p = C.ssym('p',dae.pVec().size())
        self._dvMap = NmpcMap(self.dae,self.nk,X,U,p)

        X = np.resize(np.array([None]),(dae.xVec().size(),self.nk+1))
        U = np.resize(np.array([None]),(dae.zVec().size(),self.nk))
        p = np.resize(np.array([None]),dae.pVec().size())
        self._boundMap = NmpcMap(self.dae,self.nk,X,U,p)

        X = np.resize(np.array([None]),(dae.xVec().size(),self.nk+1))
        U = np.resize(np.array([None]),(dae.zVec().size(),self.nk))
        p = np.resize(np.array([None]),dae.pVec().size())
        self._guessMap = NmpcMap(self.dae,self.nk,X,U,p)

    def __call__(self,*args,**kwargs):
        return self.lookup(*args,**kwargs)

    def lookup(self,name,timestep=None):
        return self._dvMap.lookup(name,timestep=timestep)

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
