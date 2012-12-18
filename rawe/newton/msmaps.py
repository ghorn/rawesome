import numpy as np

import casadi as C

class NmpcMap(object):
    def __init__(self,dae,nk,X,U,p):
        self._nk = nk
        self._xNames = dae.xNames()
        self._uNames = dae.uNames()
        self._pNames = dae.pNames()

        self._X = X
        self._U = U
        self._p = p
        
        self._xIdx = {}
        self._uIdx = {}
        self._pIdx = {}
        for k,name in enumerate(self._xNames):
            self._xIdx[name] = k
        for k,name in enumerate(self._uNames):
            self._uIdx[name] = k
        for k,name in enumerate(self._pNames):
            self._pIdx[name] = k

    def vectorize(self):
        return C.veccat([self._X,self._U,self._p])
        # on exception try
        return np.append(self._X.flatten(),self._U.flatten(),self._p.flatten())
    
    def lookup(self,name,timestep=None):
        if name in self._xIdx:
            assert (timestep != None), "please set timestep"
            assert (timestep <= self._nk), "timestep too large"
            return self._X[self._xIdx[name],timestep]
        elif name in self._uIdx:
            assert (timestep != None), "please set timestep"
            assert (timestep < self._nk), "timestep too large"
            return self._U[self._uIdx[name],timestep]
        elif name in self._pIdx:
            assert (timestep == None), "don't set timestep for parameter"
            return self._p[self._pIdx[name]]
        else:
            raise NameError('unrecognized name "'+name+'"')

    def setVal(self,name,val,timestep=None):
        if name in self._xIdx:
            assert (timestep != None), "please set timestep"
            assert (timestep <= self._nk), "timestep too large"
            self._X[self._xIdx[name],timestep] = val
        elif name in self._uIdx:
            assert (timestep != None), "please set timestep"
            assert (timestep < self._nk), "timestep too large"
            self._U[self._uIdx[name],timestep] = val
        elif name in self._pIdx:
            assert (timestep == None), "don't set timestep for parameter"
            self._p[self._pIdx[name]] = val
        else:
            raise NameError('unrecognized name "'+name+'"')
