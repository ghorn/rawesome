import numpy as np
from collutils import mkCollocationPoints

class CollMap(object):
    def __init__(self,ocp,name,devectorize=None):
        # grab everything needed from the ocp
        self._xNames = ocp.dae.xNames()
        self._zNames = ocp.dae.zNames()
        self._uNames = ocp.dae.uNames()
        self._pNames = ocp.dae.pNames()
        self._nk = ocp.nk
        self._nicp = ocp.nicp
        self._deg = ocp.deg
        self._collPoly = ocp.collPoly

        assert isinstance(name,str)
        self._name = name
        self._initMap()

        if devectorize is not None:
            self._devectorize(devectorize)

    def _devectorize(self,V):
        self.vec = V
        ndiff = len(self._xNames)
        nalg = len(self._zNames)
        nu = len(self._uNames)
        NP = len(self._pNames)
        nx = ndiff + nalg
        
        # Total number of variables
        NXD = self._nicp*self._nk*(self._deg+1)*ndiff # Collocated differential states
        NXA = self._nicp*self._nk*self._deg*nalg      # Collocated algebraic states
        NU = self._nk*nu                  # Parametrized controls
        NXF = ndiff                 # Final state (only the differential states)
        NV = NXD+NXA+NU+NXF+NP
        
        offset = 0
        
        # Get the parameters
        P = V[offset:offset+NP]
        for k,name in enumerate(self._pNames):
#            self._indexMap.setVal(name,k)
#            self._dvMap.setVal(name,V[offset+k])
            self.setVal(name,V[offset+k])
        offset += NP
        
        # Get collocated states and parametrized control
        XD = np.resize(np.array([],dtype=type(V)),(self._nk+1,self._nicp,self._deg+1)) # NB: same name as above
        XA = np.resize(np.array([],dtype=type(V)),(self._nk,self._nicp,self._deg+1)) # NB: same name as above
        U = np.resize(np.array([],dtype=type(V)),self._nk)
        
        for k in range(self._nk):
            # Collocated states
            for i in range(self._nicp):
                for j in range(self._deg+1):
                    # Get the expression for the state vector
                    XD[k][i][j] = V[offset:offset+ndiff]
                    for m,name in enumerate(self._xNames):
#                        self._indexMap.setVal(name,offset+m,timestep=k,nicpIdx=i,degIdx=j)
#                        self._dvMap.setVal(name,V[offset+m],timestep=k,nicpIdx=i,degIdx=j)
                        self.setVal(name,V[offset+m],timestep=k,nicpIdx=i,degIdx=j)
                        
                    if j !=0:
                        XA[k][i][j] = V[offset+ndiff:offset+ndiff+nalg]
                        for m,name in enumerate(self._zNames):
#                            self._indexMap.setVal(name,offset+ndiff+m,timestep=k,nicpIdx=i,degIdx=j)
#                            self._dvMap.setVal(name,V[offset+ndiff+m],timestep=k,nicpIdx=i,degIdx=j)
                            self.setVal(name,V[offset+ndiff+m],timestep=k,nicpIdx=i,degIdx=j)

                    index = (self._deg+1)*(self._nicp*k+i) + j
                    if k==0 and j==0 and i==0:
                        offset += ndiff
                    else:
                        if j!=0:
                            offset += nx
                        else:
                            offset += ndiff
            
            # Parametrized controls
            U[k] = V[offset:offset+nu]
            for m,name in enumerate(self._uNames):
#                self._indexMap.setVal(name,offset+m,timestep=k)
#                self._dvMap.setVal(name,V[offset+m],timestep=k)
                self.setVal(name,V[offset+m],timestep=k)
            offset += nu
        
        # State at end time
        XD[self._nk][0][0] = V[offset:offset+ndiff]
        for m,name in enumerate(self._xNames):
#            self._indexMap.setVal(name,offset+m,timestep=self.nk,degIdx=0,nicpIdx=0)
#            self._dvMap.setVal(name,V[offset+m],timestep=self.nk,degIdx=0,nicpIdx=0)
            self.setVal(name,V[offset+m],timestep=self._nk,degIdx=0,nicpIdx=0)
        
        offset += ndiff
        assert offset==NV

        self._xVec = XD
        self._zVec = XA
        self._uVec = U
        self._pVec = P

    def _initMap(self):
        self._xMap = {}
        self._zMap = {}
        self._uMap = {}
        self._pMap = {}
        for name in self._xNames:
            self._xMap[name] = np.resize(np.array([None]),(self._nk+1,self._nicp,self._deg+1))
        for name in self._zNames:
            self._zMap[name] = np.resize(np.array([None]),(self._nk,self._nicp,self._deg+1))
        for name in self._uNames:
            self._uMap[name] = np.resize(np.array([None]),(self._nk))
        for name in self._pNames:
            self._pMap[name] = None

    def fillInMissing(self,mapName,interpFun):
        assert(isinstance(mapName, str))
        import collutils
        tau_root = mkCollocationPoints(self._collPoly,self._deg)

        # all parameters should be set
        for name in self._pNames:
            val = self.lookup(name)
            if val is None:
                raise ValueError(mapName+" for parameter \""+name+"\" is not set")
        
        # all controls should be set
        for name in self._uNames:
            for k in range(self._nk):
                if self.lookup(name,timestep=k) is None:
                    raise ValueError(mapName+" for control \""+name+"\" is not set at timestep "+str(k))

        # all algebraic variables should be set
        for name in self._zNames:
            for k in range(self._nk):
                for j in range(self._nicp):
                    for d in range(1,self._deg+1):
                        if self.lookup(name,timestep=k,nicpIdx=j,degIdx=d) is None:
                            raise ValueError(mapName+" for algebraic variable \""+name+"\" is not set at timestep "+str(k)+", nicpIdx: "+str(j)+", degIdx: "+str(d))

        # states should all be set at degIdx=0, nicpIdx=0
        # if not set in between, call interpFun
        for name in self._xNames:
            for k in range(self._nk):
                # Collocated states
                val0 = self.lookup(name,timestep=k,nicpIdx=0,degIdx=0)
                val1 = self.lookup(name,timestep=k+1,nicpIdx=0,degIdx=0)
                if val0 is None:
                    raise ValueError(mapName+" for state \""+name+"\" is not set at timestep "+str(k))
                if val1 is None:
                    raise ValueError(mapName+" for state \""+name+"\" is not set at timestep "+str(k+1))
                alpha = 0
                alphaIndex = 0
                for j in range(self._nicp):
                    for d in range(self._deg+1):
                        val = self.lookup(name,timestep=k,nicpIdx=j,degIdx=d)
                        if val is None:
                            tau = j + tau_root[d]/float(self._nicp)
                            self.setVal(name,interpFun(tau,val0,val1),timestep=k,nicpIdx=j,degIdx=d)
                        alphaIndex += 1

    def vectorize(self):
        if all([hasattr(self,attr) for attr in ['_xVec','_zVec','uVec','pVec']]):
            V = [self._pVec()]
            for k in range(self._nk):
                # Collocated states
                for i in range(self._nicp):
                    for j in range(self._deg+1):
                        V.append(self.xVec(k,nicpIdx=i,degIdx=j))
                        
                        if j !=0:
                            V.append(self.zVec(k,nicpIdx=i,degIdx=j))
                
                # Parametrized controls
                V.append(self.uVec(k))
        
            # State at end time
            V.append(self.xVec(self._nk),nicpIdx=0,degIdx=0)

            if isinstance(V[0],np.array):
                return np.concatenate(V)
            else:
                import casadi
                return casadi.veccat(V)
        else:
            V = [self.lookup(name) for name in self._pNames]
            
            for k in range(self._nk):
                # Collocated states
                for i in range(self._nicp):
                    for j in range(self._deg+1):
                        V.extend([self.lookup(name,timestep=k,nicpIdx=i,degIdx=j) for name in self._xNames])
                        
                        if j !=0:
                            V.extend([self.lookup(name,timestep=k,nicpIdx=i,degIdx=j) for name in self._zNames])
                
                # Parametrized controls
                V.extend([self.lookup(name,timestep=k) for name in self._uNames])
        
            # State at end time
            V.extend([self.lookup(name,timestep=self._nk,nicpIdx=0,degIdx=0) for name in self._xNames])
            return V
            
    def xVec(self,timestep,nicpIdx=None,degIdx=None):
        assert hasattr(self,'_xVec')
        if nicpIdx is None:
            nicpIdx = 0
        if degIdx is None:
            degIdx = 0
        return self._xVec[timestep][nicpIdx][degIdx]
    def zVec(self,timestep,nicpIdx=None,degIdx=None):
        assert hasattr(self,'_zVec')
        if nicpIdx is None:
            nicpIdx = 0
        assert (degIdx is not None), "must set degIdx in zVec"
        assert (degIdx != 0), "algebraic variables not defined at tau=0"
        return self._zVec[timestep][nicpIdx][degIdx]
    def uVec(self,timestep):
        assert hasattr(self,'_uVec')
        return self._uVec[timestep]
    def pVec(self):
        assert hasattr(self,'_pVec')
        return self._pVec

    def lookup(self,name,timestep=None,nicpIdx=None,degIdx=None):
        ret = self._lookupOrSet(name,timestep,nicpIdx,degIdx)
        if type(ret) is np.ndarray:
            return float(ret)
        else:
            return ret
        
    def setVal(self,name,val,timestep=None,nicpIdx=None,degIdx=None,quiet=False,force=False):
        self._lookupOrSet(name,timestep,nicpIdx,degIdx,setVal=val,quiet=quiet,force=force)
        
    def _lookupOrSet(self,name,timestep,nicpIdx,degIdx,setVal=None,quiet=False,force=False):
        assert isinstance(name,str), "lookup key must be a string in "+self._name+" map"
        
        if name in self._xMap:
            assert timestep is not None, "must give timestep for differential state lookup ("+self._name+")"
            if nicpIdx is None:
                nicpIdx = 0
            if degIdx is None:
                degIdx = 0
            assert timestep < (self._nk+1), \
                   "timestep: "+str(timestep)+" out of range in "+self._name+" map (nk: "+str(nk)+")"
            assert degIdx >=0 and degIdx < (self._deg+1), \
                   "degIdx: "+str(deg)+" out of range in "+self._name+" map (deg: "+str(self._deg)+")"
            if timestep is self._nk:
                assert nicpIdx==0 and degIdx==0,"last timestep is only defined at nicpIdx=0,degIdx=0"
            if setVal is None:
                return self._xMap[name][timestep][nicpIdx][degIdx]
            else:
                oldval = self._xMap[name][timestep][nicpIdx][degIdx]
                if (not quiet) and (oldval is not None):
                    print "WARNING: changing \"%s\" %s at timestep %d from %s to %s" % (name,self._name,timestep,str(oldval),str(setVal))
                self._xMap[name][timestep][nicpIdx][degIdx] = setVal

        elif name in self._zMap:
            assert timestep is not None, "must give timestep for algebraic state lookup ("+self._name+")"
            if nicpIdx is None:
                nicpIdx = 0
            assert degIdx is not None, "must set degIdx for algebraic state map ("+self._name+")"
            assert timestep < self._nk, \
                   "timestep: "+str(timestep)+" out of range in "+self._name+" map (nk: "+str(nk)+")"
            assert degIdx >=0 and degIdx < (self._deg+1), \
                   "degIdx: "+str(degIdx)+" out of range in "+self._name+" map (deg: "+str(self._deg)+")"
            if setVal is None:
                return self._zMap[name][timestep][nicpIdx][degIdx]
            else:
                oldval = self._zMap[name][timestep][nicpIdx][degIdx]
                if (quiet is False) and (oldval is not None):
                    print "WARNING: changing \"%s\" %s at timestep %d from %s to %s" % (name,self._name,timestep,str(oldval),str(setVal))
                self._zMap[name][timestep][nicpIdx][degIdx] = setVal

        elif name in self._uMap:
            assert timestep is not None, "must give timestep for control input lookup ("+self._name+")"
            assert nicpIdx is None, "nicpIdx invalid for control input ("+self._name+")"
            assert degIdx is None, "degIdx invalid for control input ("+self._name+")"
            assert timestep < self._nk, \
                   "timestep: "+str(timestep)+" out of range in "+self._name+" map (nk: "+str(nk)+")"
            if setVal is None:
                return self._uMap[name][timestep]
            else:
                oldval = self._uMap[name][timestep]
                if (quiet is False) and (oldval is not None):
                    print "WARNING: changing \"%s\" %s at timestep %d from %s to %s" % (name,self._name,timestep,str(oldval),str(setVal))
                self._uMap[name][timestep] = setVal

        elif name in self._pMap:
            assert timestep is None, "timestep invalid for parameter lookup ("+self._name+")"
            assert nicpIdx is None, "nicpIdx invalid for parameter lookup ("+self._name+")"
            assert degIdx is None, "degIdx invalid for parameter lookup ("+self._name+")"
            if setVal is None:
                return self._pMap[name]
            else:
                oldval = self._pMap[name]
                if (force is False) and (oldval is not None):
                    raise ValueError("can't change \""+name+"\" "+self._name+" once it's set unless you use force=True (tried to change "+str(oldval)+" to "+str(setVal))
                self._pMap[name] = setVal

        else:
            raise KeyError("couldn't find \""+name+"\" in "+self._name+" map")
