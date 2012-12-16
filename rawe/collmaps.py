import numpy as np
from collutils import mkCollocationPoints
import casadi as C

class ReadOnlyCollMap(object):
    """
    A map of x/z/u/p handling number of timesteps, nicp, and deg.
    Useful function s are "lookup", "{x,z,u,p}Vec", and "vectorize"
    """
    def __init__(self,ocp,name):
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

    def vectorize(self):
        """
        Return all the variables in one vector
        """
        if all([hasattr(self,attr) for attr in ['_xVec','_zVec','uVec','pVec']]):
            V = [self.pVec()]
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
        
    def _lookupOrSet(self,name,timestep,nicpIdx,degIdx,setVal=None,quiet=False,force=False):
        assert isinstance(name,str), "lookup key must be a string in "+self._name+" map"
        
        if name in self._xMap:
            assert timestep is not None, "must give timestep for differential state lookup ("+self._name+")"
            if nicpIdx is None:
                nicpIdx = 0
            if degIdx is None:
                degIdx = 0
            assert timestep <= self._nk, \
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
                    print "WARNING: changing \"%s\" %s at timestep %d from %s to %s" % \
                        (name,self._name,timestep,str(oldval),str(setVal))
                self._xMap[name][timestep][nicpIdx][degIdx] = setVal

        elif name in self._zMap:
            assert timestep is not None, "must give timestep for algebraic state lookup ("+self._name+")"
            if nicpIdx is None:
                nicpIdx = 0
            assert degIdx is not None, "must set degIdx for algebraic state map ("+self._name+")"
            assert degIdx != 0, "algebraic variable ("+self._name+") not defined at degIdx 0"
            assert timestep < self._nk, \
                "timestep: "+str(timestep)+" out of range in "+self._name+" map (nk: "+str(nk)+")"
            assert degIdx > 0 and degIdx <= self._deg, \
                "degIdx: "+str(degIdx)+" out of range in "+self._name+" map (deg: "+str(self._deg)+")"
            if setVal is None:
                return self._zMap[name][timestep][nicpIdx][degIdx]
            else:
                oldval = self._zMap[name][timestep][nicpIdx][degIdx]
                if (quiet is False) and (oldval is not None):
                    print "WARNING: changing \"%s\" %s at timestep %d from %s to %s" % \
                        (name,self._name,timestep,str(oldval),str(setVal))
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
                    print "WARNING: changing \"%s\" %s at timestep %d from %s to %s" % \
                        (name,self._name,timestep,str(oldval),str(setVal))
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
                    msg = "can't change \""+name+"\" "+self._name+" once it's set unless " + \
                        "you use force=True (tried to change "+str(oldval)+" to "+str(setVal)
                    raise ValueError(msg)
                self._pMap[name] = setVal

        else:
            raise NameError("couldn't find \""+name+"\" in "+self._name+" map")


class WriteableCollMap(ReadOnlyCollMap):
    """
    A ReadOnlyCollMap that also has the "setVal" and "fillInMissing" methods
    """
    def setVal(self,name,val,timestep=None,nicpIdx=None,degIdx=None,quiet=False,force=False):
        self._lookupOrSet(name,timestep,nicpIdx,degIdx,setVal=val,quiet=quiet,force=force)
        
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


class VectorizedReadOnlyCollMap(ReadOnlyCollMap):
    """
    A ReadOnlyCollMap meant to play more nicely with the MX class.
    Everything is stored in one vector so you can call
    {x,z,u,p}Vec and get efficient slices, and "vectorize"
    returns the original vector instead of concatenating all the individual elements.
    This is 
    """
    def __init__(self,ocp,name,vec):
        ReadOnlyCollMap.__init__(self,ocp,name)
        self._devectorize(vec)
        
    def vectorize(self):
        return self._vec
            
    def xVec(self,timestep,nicpIdx=None,degIdx=None):
        if nicpIdx is None:
            nicpIdx = 0
        if degIdx is None:
            degIdx = 0
        return self._xVec[timestep][nicpIdx][degIdx]
    def zVec(self,timestep,nicpIdx=None,degIdx=None):
        if nicpIdx is None:
            nicpIdx = 0
        assert (degIdx is not None), "must set degIdx in zVec"
        assert (degIdx != 0), "algebraic variables not defined at tau=0"
        return self._zVec[timestep][nicpIdx][degIdx]
    def uVec(self,timestep):
        return self._uVec[timestep]
    def pVec(self):
        return self._pVec

    def _devectorize(self,V):
        """
        Take a vector and populate internal _{x,z,u,p}Vec, then
        look through and call setVal so that lookup() works as normal
        """
        self._vec = V
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
        
        def setVal(name,val,timestep=None,nicpIdx=None,degIdx=None):
            self._lookupOrSet(name,timestep,nicpIdx,degIdx,setVal=val)

        offset = 0
        
        # Get the parameters
        P = V[offset:offset+NP]
        for k,name in enumerate(self._pNames):
#            self._indexMap.setVal(name,k)
#            self._dvMap.setVal(name,V[offset+k])
            setVal(name,V[offset+k])
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
                        setVal(name,V[offset+m],timestep=k,nicpIdx=i,degIdx=j)
                        
                    if j !=0:
                        XA[k][i][j] = V[offset+ndiff:offset+ndiff+nalg]
                        for m,name in enumerate(self._zNames):
#                            self._indexMap.setVal(name,offset+ndiff+m,timestep=k,nicpIdx=i,degIdx=j)
#                            self._dvMap.setVal(name,V[offset+ndiff+m],timestep=k,nicpIdx=i,degIdx=j)
                            setVal(name,V[offset+ndiff+m],timestep=k,nicpIdx=i,degIdx=j)

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
                setVal(name,V[offset+m],timestep=k)
            offset += nu
        
        # State at end time
        XD[self._nk][0][0] = V[offset:offset+ndiff]
        for m,name in enumerate(self._xNames):
#            self._indexMap.setVal(name,offset+m,timestep=self.nk,degIdx=0,nicpIdx=0)
#            self._dvMap.setVal(name,V[offset+m],timestep=self.nk,degIdx=0,nicpIdx=0)
            setVal(name,V[offset+m],timestep=self._nk,degIdx=0,nicpIdx=0)
        
        offset += ndiff
        assert offset==NV

        self._xVec = XD
        self._zVec = XA
        self._uVec = U
        self._pVec = P


class OutputMapGenerator(object):
    """
    Something which will efficiently generate a map of all outputs.
    The outputs are all computed all at once to ensure no (additional) CSEs are generated.

    On initialization, the function which creates all the outputs from a dv vector is created.
    Then you use it to initialize an OutputMap object
    """
    def __init__(self,ocp,xDot):
        (fAll,(f0,outputNames0)) = ocp.dae.outputsFun()
        self._outputNames0 = outputNames0
        self._outputNames = ocp.dae.outputNames()

        assert (len(self._outputNames0) == f0.getNumOutputs())
        assert (len(self._outputNames) == fAll.getNumOutputs())

        self._nk = ocp.nk
        self._nicp = ocp.nicp
        self._deg = ocp.deg

        outs = []
        for timestepIdx in range(self._nk):
            for nicpIdx in range(self._nicp):
                # outputs defined at tau_i0
                if f0 is not None:
                    outs += f0.call([ocp._dvMap.xVec(timestepIdx,nicpIdx=nicpIdx,degIdx=0),
                                     ocp._dvMap.uVec(timestepIdx),
                                     ocp._dvMap.pVec()])
                # all outputs
                for degIdx in range(1,self._deg+1):
                    if fAll is not None:
                        outs += fAll.call([xDot[timestepIdx,nicpIdx,degIdx],
                                           ocp._dvMap.xVec(timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx),
                                           ocp._dvMap.zVec(timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx),
                                           ocp._dvMap.uVec(timestepIdx),
                                           ocp._dvMap.pVec()])
        # make the function
        self.fEveryOutput = C.MXFunction([ocp._dvMap.vectorize()],outs)
        self.fEveryOutput.init()

class OutputMap(object):
    """
    Initialize this with an outputMapGenerator and a vector of design vars.
    If you pass a symbolic vector you get symbolic outputs with MXFunction.call().
    If you pass a numeric vector you get numeric outputs with MXFunction.setInput(); MXFunction.evaluate(); ..
    """
    def __init__(self,outputMapGenerator,dvs):
        if type(dvs) in [C.MX,C.SXMatrix]:
            allOutputs = outputMapGenerator.fEveryOutput.call([dvs])
        elif type(dvs) in [np.ndarray,C.DMatrix]:
            outputMapGenerator.fEveryOutput.setInput(dvs,0)
            outputMapGenerator.fEveryOutput.evaluate()
            allOutputs = [np.array(outputMapGenerator.fEveryOutput.output(k))
                          for k in range(outputMapGenerator.fEveryOutput.getNumOutputs())]
        else:
            raise TypeError("OutputMap got unrecognized design vector type: "+str(type(dvs)))

        self._outputNames0 = outputMapGenerator._outputNames0
        self._outputNames = outputMapGenerator._outputNames

        self._numOutputs0 = len(self._outputNames0)
        self._numOutputs  = len(self._outputNames)

        self._nk = outputMapGenerator._nk
        self._nicp = outputMapGenerator._nicp
        self._deg = outputMapGenerator._deg

        self._outputs = {}
        self._outputs0 = {}

        for name in self._outputNames0:
            self._outputs0[name] = np.resize(np.array([None]),(self._nk,self._nicp))
        for name in self._outputNames:
            self._outputs[name] = np.resize(np.array([None]),(self._nk,self._nicp,self._deg+1))

        outs = []
        k = 0
        for timestepIdx in range(self._nk):
            for nicpIdx in range(self._nicp):
                # outputs defined at tau_i0
                outs = allOutputs[k:k+self._numOutputs0]
                k += self._numOutputs0
                for name,val in zip(self._outputNames0,outs):
                    self._outputs0[name][timestepIdx,nicpIdx] = val
                
                # all outputs
                for degIdx in range(1,self._deg+1):
                    outs = allOutputs[k:k+self._numOutputs]
                    k += self._numOutputs
                    for name,val in zip(self._outputNames,outs):
                        self._outputs[name][timestepIdx,nicpIdx,degIdx] = val


    def lookup(self,name,timestep=None,nicpIdx=None,degIdx=None):
        if not (name in self._outputs0 or name in self._outputs):
            raise NameError("couldn't find \""+name+"\"")
            
        #assert (timestep is not None), "please specify timestep at which you want to look up output \""+name+"\""
        if nicpIdx is None:
            nicpIdx = 0
        if degIdx is None:
            degIdx = 0
        
        # if degIdx == 0, return value with no algebraic inputs
        if degIdx == 0:
            assert (name in self._outputs0), "output \""+name+"\" is not an explicit function of x/u/p and is not defined at the beginning of the collocation interval, specify degIdx > 0"
            if timestep is None:
                return [self._outputs0[name][k][nicpIdx] for k in range(self.nk)]
            else:
                return self._outputs0[name][timestep][nicpIdx]
                
        # if degIdx != 0, return value which may or may not have algebraic inputs
        if timestep is None:
            return [self._outputs[name][k][nicpIdx][degIdx] for k in range(self.nk)]
        else:
            return self._outputs[name][timestep][nicpIdx][degIdx]


class QuadratureManager(object):
    """
    This is a thing which manages the quadrature states.
    Call setQuadratureDdt with the name of a quadrature state and the name of it's derivative.
    The derivative name must be something that ocp.lookup() will find.
    """
    def __init__(self,ocp):
        self._quadratures = {}
        self._nk = ocp.nk
        self._nicp = ocp.nicp
        self._deg = ocp.deg

    def setQuadratureDdt(self,quadratureStateName,quadratureStateDotName,lookup,lagrangePoly,h,symbolicDvs):
        if quadratureStateName in self._quadratures:
            raise ValueError(quadratureStateName+" is not unique")

        qStates = np.resize(np.array([None],dtype=C.MX),(self._nk+1,self._nicp,self._deg+1))
        # set the dimension of the initial quadrature state correctly, and make sure the derivative name is valid
        try:
            qStates[0,0,0] = 0*lookup(quadratureStateDotName,timestep=0,nicpIdx=0,degIdx=1)
        except ValueError:
            raise ValueError('quadrature derivative name \"'+quadratureStateDotName+
                             '\" is not a valid x/z/u/p/output')

        ldInv = np.linalg.inv(lagrangePoly.lDotAtTauRoot[1:,1:])
        ld0 = lagrangePoly.lDotAtTauRoot[1:,0]
        l1 = lagrangePoly.lAtOne
#        print -C.mul(C.mul(ldInv, ld0).T, l1[1:]) + l1[0]

        breakQuadratureIntervals = True

        if breakQuadratureIntervals:
            for k in range(self._nk):
                for i in range(self._nicp):
                    dQs = h*C.veccat([lookup(quadratureStateDotName,timestep=k,nicpIdx=i,degIdx=j)
                                      for j in range(1,self._deg+1)])
                    Qs = C.mul( ldInv, dQs)
                    for j in range(1,self._deg+1):
                        qStates[k,i,j] = qStates[k,i,0] + Qs[j-1]
                    
                    m = C.mul( Qs.T, l1[1:])
                    if i < self._nicp - 1:
                        qStates[k,i+1,0] = qStates[k,i,0] + m
                    else:
                        qStates[k+1,0,0] = qStates[k,i,0] + m
        else:
            for k in range(self._nk):
                for i in range(self._nicp):
                    dQs = h*C.veccat([lookup(quadratureStateDotName,timestep=k,nicpIdx=i,degIdx=j)
                                      for j in range(1,self._deg+1)])
                    Qs = C.mul( ldInv, dQs - ld0*qStates[k,i,0] )
                    for j in range(1,self._deg+1):
                        qStates[k,i,j] = Qs[j-1]
                    
                    m = C.veccat( [qStates[k,i,0], Qs] )
                    m = C.mul( m.T, l1)
                    if i < self._nicp - 1:
                        qStates[k,i+1,0] = m
                    else:
                        qStates[k+1,0,0] = m
                
        self._quadratures[quadratureStateName] = qStates
        self._setupQuadratureFunctions(symbolicDvs)

    def lookup(self,name,timestep,nicpIdx,degIdx):
        if name not in self._quadratures:
            raise NameError('unrecognized name "'+name+'"')
        assert (timestep is not None), "please specify timestep at which you want to look up output \""+name+"\""
        if nicpIdx is None:
            nicpIdx = 0
        if degIdx is None:
            degIdx = 0
        return self._quadratures[name][timestep][nicpIdx][degIdx]

    def _setupQuadratureFunctions(self,symbolicDvs):
        quadouts = []
        for name in self._quadratures:
            quadouts.extend(list(self._quadratures[name].flatten()[:-self._deg]))
        self.quadratureFun = C.MXFunction([symbolicDvs],quadouts)
        self.quadratureFun.init()
        
class QuadratureMap(object):
    def __init__(self,quadratureManager,numericDvs):
        self._nk = quadratureManager._nk
        self._nicp = quadratureManager._nicp
        self._deg = quadratureManager._deg

        f = quadratureManager.quadratureFun
        f.setInput(numericDvs,0)
        f.evaluate()
        allOutputs = [np.array(f.output(k)) for k in range(f.getNumOutputs())]

        k = 0
        self._quadMap = {}
        numPerState = self._nk*self._nicp*(self._deg+1) + 1

        for name in quadratureManager._quadratures:
            self._quadMap[name] = np.resize(allOutputs[k:k+numPerState],(self._nk+1,self._nicp,self._deg+1))
#            self.quadMap[name][-1,1:,:] = None
#            self.quadMap[name][-1,0,1:] = None
            k += numPerState

    def lookup(self,name,timestep,nicpIdx,degIdx):
        if name not in self._quadMap:
            raise NameError('unrecognized name "'+name+'"')
        assert (timestep is not None), "please specify timestep at which you want to look up output \""+name+"\""
        if nicpIdx is None:
            nicpIdx = 0
        if degIdx is None:
            degIdx = 0
        if (timestep == -1) or (timestep == self._nk):
            assert (degIdx==0 and nicpIdx==0), "quadrature state undefined at last timestep, degIdx>0 and/or nicpIdx>0"
        return self._quadMap[name][timestep][nicpIdx][degIdx]
