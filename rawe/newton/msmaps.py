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

class OutputMapGenerator(object):
    """
    Something which will efficiently generate a map of all outputs.
    The outputs are all computed all at once to ensure no (additional) CSEs are generated.

    On initialization, the function which creates all the outputs from a dv vector is created.
    Then you use it to initialize an OutputMap object
    """
    def __init__(self,ocp):
        (fAll,(f0,outputNames0)) = ocp.dae.outputsFun()
        self._outputNames0 = outputNames0
        self._outputNames = ocp.dae.outputNames()

        assert (len(self._outputNames0) == f0.getNumOutputs())
        assert (len(self._outputNames) == fAll.getNumOutputs())

        self._nk = ocp.nk
        self._nicp = ocp.nicp

        outs = []
        for timestepIdx in range(self._nk):
            for nicpIdx in range(self._nicp):
                # outputs defined at tau_i0
                if f0 is not None:
                    outs += f0.call([ocp._dvMap.xVec(timestepIdx,nicpIdx=nicpIdx),
                                     ocp._dvMap.uVec(timestepIdx),
                                     ocp._dvMap.pVec()])
        # make the function
        self.fEveryOutput = C.SXFunction([ocp._dvMap.vectorize()],outs)
        self.fEveryOutput.init()

class OutputMap(object):
    """
    Initialize this with an outputMapGenerator and a vector of design vars.
    If you pass a symbolic vector you get symbolic outputs with MXFunction.call().
    If you pass a numeric vector you get numeric outputs with MXFunction.setInput(); MXFunction.evaluate(); ..
    """
    def __init__(self,outputMapGenerator,dvs):
        if type(dvs) == C.MX:
            allOutputs = outputMapGenerator.fEveryOutput.call([dvs])
        if type(dvs) == C.SXMatrix:
            allOutputs = outputMapGenerator.fEveryOutput.eval([dvs])
        elif type(dvs) in [np.ndarray,C.DMatrix]:
            outputMapGenerator.fEveryOutput.setInput(dvs,0)
            outputMapGenerator.fEveryOutput.evaluate()
            allOutputs = [np.array(outputMapGenerator.fEveryOutput.output(k)).squeeze()
                          for k in range(outputMapGenerator.fEveryOutput.getNumOutputs())]
        else:
            raise TypeError("OutputMap got unrecognized design vector type: "+str(type(dvs)))

        self._outputNames0 = outputMapGenerator._outputNames0
        self._outputNames = outputMapGenerator._outputNames

        self._numOutputs0 = len(self._outputNames0)
        self._numOutputs  = len(self._outputNames)

        self._nk = outputMapGenerator._nk
        self._nicp = outputMapGenerator._nicp

        self._outputs0 = {}

        for name in self._outputNames0:
            self._outputs0[name] = np.resize(np.array([None]),(self._nk,self._nicp))

        outs = []
        k = 0
        for timestepIdx in range(self._nk):
            for nicpIdx in range(self._nicp):
                # outputs defined at tau_i0
                outs = allOutputs[k:k+self._numOutputs0]
                k += self._numOutputs0
                for name,val in zip(self._outputNames0,outs):
                    self._outputs0[name][timestepIdx,nicpIdx] = val

    def lookup(self,name,timestep=None,nicpIdx=None):
        if not (name in self._outputs0 or name in self._outputs):
            raise NameError("couldn't find \""+name+"\"")

        if not (name in self._outputs):
            raise NameError("sorry, \""+name+"\" depends on algebraic variable or ddt(differential variable) \
                            and Multiple Shooting cannot access it")

        #assert (timestep is not None), "please specify timestep at which you want to look up output \""+name+"\""
        if nicpIdx is None:
            nicpIdx = 0

        assert (timestep != None), "please set timestep"
        return self._outputs0[name][timestep][nicpIdx]


    def makesolver(self):
        V = self._dvMap.vectorize()
        lb,ub = unzip(self._boundMap.vectorize()
        guess = self._guessMap.vectorize()

        ffcn = SXFunction([V],[self.obj])
        gfcn = SXFunction([V],[])
        IpoptSolver(

