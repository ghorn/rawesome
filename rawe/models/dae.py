
import casadi as C

class Dae():
    def __init__(self):
        self._frozen = set()
        
        self._xNames = []
        self._zNames = []
        self._uNames = []
        self._pNames = []
        self._outputNames = []

        self._xDict = {}
        self._zDict = {}
        self._uDict = {}
        self._pDict = {}
        self._outputDict = {}

    def _freeze(self,msg):
        self._frozen.add(msg)

    def assertNotFrozen(self):
        if len(self._frozen) > 0:
            raise ValueError("can't perform this operation because Dae has been frozen by: "+str([n for n in self._frozen]))

    def assertUniqueName(self, name):
        allNames = self._xNames + self._zNames + self._uNames + self._pNames + self._outputNames
        if name in allNames:
            raise ValueError('name "'+name+'" is not unique')

    def _getVar(self,name,namelist,namedict):
        if isinstance(name,list):
            return map(lambda n: self._getVar(n,namelist,namedict), name)

        if isinstance(name,tuple):
            return tuple(self._getVar(list(name),namelist,namedict))
                
        assert(isinstance(name,str))

        if name not in namedict:
            raise KeyError(name+' is not in '+str(namelist))
        
        return namedict[name]

    def _addVar(self,name,namelist,namedict):
        self.assertNotFrozen()
        if isinstance(name,list):
            return [self._addVar(n,namelist,namedict) for n in name]

        assert(isinstance(name,str))
        
        if name in namedict:
            raise KeyError(name+' is already in in '+str(namelist))
        
        self.assertUniqueName(name)
        namelist.append(name)
        namedict[name] = C.ssym(name)
        return namedict[name]

    def setAlgRes(self,res):
        if hasattr(self,'_algRes'):
            raise ValueError('algebraic residual already set')
        if isinstance(res,list):
            res = C.veccat(res)
        self._algRes = res

    def setOdeRes(self,res):
        if hasattr(self,'_odeRes'):
            raise ValueError('ode residual already set')
        if isinstance(res,list):
            res = C.veccat(res)
        self._odeRes = res

    def addX(self,name):
        return self._addVar(name,self._xNames,self._xDict)

    def addZ(self,name):
        return self._addVar(name,self._zNames,self._zDict)

    def addU(self,name):
        return self._addVar(name,self._uNames,self._uDict)

    def addP(self,name):
        return self._addVar(name,self._pNames,self._pDict)
    
    def x(self,name):
        return self._getVar(name,self._xNames,self._xDict)

    def z(self,name):
        return self._getVar(name,self._zNames,self._zDict)

    def u(self,name):
        return self._getVar(name,self._uNames,self._uDict)

    def p(self,name):
        return self._getVar(name,self._pNames,self._pDict)

    def output(self,name):
        return self._getVar(name,self._outputNames,self._outputDict)

    def xVec(self):
        self._freeze('xVec()')
        return C.veccat([self._xDict[n] for n in self._xNames])
    def zVec(self):
        self._freeze('zVec()')
        return C.veccat([self._zDict[n] for n in self._zNames])
    def uVec(self):
        self._freeze('uVec()')
        return C.veccat([self._uDict[n] for n in self._uNames])
    def pVec(self):
        self._freeze('pVec()')
        return C.veccat([self._pDict[n] for n in self._pNames])

    def xNames(self):
        self._freeze('xNames()')
        return self._xNames
    def zNames(self):
        self._freeze('zNames()')
        return self._zNames
    def uNames(self):
        self._freeze('uNames()')
        return self._uNames
    def pNames(self):
        self._freeze('pNames()')
        return self._pNames
    def outputNames(self):
        self._freeze('outputNames()')
        return self._outputNames

    def addOutput(self,name,val):
        self.assertNotFrozen()
        assert( isinstance(name, str) )
        assert( isinstance(val, C.SXMatrix) )
        self.assertUniqueName(name)
        
        self._outputNames.append(name)
        self._outputDict[name] = val

    def outputsFun(self):
        self._freeze('outputsFun()')

        outputsNoZ = [] # only outputs with no algebraic variables

        # which outputs have no algebraic vars
        for name in self.outputNames():
            # try to make a function without any algebraic variables
            f = C.SXFunction([self.xVec(),self.uVec(),self.pVec()],
                             [self.output(name)]
                             )
            f.init()
            if len(f.getFree()) == 0:
                # only add if there are no algebraic variables
                outputsNoZ.append(name)
            
        # function with only outputs with no algebraic vars
        if len(outputsNoZ)>0:
            fNoZ = C.SXFunction([self.xVec(), self.uVec(), self.pVec()],
                                [self.output(n) for n in outputsNoZ])
            fNoZ.setOption('name','outputs with no algebraic vars')
            fNoZ.init()
        else:
            fNoZ = None

        # function which outputs everything
        if len(self.outputNames())>0:
            fAll = C.SXFunction([self.xVec(), self.zVec(), self.uVec(), self.pVec()],
                                [self.output(n) for n in self.outputNames()])
            fAll.init()
            fAll.setOption('name','all outputs')
        else:
            fAll = None

        return (fAll,(fNoZ,outputsNoZ))

    def sxFun(self):
        self._freeze('sxFun()')
        algRes = None
        if hasattr(self,'_algRes'):
            algRes = self._algRes
        odeRes = self._odeRes

        if isinstance(odeRes,list):
            odeRes = C.veccat(odeRes)
        if isinstance(algRes,list):
            algRes = C.veccat(algRes)
            
        return C.SXFunction( C.daeIn( x=self.xVec(), z=self.zVec(), p=C.veccat([self.uVec(),self.pVec()]), xdot=self.stateDotDummy ),
                             C.daeOut( alg=algRes, ode=odeRes) )

if __name__=='__main__':
    dae = Dae()
    print dae.xNames
    print dae.xDict
    dae.x('x')
    print dae.xNames
    print dae.xDict
    dae.x('y')
    print dae.xNames
    print dae.xDict
