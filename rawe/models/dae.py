
import casadi as C

class Dae(object):
    """
    Class to hold represent a differential-algebraic or ordinary differential equation
    """
    def __init__(self):
        # set of reasons (strings) why the Dae cannot be modified (add new x/z/u/p/output)
        self._frozen = set()

        # lists of names (strings)
        self._xNames = []
        self._zNames = []
        self._uNames = []
        self._pNames = []
        self._outputNames = []

        # dictionaries of SXMatrix symbols
        self._syms = {}

        # map of derivatives
        self._dummyDdtMap = {}
        
    def _freeze(self,msg):
        self._frozen.add(msg)

    def assertNotFrozen(self):
        if len(self._frozen) > 0:
            raise ValueError("can't perform this operation because Dae has been frozen by: "+str([n for n in self._frozen]))

    def assertUniqueName(self, name):
        allNames = self._xNames + self._zNames + self._uNames + self._pNames + self._outputNames
        if name in allNames:
            raise ValueError('name "'+name+'" is not unique')

    def _addVar(self,name,namelist):
        self.assertNotFrozen()
        if isinstance(name,list):
            return [self._addVar(n,namelist) for n in name]

        assert(isinstance(name,str))
        
        self.assertUniqueName(name)
        namelist.append(name)

        ret = C.ssym(name)
        self._syms[name] = ret
        return ret

    def ddt(self,name):
        if name not in self._xNames:
            raise ValueError("unrecognized state \""+name+"\"")
        try:
            return self._dummyDdtMap[name]
        except KeyError:
            self._dummyDdtMap[name] = C.ssym(name+"DotDummy_____")
            return self.ddt(name)

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
        """
        Add a differential state
        """
        return self._addVar(name,self._xNames)

    def addZ(self,name):
        """
        Add an algebraic variable
        """
        return self._addVar(name,self._zNames)

    def addU(self,name):
        """
        Add a control input
        """
        return self._addVar(name,self._uNames)

    def addP(self,name):
        """
        Add a parameter
        """
        return self._addVar(name,self._pNames)
    
    def xVec(self):
        self._freeze('xVec()')
        return C.veccat([self._syms[n] for n in self._xNames])
    def zVec(self):
        self._freeze('zVec()')
        return C.veccat([self._syms[n] for n in self._zNames])
    def uVec(self):
        self._freeze('uVec()')
        return C.veccat([self._syms[n] for n in self._uNames])
    def pVec(self):
        self._freeze('pVec()')
        return C.veccat([self._syms[n] for n in self._pNames])

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

    def __getitem__(self,name):
        """
        Get a differential state, algebraic var, control, param, or output
        """
        if not isinstance(name,str):
            raise KeyError('Dae key must be a string')

        try:
            return self._syms[name]
        except KeyError:
            raise KeyError(name+' is not a symbol in Dae')
        
    def __setitem__(self,name,val):
        """
        Add an output
        """
        self.assertNotFrozen()
        if not isinstance(name,str):
            raise KeyError('Output name must be a string')

#        assert( isinstance(val, C.SXMatrix) )

        self.assertUniqueName(name)
        self._outputNames.append(name)
        self._syms[name] = val

    def outputsFun(self):
        self._freeze('outputsFun()')

        outputsNoZ = [] # only outputs with no algebraic variables

        # which outputs have no algebraic vars
        for name in self.outputNames():
            # try to make a function without any algebraic variables
            f = C.SXFunction([self.xVec(),self.uVec(),self.pVec()],
                             [self[name]]
                             )
            f.init()
            if len(f.getFree()) == 0:
                # only add if there are no algebraic variables
                outputsNoZ.append(name)
            
        # function with only outputs with no algebraic vars
        if len(outputsNoZ)>0:
            fNoZ = C.SXFunction([self.xVec(), self.uVec(), self.pVec()],
                                [self[name] for name in outputsNoZ])
            fNoZ.setOption('name','outputs with no algebraic vars')
            fNoZ.init()
        else:
            fNoZ = None

        # function which outputs everything
        if len(self.outputNames())>0:
            fAll = C.SXFunction([self.xVec(), self.zVec(), self.uVec(), self.pVec()],
                                [self[name] for name in self.outputNames()])
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
            
        xdot = C.veccat([self.ddt(name) for name in self.xNames()])
        return C.SXFunction( C.daeIn( x=self.xVec(),
                                      z=self.zVec(),
                                      p=C.veccat([self.uVec(),self.pVec()]),
                                      xdot=xdot ),
                             C.daeOut( alg=algRes, ode=odeRes) )

if __name__=='__main__':
    dae = Dae()

    [pos,vel,mass] = dae.addX( ["pos","vel","mass"] )
    [zdummy] = dae.addZ( ["zdummy"] )
    thrust = dae.addU( "thrust" )
    
    dae['pos*vel'] = pos*vel

    dae.setOdeRhs([vel, (thrust - 0.05*vel*vel)/mass, -0.1*thrust*thrust])
    dae.setAlgRes([zdummy])
