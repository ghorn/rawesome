
import casadi as C
import acadoSimExport
import acadoModelExport
from octaveSimExport import generateOctaveSim

class Dae(object):
    """
    Class to hold represent a differential-algebraic or ordinary differential equation
    """
    def __init__(self):
        # set of reasons (strings) why the Dae cannot be modified (add new x/z/u/p/output)
        self._xzupFrozen = set()
        self._outputsFrozen = set()

        # lists of names (strings)
        self._xNames = []
        self._zNames = []
        self._uNames = []
        self._pNames = []
        self._outputNames = []

        # list of illegal names
        self._illegalNames = ['']

        # dictionaries of SXMatrix symbols
        self._syms = {}

        # map of derivatives
        self._dummyDdtMap = {}
        
    def _freezeXzup(self,msg):
        self._xzupFrozen.add(msg)

    def _freezeOutputs(self,msg):
        self._outputsFrozen.add(msg)

    def assertXzupNotFrozen(self):
        if len(self._xzupFrozen) > 0:
            raise ValueError("can't perform this operation because Dae has been frozen by: "+str([n for n in self._xzupFrozen]))

    def assertOutputsNotFrozen(self):
        if len(self._outputsFrozen) > 0:
            raise ValueError("can't perform this operation because Dae has been frozen by: "+str([n for n in self._outputsFrozen]))

    def assertUniqueName(self, name):
        allNames = self._xNames + self._zNames + self._uNames + self._pNames + self._outputNames + self._illegalNames
        if name in allNames:
            raise ValueError('name "'+name+'" is not unique or illegal')

        if name[0]=='_':
            raise ValueError("underscores before names are illegal")

    def _addVar(self,name,namelist):
        self.assertXzupNotFrozen()
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
            self._dummyDdtMap[name] = C.ssym('_DotDummy_'+name)
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
        self._freezeXzup('xVec()')
        return C.veccat([self._syms[n] for n in self._xNames])
    def zVec(self):
        self._freezeXzup('zVec()')
        return C.veccat([self._syms[n] for n in self._zNames])
    def uVec(self):
        self._freezeXzup('uVec()')
        return C.veccat([self._syms[n] for n in self._uNames])
    def pVec(self):
        self._freezeXzup('pVec()')
        return C.veccat([self._syms[n] for n in self._pNames])
    def xDotVec(self):
        self._freezeXzup('xDotVec()')
        return C.veccat([self.ddt(n) for n in self._xNames])

    def xNames(self):
        self._freezeXzup('xNames()')
        return self._xNames
    def zNames(self):
        self._freezeXzup('zNames()')
        return self._zNames
    def uNames(self):
        self._freezeXzup('uNames()')
        return self._uNames
    def pNames(self):
        self._freezeXzup('pNames()')
        return self._pNames
    def outputNames(self):
        self._freezeOutputs('outputNames()')
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
        self.assertOutputsNotFrozen()
        if not isinstance(name,str):
            raise KeyError('Output name must be a string')

        self.assertUniqueName(name)
        self._outputNames.append(name)
        self._syms[name] = val

    def outputsFun(self):
        self._freezeOutputs('outputsFun()')

        # which outputs are defined at tau_i0
        outputs0 = [] # only outputs which are defined at tau_i0
        for name in self.outputNames():
            # try to make a function without any algebraic or ddt(x) variables
            f = C.SXFunction([self.xVec(),self.uVec(),self.pVec()],
                             [self[name]]
                             )
            f.init()
            if len(f.getFree()) == 0:
                # only add if there are no algebraic or ddt(x) variables
                outputs0.append(name)
            
        # function with outputs defined at tau_i0
        if len(outputs0)>0:
            f0 = C.SXFunction([self.xVec(), self.uVec(), self.pVec()],
                              [self[name] for name in outputs0])
            f0.setOption('name','outputs defined at tau_i0')
            f0.init()
        else:
            f0 = None

        # function with all outputs (defined at collocation points)
        if len(self.outputNames())>0:
            xdot = C.veccat([self.ddt(name) for name in self.xNames()])
            fAll = C.SXFunction([xdot, self.xVec(), self.zVec(), self.uVec(), self.pVec()],
                                [self[name] for name in self.outputNames()])
            fAll.init()
            fAll.setOption('name','all outputs')
        else:
            fAll = None

        return (fAll,(f0,outputs0))

    def getResidual(self):
        self._freezeXzup('getResidual()')
        if not hasattr(self,'_odeRes'):
            raise ValueError('need to set the residual')
        f = self._odeRes
        if isinstance(f,list):
            f = C.veccat(f)

        if hasattr(self,'_algRes'):
            algRes = self._algRes
            if isinstance(algRes,list):
                algRes = C.veccat(algRes)
            f = C.veccat([f,algRes])
        return f

    def casadiDae(self):
#       I have:
#       0 = fg(xdot,x,z)
#
#       I need dot(x') = f(x',z')
#                    0 = g(x',z')
#
#       let x' = x
#           z' = [z,xdot]
#           dot(x') = xdot
#           0 = g(x',z') = g(x,[z,xdot]) = fg(xdot,x,z)
        self._freezeXzup('casadiDae()')
        f = self.getResidual()
            
        xdot = C.veccat([self.ddt(name) for name in self.xNames()])
        return C.SXFunction( C.daeIn( x=self.xVec(),
                                      z=C.veccat([self.zVec(),xdot]),
                                      p=C.veccat([self.uVec(),self.pVec()])
                                      ),
                             C.daeOut( alg=f, ode=xdot) )

    def octaveSimGen(self,functionName):
        self._freezeXzup('octaveSimGen()')
        return generateOctaveSim(self,functionName)

    def acadoSimGen(self):
        self._freezeXzup('agadoSimGen()')

        f = self._odeRes
        if isinstance(f,list):
            f = C.veccat(f)
        if hasattr(self,'_algRes'):
            algRes = self._algRes
            if isinstance(algRes,list):
                algRes = C.veccat(algRes)
            f = C.veccat([f,algRes])
            
        xdot = C.veccat([self.ddt(name) for name in self.xNames()])

        info = { 'x':self.xVec(),
                 'z':self.zVec(),
                 'p':self.pVec(),
                 'u':self.uVec(),
                 'xdot':xdot,
                 'f':f }
        return acadoSimExport.simExport(self, info)

    def acadoModelGen(self):
        self._freezeXzup('agadoModelGen()')
        self._freezeOutputs('agadoModelGen()')

        f = self._odeRes
        if isinstance(f,list):
            f = C.veccat(f)
        if hasattr(self,'_algRes'):
            algRes = self._algRes
            if isinstance(algRes,list):
                algRes = C.veccat(algRes)
            f = C.veccat([f,algRes])
            
        return acadoModelExport.generateAcadoCodegenModel(self,f)

    def convertToOde(self):
        # get the residual fg(xdot,x,z)
        fg = self.getResidual()

        # take the jacobian w.r.t. xdot,z
        jac = C.jacobian(fg,C.veccat([self.xDotVec(), self.zVec()]))

        # make sure that it was linear in {xdot,z}, i.e. the jacobian is not a function of {xdot,z}
        testJac = C.SXFunction([self.xVec(),self.uVec(),self.pVec()], [jac])
        testJac.init()
        assert len(testJac.getFree()) == 0, \
            "can't convert dae to ode, residual jacobian is a function of {xdot,z}"

        # get the constant term
        fg_fun = C.SXFunction([self.xVec(),self.zVec(),self.uVec(),self.pVec(),self.xDotVec()], [fg])
        fg_fun.init()
        [fg_zero] = fg_fun.eval([self.xVec(),0*self.zVec(),self.uVec(),self.pVec(),0*self.xDotVec()])
        testFun = C.SXFunction([self.xVec(),self.uVec(),self.pVec()], [fg_zero])
        testFun.init()
        assert len(testFun.getFree()) == 0, \
            "the \"impossible\" happened in Dae -> Ode conversion"

        xDotAndZ = C.solve(jac, -fg_zero)
        xDot = xDotAndZ[0:len(self.xNames())]
        z = xDotAndZ[len(self.xNames()):]

        # construct outputs and residual as a function of x/z/u/p/xdot
        oldFunctions = C.SXFunction([self.xVec(), self.zVec(), self.uVec(), self.pVec(), self.xDotVec()],
                                    [xDot,z]+[self[name] for name in self.outputNames()])
        oldFunctions.init()

        # construct outputs and residual as a function of x/u/p only
        newOutputs = oldFunctions.eval([self.xVec(), z, self.uVec(), self.pVec(), xDot])
        newFunctions = C.SXFunction([self.xVec(),self.uVec(),self.pVec()], newOutputs)
        newFunctions.init()

        # use this to construct a new Dae which has no algebraic states
        dae = Dae()
        xs = C.veccat([dae.addX(name) for name in self.xNames()])
        us = C.veccat([dae.addU(name) for name in self.uNames()])
        ps = C.veccat([dae.addP(name) for name in self.pNames()])
        outs = newFunctions.eval([xs,us,ps])
        newXdot = outs[0]
        newZ = outs[1]
        odeRes = []
        for k,name in enumerate(self.xNames()):
            odeRes.append( dae.ddt(name) - newXdot[k] )
        dae.setOdeRes(C.veccat(odeRes))
        for k,name in enumerate(self.zNames()):
            dae[name] = newZ[k]
        for k,name in enumerate(self.outputNames()):
            dae[name] = outs[k+2]
        return dae

if __name__=='__main__':
    dae = Dae()

    [pos,vel,mass] = dae.addX( ["pos","vel","mass"] )
    [zdummy] = dae.addZ( ["zdummy"] )
    thrust = dae.addU( "thrust" )
    
    dae['pos*vel'] = pos*vel

    dae.setOdeRhs([vel, (thrust - 0.05*vel*vel)/mass, -0.1*thrust*thrust])
    dae.setAlgRes([zdummy])
