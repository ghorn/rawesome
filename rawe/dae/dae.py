# Copyright 2012-2013 Greg Horn
#
# This file is part of rawesome.
#
# rawesome is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rawesome is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with rawesome.  If not, see <http://www.gnu.org/licenses/>.

import re
import casadi as C

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
        '''
        call this when it's illegal for user to add any more states/controls/params
        to lock the model
        '''
        self._xzupFrozen.add(msg)

    def _freezeOutputs(self,msg):
        '''
        call this when it's illegal for user to add any more outputs
        to lock the model
        '''
        self._outputsFrozen.add(msg)

    def _getAllNames(self):
        return self._xNames + self._zNames + self._uNames + self._pNames + self._outputNames

    def assertXzupNotFrozen(self):
        if len(self._xzupFrozen) > 0:
            raise ValueError("can't perform this operation because Dae has been frozen by: "+str([n for n in self._xzupFrozen]))

    def assertOutputsNotFrozen(self):
        if len(self._outputsFrozen) > 0:
            raise ValueError("can't perform this operation because Dae has been frozen by: "+str([n for n in self._outputsFrozen]))

    def assertUniqueName(self, name):
        '''
        raise an error if given name is already used, or is illegal
        '''
        if name in (self._getAllNames() + self._illegalNames):
            raise ValueError('name "'+name+'" is not unique or illegal')

        # make sure it's a valid variable name (useful for codegen)
        if not re.match("^[A-Za-z0-9_]+$", name):
            raise Exception('invalid name: "'+name+'"')

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
        '''
        get the time derivative of a differential state
        '''
        if name not in self._xNames:
            raise ValueError("unrecognized state \""+name+"\"")
        try:
            return self._dummyDdtMap[name]
        except KeyError:
            self._dummyDdtMap[name] = C.ssym('_DotDummy_'+name)
            return self.ddt(name)

    def setResidual(self,res):
        '''
        set the dae implicit residual vector f
        where f(x,xdot,z,u,p) == 0
        '''
        if hasattr(self,'_residual'):
            raise ValueError('residual already set')
        if isinstance(res,list):
            res = C.veccat(res)
        self._residual = res

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
        '''
        get differential state vector
        '''
        self._freezeXzup('xVec()')
        return C.veccat([self._syms[n] for n in self._xNames])
    def zVec(self):
        '''
        get algebraic variable
        '''
        self._freezeXzup('zVec()')
        return C.veccat([self._syms[n] for n in self._zNames])
    def uVec(self):
        '''
        get control vector
        '''
        self._freezeXzup('uVec()')
        return C.veccat([self._syms[n] for n in self._uNames])
    def pVec(self):
        '''
        get parameter vector
        '''
        self._freezeXzup('pVec()')
        return C.veccat([self._syms[n] for n in self._pNames])
    def xDotVec(self):
        '''
        get differential state derivative vector
        '''
        self._freezeXzup('xDotVec()')
        return C.veccat([self.ddt(n) for n in self._xNames])

    def xNames(self):
        '''
        get list of names of all differential states
        '''
        self._freezeXzup('xNames()')
        return self._xNames
    def zNames(self):
        '''
        get list of names of all algebraic variables
        '''
        self._freezeXzup('zNames()')
        return self._zNames
    def uNames(self):
        '''
        get list of names of all controls
        '''
        self._freezeXzup('uNames()')
        return self._uNames
    def pNames(self):
        '''
        get list of names of all parameters
        '''
        self._freezeXzup('pNames()')
        return self._pNames
    def outputNames(self):
        '''
        get list of names of all outputs
        '''
        self._freezeOutputs('outputNames()')
        return self._outputNames

    def __contains__(self,name):
        if not isinstance(name,str):
            raise KeyError('Dae key must be a string')
        if name in self._getAllNames():
            return True
        return False

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

        try:
            self._syms[name] = C.DMatrix(val)
        except:
            self._syms[name] = val

    def outputsFun(self):
        '''
        return (fAll, (f0,outputs0)) where
        fAll = all outputs as fAll(xdot, x, z, u, p)
        f0 = outputs defined by f0(x, u, p)
        outputs0 is list of names of f0 outputs
        '''
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

    def outputsFunWithSolve(self):
        '''
        this solves for unknown xdot and z symbolically
        then gives all outputs as function of only f(x,u,p)
        '''
        # get output fun as fcn of [xdot, x, z, u, p]
        (fAll, _) = self.outputsFun()
        if fAll == None:
            f = C.SXFunction([self.xVec(), self.uVec(), self.pVec()], [0])
            return f
        # solve for xdot, z
        (xDotDict, zDict) = self.solveForXDotAndZ()
        xDot = C.veccat([xDotDict[name] for name in self.xNames()])
        z    = C.veccat([zDict[name] for name in self.zNames()])
        # plug in xdot, z solution to outputs fun
        fAll.init()
        outputs = fAll.eval([xDot, self.xVec(), z, self.uVec(), self.pVec()])
        # make new SXFunction that is only fcn of [x, u, p]
        f = C.SXFunction([self.xVec(), self.uVec(), self.pVec()], outputs)

        f.init()
        assert len(f.getFree()) == 0, 'the "impossible" happened >_<'
        return f

    def getResidual(self):
        '''
        return the dae residual vector
        '''
        self._freezeXzup('getResidual()')
        if not hasattr(self,'_residual'):
            raise ValueError('need to set the residual')
        return self._residual

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

    def solveForXDotAndZ(self):
        '''
        returns (xDotDict,zDict) where these dictionaries contain symbolic
        xdot and z which are only a function of x,u,p        
        '''
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
            "the \"impossible\" happened in solveForXDotAndZ"

        xDotAndZ = C.solve(jac, -fg_zero)
        xDot = xDotAndZ[0:len(self.xNames())]
        z = xDotAndZ[len(self.xNames()):]

        xDotDict = {}
        for k,name in enumerate(self.xNames()):
            xDotDict[name] = xDot[k]
        zDict = {}
        for k,name in enumerate(self.zNames()):
            zDict[name] = z[k]
        return (xDotDict, zDict)

    def convertToOde(self):
        '''
        EXPERIMENTAL!
        Solve for xdot and z as fcns of x,u,p, then backsubstitude these in.
        return dae with no xdot or z, and residual/outputs are same as before.
        '''
        (xDotDict,zDict) = self.solveForXDotAndZ()
        xDot = C.veccat([xDotDict[name] for name in self.xNames()])
        z = C.veccat([zDict[name] for name in self.zNames()])

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
        res = []
        for k,name in enumerate(self.xNames()):
            res.append( dae.ddt(name) - newXdot[k] )
        dae.setResidual(res)
        for k,name in enumerate(self.zNames()):
            dae[name] = newZ[k]
        for k,name in enumerate(self.outputNames()):
            dae[name] = outs[k+2]
        return dae

    def assertNoFreeParams(self):
        '''
        throw exception if there are some symbolics in the dae which were not added
        with add{X,Z,U,P} or dae.ddt(blah)
        '''
        # get the residual fg(xdot,x,z,u,p)
        fg = self.getResidual()
        alloutputs = [fg] + [self[name] for name in self.outputNames()]

        testFun = C.SXFunction([self.xDotVec(),self.xVec(),self.zVec(),self.uVec(),self.pVec()],
                               alloutputs)
        testFun.init()
        assert len(testFun.getFree()) == 0, "oh noes, dae has free parameters: "+str(testFun.getFree())
