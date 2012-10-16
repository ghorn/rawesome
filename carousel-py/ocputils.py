import numpy
import numbers

from dae import Dae
from dvmap import DesignVarMap

import casadi as C

def setFXOptions(fun, options):
    assert(isinstance(options,list))
    for intOpt in options:
        assert(isinstance(intOpt,tuple))
        assert(len(intOpt)==2)
        assert(isinstance(intOpt[0],str))
        optName,optVal = intOpt
        fun.setOption(optName, optVal)

class MultipleShootingStage():
    def __init__(self, dae, nSteps):
        # check inputs
        assert(isinstance(dae, Dae))
        assert(isinstance(nSteps, int))

        # make sure dae has everything
        assert(hasattr(dae,'odeRes'))
        assert(hasattr(dae,'algRes'))
        assert(hasattr(dae,'stateDotDummy'))
        
        self.dae = dae
        self.dae._freeze('MultipleShootingStage(dae)')

        # set up design vars
        self.nSteps = nSteps
        self.states = C.msym("x" ,self.nStates(),self.nSteps)
        self.actions = C.msym("u",self.nActions(),self.nSteps)
        self.params = C.msym("p",self.nParams())
        
#        self._dvs = C.msym("dv",self.nStates()*self.nSteps+self.nActions()*self.nSteps+self.nParams())
#
#        numXVars = self.nStates()*self.nSteps
#        numUVars = self.nActions()*self.nSteps
#        numPVars = self.nParams()
#        
#        self.states  = C.reshape(self._dvs[:numXVars], [self.nStates(), self.nSteps])
#        self.actions = C.reshape(self._dvs[numXVars:numXVars+numUVars], [self.nActions(), self.nSteps])
#        self.params = self._dvs[numXVars+numUVars:]
#
#        assert( self.states.size() == numXVars )
#        assert( self.actions.size() == numUVars )
#        assert( self.params.size() == numPVars )

        # set up interface
        self._constraints = Constraints()
        self._bounds = Bounds(self.dae._xNames, self.dae._uNames, self.dae._pNames, self.nSteps)
        self._initialGuess = InitialGuess(self.dae._xNames, self.dae._uNames, self.dae._pNames, self.nSteps)
        self._designVars = DesignVars((self.dae._xNames,self.states),
                                      (self.dae._uNames,self.actions),
                                      (self.dae._pNames,self.params),
                                      self.nSteps)

    def nStates(self):
        return self.dae.xVec().size()
    def nActions(self):
        return self.dae.uVec().size()
    def nParams(self):
        return self.dae.pVec().size()

    def setIdasIntegrator(self, integratorOptions=[]):
        # make dae input fun
        daeSXFun = self.dae.sxFun()
        assert(self.nStates()==daeSXFun.inputSX(C.DAE_X).size())
        assert(self.nActions()+self.nParams()==daeSXFun.inputSX(C.DAE_P).size())

        # make integrator
        self.integrator = C.IdasIntegrator(daeSXFun)
        setFXOptions(self.integrator, integratorOptions)
        self.integrator.init()

        # set up dynamics constraints
        for k in range(0,self.nSteps-1):
            uk   = self.actions[:,k]
            ukp1 = self.actions[:,k+1]
            xk   = self.states[:,k]
            xkp1 = self.states[:,k+1]
            p = self.params
            upk = C.veccat([uk,p])
            self.addConstraint(self.integrator.call([xk,upk])[C.INTEGRATOR_XF],'==',xkp1)

    # constraints
    def addConstraint(self,lhs,comparison,rhs):
        if hasattr(self, '_solver'):
            raise ValueError("Can't add a constraint once the solver has been set")
        self._constraints.add(lhs,comparison,rhs)

    # bounds
    def setBound(self,name,val,**kwargs):
        self._bounds.setBound(name,val,**kwargs)

    # initial guess
    def setGuess(self,name,val,**kwargs):
        self._initialGuess.setGuess(name,val,**kwargs)
    def setXGuess(self,*args,**kwargs):
        self._initialGuess.setXVec(*args,**kwargs)
    def setUGuess(self,*args,**kwargs):
        self._initialGuess.setUVec(*args,**kwargs)

    def getDesignVars(self):
        return C.veccat( [C.flatten(self.states), C.flatten(self.actions), C.flatten(self.params)] )
        #return self._dvs

    # design vars
    def lookup(self,name,timestep=None):
        return self._designVars.lookup(name,timestep)
    def devectorize(self,xup):
        return self._designVars.devectorize(xup)
    def getTimestepsFromDvs(self,dvs):
        return self._designVars.getTimestepsFromDvs(dvs)

    # solver
    def setObjective(self, objective):
        if hasattr(self, '_objective'):
            raise ValueError("You've already set an objective and you can't change it")
        self._objective = objective
        
    def setSolver(self, solver, solverOptions=[], objFunOptions=[], constraintFunOptions=[]):
        if hasattr(self, '_solver'):
            raise ValueError("You've already set a solver and you can't change it")
        if not hasattr(self, '_objective'):
            raise ValueError("You need to set an objective")

        # make objective function
        f = C.MXFunction([self.getDesignVars()], [self._objective])
        setFXOptions(f, objFunOptions)
        f.init()

        # make constraint function
        g = C.MXFunction([self.getDesignVars()], [self._constraints.getG()])
        setFXOptions(g, constraintFunOptions)
        g.init()

        def mkParallelG():
            gs = [C.MXFunction([self.getDesignVars()],[gg]) for gg in self._constraints._g]
            for gg in gs:
                gg.init()
            
            pg = C.Parallelizer(gs)
#            pg.setOption("parallelization","openmp")
            pg.setOption("parallelization","serial")
#            pg.setOption("parallelization","expand")
            pg.init()
    
            dvsDummy = C.msym('dvs',(self.nStates()+self.nActions())*self.nSteps+self.nParams())
            g_ = C.MXFunction([dvsDummy],[C.veccat(pg.call([dvsDummy]*len(gs)))])
            g_.init()
            return g_

#        parallelG = mkParallelG()
#        guess = self._initialGuess.vectorize()
#        parallelG.setInput([x*1.1 for x in guess])
#        g.setInput([x*1.1 for x in guess])
#    
#        g.evaluate()
#        parallelG.evaluate()
    
#        print parallelG.output()-g.output()
#        exit(0)

        # make solver function
        # self._solver = solver(f, mkParallelG())
        self._solver = solver(f, g)
        setFXOptions(self._solver, solverOptions)
        self._solver.init()

        # set constraints
        self._solver.setInput(self._constraints.getLb(), C.NLP_LBG)
        self._solver.setInput(self._constraints.getUb(), C.NLP_UBG)

        self.setBounds()

    def setBounds(self):
        lb,ub = self._bounds.get()
        self._solver.setInput(lb, C.NLP_LBX)
        self._solver.setInput(ub, C.NLP_UBX)

    def solve(self):
        self._solver.setInput(self._initialGuess.vectorize(), C.NLP_X_INIT)
        self._solver.solve()


class Constraints():
    def __init__(self):
        self._g = []
        self._glb = []
        self._gub = []
        
    def add(self,lhs,comparison,rhs):
        #print "\n\nadding constraint\nlhs: "+str(lhs)+"\ncomparison: "+comparison+"\nrhs: "+str(rhs)
        if comparison=="==":
            g = lhs - rhs
            self._g.append(g)
            self._glb.append(numpy.zeros(g.size()))
            self._gub.append(numpy.zeros(g.size()))
        elif comparison=="<=":
            g = lhs - rhs
            self._g.append(g)
            self._glb.append(-numpy.inf*numpy.ones(g.size()))
            self._gub.append(numpy.zeros(g.size()))
        elif comparison==">=":
            g = rhs - lhs
            self._g.append(g)
            self._glb.append(-numpy.inf*numpy.ones(g.size()))
            self._gub.append(numpy.zeros(g.size()))
        else:
            raise ValueError('Did not recognize comparison \"'+str(comparison)+'\"')

    def getG(self):
        return C.veccat(self._g)
    def getLb(self):
        return C.veccat(self._glb)
    def getUb(self):
        return C.veccat(self._gub)


class Bounds(DesignVarMap):
    descriptor = "bound"
        
    def setBound(self,name,val,**kwargs):
        assert(isinstance(name,str))
        assert(isinstance(val,tuple))
        assert(len(val)==2)
        assert(isinstance(val[0],numbers.Real))
        assert(isinstance(val[1],numbers.Real))
        self.dvmapSet(name,val,**kwargs)

    def get(self):
        return zip(*DesignVarMap.vectorize(self))

class InitialGuess(DesignVarMap):
    descriptor = "initial guess"

    def setGuess(self,name,val,**kwargs):
        assert(isinstance(name,str))
        assert(isinstance(val,numbers.Real))
        self.dvmapSet(name,val,**kwargs)

class DesignVars(DesignVarMap):
    descriptor = "design variables"
    def __init__(self, (xNames,states), (uNames,actions), (pNames,params), nSteps):
        DesignVarMap.__init__(self,xNames,uNames,pNames,nSteps)
        for k in range(0,self.nSteps):
            self.setXVec(states[:,k],timestep=k)
            self.setUVec(actions[:,k],timestep=k)
        for k,pname in enumerate(pNames):
            self.dvmapSet(pname,params[k])
