import numpy as np

import casadi as C

import nmheMaps
from ocputils import Constraints

from newton import Newton
from collocation import LagrangePoly

class Nmhe(object):
    def __init__(self,dae,nk):
        self.dae = dae
        self.nk = nk

        self._gaussNewtonObjF = []

        mapSize = len(self.dae.xNames())*(self.nk+1) + len(self.dae.pNames())
        V = C.msym('dvs',mapSize)
        self._dvMap = nmheMaps.VectorizedReadOnlyNmheMap(self.dae,self.nk,V)

        self._boundMap = nmheMaps.WriteableNmheMap(self.dae,self.nk)
        self._guessMap = nmheMaps.WriteableNmheMap(self.dae,self.nk)

        self._U = C.msym('u',self.nk,len(self.dae.uNames()))
        self._outputMapGenerator = nmheMaps.NmheOutputMapGenerator(self,self._U)
        self._outputMap = nmheMaps.NmheOutputMap(self._outputMapGenerator, self._dvMap.vectorize(), self._U)

        self._constraints = Constraints()
        
    def __call__(self,*args,**kwargs):
        return self.lookup(*args,**kwargs)

    def lookup(self,name,timestep=None):
        try:
            return self._dvMap.lookup(name,timestep=timestep)
        except NameError:
            pass
        try:
            return self._outputMap.lookup(name,timestep)
        except NameError:
            pass
        raise NameError("unrecognized name \""+name+"\"")
        
    def bound(self,name,(lb,ub),timestep=None):
        self._boundMap.setVal(name,(lb,ub),timestep=timestep)

    def guess(self,name,val,timestep=None):
        self._guessMap.setVal(name,val,timestep=timestep)

    def constrain(self,lhs,comparison,rhs,tag='unnamed_constraint'):
        self._constraints.add(lhs,comparison,rhs,tag)
        
    def setObj(self,obj):
        if hasattr(self,'_obj'):
            raise ValueError("don't change the objective function")
        self._obj = obj

    def addGaussNewtonObjF(self,gnF):
        self._gaussNewtonObjF.append(gnF)

    def _setupDynamicsConstraints(self):
        # Todo: add parallelization
        # Todo: add initialization 
        g = []
        nicp = 10
        deg = 4
        p = self._dvMap.pVec()
        for k in range(self.nk):
            newton = Newton(LagrangePoly,self.dae,1,nicp,deg,'RADAU')
            endTime = 0.05
            newton.setupStuff(endTime)
            
            X0_i = self._dvMap.xVec(k)
            U_i  = self._U[k,:].T
            _, Xf_i = newton.isolver.call([X0_i,U_i,p])
            X0_i_plus = self._dvMap.xVec(k+1)
            g.append(Xf_i-X0_i_plus)
        return g
            
    def makeSolver(self):
        # make sure all bounds are set
        (xMissing,pMissing) = self._boundMap.getMissing()
        msg = []
        for name in xMissing:
            msg.append("you forgot to set a bound on \""+name+"\" at timesteps: "+str(xMissing[name]))
        for name in pMissing:
            msg.append("you forgot to set a bound on \""+name+"\"")
        if len(msg)>0:
            raise ValueError('\n'.join(msg))

        # constraints:
        constraints = self._constraints._g
        constraintLbgs = self._constraints._glb
        constraintUbgs = self._constraints._gub

        g = [self._setupDynamicsConstraints()]
        g = []
        h = []
        hlbs = []
        hubs = []
        for k in range(len(constraints)):
            lb = constraintLbgs[k]
            ub = constraintUbgs[k]
            if all(lb==ub):
                g.append(constraints[k]-lb) # constrain to be zero
            else:
                h.append(constraints[k])
                hlbs.append(lb)
                hubs.append(ub)
        g = C.veccat(g)

        h = C.veccat(h)
        hlbs = C.veccat(hlbs)
        hubs = C.veccat(hubs)

        # design vars
        V = self._dvMap.vectorize()

        # gradient of arbitraryObj
        if hasattr(self,'_obj'):
            arbitraryObj = self._obj
        else:
            arbitraryObj = 0
        gradF = C.gradient(arbitraryObj,V)
        
        # hessian of lagrangian:
        J = 0
        for gnf in self._gaussNewtonObjF:
            J += C.jacobian(gnf,V)
        hessL = C.mul(J.T,J) + C.jacobian(gradF,V)
        
        # equality constraint jacobian
        jacobG = C.jacobian(g,V)

        # inequality constraint jacobian
        jacobH = C.jacobian(h,V)

        # function which generates everything needed
        masterFun = C.MXFunction([V,self._U],[hessL, gradF, g, jacobG, h, jacobH])
        masterFun.init()

        ##########  need to setup/solve the following qp:  #########
        #     min   1/2*x.T*hessL*x + gradF.T*x
        #      x
        #     
        #     S.T           g + jacobG*x == 0
        #           hlbs <= h + jacobH*x <= hubs
        ############################################################
