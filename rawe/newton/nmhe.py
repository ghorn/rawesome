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

    def _setupDynamicsConstraints(self,endTime):
        # Todo: add parallelization
        # Todo: get endTime right
        g = []
        nicp = 1
        deg = 4
        p = self._dvMap.pVec()
        for k in range(self.nk):
            newton = Newton(LagrangePoly,self.dae,1,nicp,deg,'RADAU')
            newton.setupStuff(endTime)
            
            X0_i = self._dvMap.xVec(k)
            U_i  = self._U[k,:].T
            newton.isolver.setOutput(1,0)
            _, Xf_i = newton.isolver.call([X0_i,U_i,p])
            X0_i_plus = self._dvMap.xVec(k+1)
            g.append(Xf_i-X0_i_plus)
        return g
            
    def makeSolver(self,endTime):
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
        g   = self._constraints.getG()
        glb = self._constraints.getLb()
        gub = self._constraints.getUb()

        gDyn = self._setupDynamicsConstraints(endTime)
        gDynLb = gDynUb = [C.DMatrix.zeros(gg.shape) for gg in gDyn]
        
        g = C.veccat([g]+gDyn)
        glb = C.veccat([glb]+gDynLb)
        gub = C.veccat([gub]+gDynUb)

        self.glb = glb
        self.gub = gub

        # design vars
        V = self._dvMap.vectorize()

        # gradient of arbitraryObj
        if hasattr(self,'_obj'):
            arbitraryObj = self._obj
        else:
            arbitraryObj = 0
        gradF = C.gradient(arbitraryObj,V)
        
        # hessian of lagrangian:
        Js = [C.jacobian(gnf,V) for gnf in self._gaussNewtonObjF]
        gradFgns = [C.mul(J.T,F) for (F,J) in zip(self._gaussNewtonObjF, Js)]
        gaussNewtonHess = sum([C.mul(J.T,J) for J in Js])
        hessL = gaussNewtonHess + C.jacobian(gradF,V)
        
        gradF += sum(gradFgns)
        
        # equality/inequality constraint jacobian
        gfcn = C.MXFunction([V,self._U],[g])
        gfcn.init()
        jacobG = gfcn.jacobian(0,0)
        jacobG.init()

        # function which generates everything needed
        f = sum(self._gaussNewtonObjF)
        if hasattr(self,'_obj'):
            f += self._obj
        
        self.masterFun = C.MXFunction([V,self._U],[hessL, gradF, g, jacobG.call([V,self._U])[0], f])
        self.masterFun.init()

        self.qp = C.CplexSolver(hessL.sparsity(),jacobG.output(0).sparsity())
        self.qp.init()

    def runSolver(self,U):
        # make sure all bounds are set
        (xMissing,pMissing) = self._guessMap.getMissing()
        msg = []
        for name in xMissing:
            msg.append("you forgot to set a guess for \""+name+"\" at timesteps: "+str(xMissing[name]))
        for name in pMissing:
            msg.append("you forgot to set a guess for \""+name+"\"")
        if len(msg)>0:
            raise ValueError('\n'.join(msg))


        lbx,ubx = zip(*(self._boundMap.vectorize()))
        xk = C.DMatrix(list(self._guessMap.vectorize()))

        import time
#        from pylab import *
        for k in range(1000):
            ############# plot stuff ###############
            print "iteration: ",k
#            import nmheMaps
#            xOpt = np.array(xk).squeeze()
#            traj = nmheMaps.VectorizedReadOnlyNmheMap(self.dae,self.nk,xOpt)
#            
#            xs =  np.array([traj.lookup('x',timestep=kk) for kk in range(self.nk+1)] )
#            zs =  np.array([traj.lookup('z',timestep=kk) for kk in range(self.nk+1)] )
#            dxs = np.array([traj.lookup('dx',timestep=kk) for kk in range(self.nk+1)])
#            dzs = np.array([traj.lookup('dz',timestep=kk) for kk in range(self.nk+1)])
#            m = traj.lookup('m')
#            
#            outputMap = nmheMaps.NmheOutputMap(self._outputMapGenerator, xOpt, U)
#            c = np.array([outputMap.lookup('c',timestep=kk) for kk in range(self.nk)])
#            cdot = np.array([outputMap.lookup('cdot',timestep=kk) for kk in range(self.nk)])
#
#            figure()
#            title(str(float(k)))
#            subplot(2,2,1)
#            plot(xs,-zs)
#            ylabel('pos '+str(k))
#            axis('equal')
#            
#            subplot(2,2,2)
#            plot(dxs,-dzs)
#            ylabel('vel')
#            axis('equal')
#
#            subplot(2,2,3)
#            plot(c)
#            ylabel('c')
#
#            subplot(2,2,4)
#            plot(cdot)
#            ylabel('cdot')
#            ##########################################


            self.masterFun.setInput(xk,0)
            self.masterFun.setInput(U,1)
            t0 = time.time()
            self.masterFun.evaluate()
            t1 = time.time()
            print "masterFun delta time: %.3f ms" % ((t1-t0)*1000)
            hessL  = self.masterFun.output(0)
            gradF  = self.masterFun.output(1)
            g      = self.masterFun.output(2)
            jacobG = self.masterFun.output(3)
            f      = self.masterFun.output(4)

            self.qp.setInput(0,      C.QP_X_INIT)
            self.qp.setInput(hessL,  C.QP_H)
            self.qp.setInput(jacobG, C.QP_A)
            self.qp.setInput(gradF,  C.QP_G)

            self.qp.setInput(lbx-xk,C.QP_LBX)
            self.qp.setInput(ubx-xk,C.QP_UBX)
            
            self.qp.setInput(self.glb-g, C.QP_LBA)
            self.qp.setInput(self.gub-g, C.QP_UBA)

            t0 = time.time()
            self.qp.evaluate()
            t1 = time.time()
            print "qp delta time: %.3f ms" % ((t1-t0)*1000)
            print ""
            deltaX = self.qp.output(C.QP_PRIMAL)
#            import scipy.io
#            scipy.io.savemat('hessL.mat',{'hessL':np.array(hessL),
#                                          'gradF':np.array(gradF),
#                                          'x0':0*np.array(deltaX),
#                                          'xopt':np.array(deltaX),
#                                          'lbx':np.array(lbx-xk),
#                                          'ubx':np.array(ubx-xk),
#                                          'jacobG':np.array(jacobG),
#                                          'lba':np.array(self.glb-g),
#                                          'uba':np.array(self.gub-g)})
#            import sys; sys.exit()

#            print deltaX
            xk += deltaX
#        show()
