import casadi as C
import numpy as np

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

class Nmpc(object):
    def __init__(self,dae,nk):
        self.dae = dae
        self.nk = nk

        X = C.ssym('x',dae.xVec().size(),self.nk+1)
        U = C.ssym('u',dae.zVec().size(),self.nk)
        p = C.ssym('p',dae.pVec().size())
        self._dvMap = NmpcMap(self.dae,self.nk,X,U,p)

        X = np.resize(np.array([None]),(dae.xVec().size(),self.nk+1))
        U = np.resize(np.array([None]),(dae.zVec().size(),self.nk))
        p = np.resize(np.array([None]),dae.pVec().size())
        self._boundMap = NmpcMap(self.dae,self.nk,X,U,p)

        X = np.resize(np.array([None]),(dae.xVec().size(),self.nk+1))
        U = np.resize(np.array([None]),(dae.zVec().size(),self.nk))
        p = np.resize(np.array([None]),dae.pVec().size())
        self._guessMap = NmpcMap(self.dae,self.nk,X,U,p)

    def __call__(self,*args,**kwargs):
        return self.lookup(*args,**kwargs)

    def lookup(self,name,timestep=None):
        return self._dvMap.lookup(name,timestep=timestep)

    def bound(self,name,(lb,ub),timestep=None):
        self._boundMap.setVal(name,(lb,ub),timestep=timestep)
        
#    def setObj(self,obj):
#        if hasattr(self,'_obj'):
#            raise ValueError("don't change the objective function")
#        self._obj = obj

    def setGaussNewtonObjF(self,gnF):
        if hasattr(self,'_gnF'):
            raise ValueError("don't change the objective function")
        self._gnF = gnF
        
class Newton(object):
    def __init__(self,LagrangePoly,dae,nk,nicp,deg,collPoly):
        self.dae = dae
        self.nk = nk
        self.nicp = nicp
        self.deg = deg
        self.collPoly = collPoly
        self.lagrangePoly = LagrangePoly(deg=self.deg,collPoly=self.collPoly)

    def setupStuff(self,endTime):
        self.h = endTime/float(self.nk*self.nicp)
        assert isinstance(self.h,float), "gauss newton doesn't support free end time yet"

        ifcn = self._makeImplicitFunction()
        self.isolver = self._makeImplicitSolver(ifcn)

    def _makeImplicitSolver(self,ifcn):
#        self.implicitSolver = C.KinsolSolver(ifcn)
        #implicitSolver = C.NLPImplicitSolver(ifcn)
        #implicitSolver.setOption("nlp_solver",C.IpoptSolver)
        implicitSolver = C.NewtonImplicitSolver(ifcn)
        implicitSolver.setOption("linear_solver",C.CSparse)
        implicitSolver.init()
        return implicitSolver

    def _makeImplicitFunction(self):
        ffcn = self._makeResidualFunction()
        x0 = C.ssym('x0',self.dae.xVec().size())
        X =  C.ssym('X',self.dae.xVec().size(),self.deg*self.nicp)
        Z =  C.ssym('Z',self.dae.zVec().size(),self.deg*self.nicp)
        u =  C.ssym('u',self.dae.uVec().size())
        p =  C.ssym('p',self.dae.pVec().size())
        constraints = []
        ############################################################
        
        ndiff = self.dae.xVec().size()

        x0_ = x0
        for nicpIdx in range(self.nicp):
            X_ = X[:,nicpIdx*self.deg:(nicpIdx+1)*self.deg]
            Z_ = Z[:,nicpIdx*self.deg:(nicpIdx+1)*self.deg]
            for j in range(1,self.deg+1):
                # Get an expression for the state derivative at the collocation point
                xp_jk = self.lagrangePoly.lDotAtTauRoot[j,0]*x0_
                for j2 in range (1,self.deg+1):
                    # get the time derivative of the differential states (eq 10.19b)
                    xp_jk += self.lagrangePoly.lDotAtTauRoot[j,j2]*X_[:,j2-1]
                # Add collocation equations to the NLP
                [fk] = ffcn.eval([xp_jk/self.h,
                                  X_[:,j-1],
                                  Z_[:,j-1],
                                  u,
                                  p])
                
                # impose system dynamics (for the differential states (eq 10.19b))
                constraints.append(fk[:ndiff])
                
                # impose system dynamics (for the algebraic states (eq 10.19b))
                constraints.append(fk[ndiff:])
                
            # Get an expression for the state at the end of the finite element
            xf = self.lagrangePoly.lAtOne[0]*x0_
            for j in range(1,self.deg+1):
                xf += self.lagrangePoly.lAtOne[j]*X_[:,j-1]
            x0_ = xf
     
        ifcn = C.SXFunction([C.veccat([X,Z]),x0,u,p],[C.veccat(constraints),xf])
        ifcn.init()
        return ifcn


    def _makeResidualFunction(self):
        if not hasattr(self.dae, '_odeRes'):
            raise ValueError("need to set ode residual")

        residual = self.dae._odeRes

        xSize = self.dae.xVec().size()
        zSize = self.dae.zVec().size()
        
        if hasattr(self.dae,'_algRes'):
            residual = C.veccat([residual, self.dae._algRes])

        assert (residual.size() == zSize+xSize)
    
        # residual function
        u = self.dae.uVec()
        xd = self.dae.xVec()
        xa = self.dae.zVec()
        xddot = C.veccat([self.dae.ddt(name) for name in self.dae.xNames()])
        p  = self.dae.pVec()
        
        ffcn = C.SXFunction([xddot,xd,xa,u,p],[residual])
        ffcn.init()

        return ffcn
