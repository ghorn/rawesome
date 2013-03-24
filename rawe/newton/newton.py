import casadi as C
import numpy as np

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
        residual = self.dae.getResidual()

        xSize = self.dae.xVec().size()
        zSize = self.dae.zVec().size()
        
        if (residual.size() != zSize+xSize):
            print "WARNING: residual.size() != zSize+xSize)"
    
        # residual function
        u = self.dae.uVec()
        xd = self.dae.xVec()
        xa = self.dae.zVec()
        xddot = C.veccat([self.dae.ddt(name) for name in self.dae.xNames()])
        p  = self.dae.pVec()
        
        ffcn = C.SXFunction([xddot,xd,xa,u,p],[residual])
        ffcn.init()

        return ffcn
