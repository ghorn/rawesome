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

import casadi as CS
import numpy as np
import numbers
import sys

from rawe.ocputils import Constraints,setFXOptions
import collmaps
from collpoints import mkCollocationPoints
from rawe.dae import Dae
import trajectory

class LagrangePoly(object):
    def __init__(self,deg=None,collPoly=None):
        assert deg is not None
        assert collPoly is not None

        self.deg = deg
        self.collPoly = collPoly
        self.tau_root = mkCollocationPoints(self.collPoly,self.deg)

        self._mkLagrangePolynomials()

    def _mkLagrangePolynomials(self):
        # Collocation point
        tau = CS.ssym("_tau")

        # lagrange polynomials
        self.lfcns = []

        # lagrange polynomials evaluated at collocation points
        self.lAtOne = np.zeros(self.deg+1)

        # derivative of lagrange polynomials evaluated at collocation points
        self.lDotAtTauRoot = np.zeros((self.deg+1,self.deg+1))

        # For all collocation points: eq 10.4 or 10.17 in Biegler's book
        # Construct Lagrange polynomials to get the polynomial basis at the collocation point
        for j in range(self.deg+1):
            L = 1
            for k in range(self.deg+1):
                if k != j:
                    L *= (tau-self.tau_root[k])/(self.tau_root[j]-self.tau_root[k])
            lfcn = CS.SXFunction([tau],[L])
            lfcn.init()
            self.lfcns.append(lfcn)
            # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
            lfcn.setInput(1.0)
            lfcn.evaluate()
            self.lAtOne[j] = lfcn.output()
            # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the collocation equation
            for k in range(self.deg+1):
                lfcn.setInput(self.tau_root[k])
                lfcn.setFwdSeed(1.0)
                lfcn.evaluate(1,0)
                self.lDotAtTauRoot[k,j] = lfcn.fwdSens()

    def interp(self,tau,zs):
        ret = 0.0
        for j in range(deg+1):
            self.lfcns[j].setInput(tau)
            self.lfcns[j].evaluate()
            ret += self.lfcn.output()*zs[j]

        return ret

class BoundsMap(collmaps.WriteableCollMap):
    def __init__(self,ocp,name):
        collmaps.WriteableCollMap.__init__(self,ocp,name)

        # tags for the bounds
        bndtags = []
        for name in ocp.dae.pNames():
            bndtags.append((name,None))
        for k in range(ocp.nk):
            for i in range(ocp.nicp):
                for j in range(ocp.deg+1):
                    if k==0 and j==0 and i==0:
                        for name in ocp.dae.xNames():
                            bndtags.append((name,(k,i,j)))
                    else:
                        if j!=0:
                            for name in ocp.dae.xNames()+ocp.dae.zNames():
                                bndtags.append((name,(k,i,j)))
                        else:
                            for name in ocp.dae.xNames():
                                bndtags.append((name,(k,i,j)))
            for name in ocp.dae.uNames():
                bndtags.append((name,(k,i,j)))
        for name in ocp.dae.xNames():
            bndtags.append((name,(ocp.nk,0,0)))

        self.bndtags = bndtags

    def boundsFeedback(self,x,lbx,ubx,reportThreshold=0):
        """
        Tests if x >= ub + reportThreshold
                 x <= lb - reportThreshold
        Positive reportThreshold supresses barely active bounds
        Negative reportThreshold reports not-quite-active bounds
        """
        violations = {}

        ubviols = x - ubx
        lbviols = lbx - x
        ubviolsIdx = np.where(CS.logic_and(ubviols >= reportThreshold, ubx > lbx))[0]
        lbviolsIdx = np.where(CS.logic_and(lbviols >= reportThreshold, ubx > lbx))[0]
        violations = {}
        for k in ubviolsIdx:
            (name,time) = self.bndtags[k]
            viol = ('ub',time,float(ubviols[k]))
            if name not in violations:
                violations[name] = [viol]
            else:
                violations[name].append(viol)
        for k in lbviolsIdx:
            (name,time) = self.bndtags[k]
            viol = ('lb',time,float(lbviols[k]))
            if name not in violations:
                violations[name] = [viol]
            else:
                violations[name].append(viol)
        return violations

    def boundsFeedbackStr(self,*args,**kwargs):
        viols = self.boundsFeedback(*args,**kwargs)
        ret = []
        for name in viols:
            ret.append("bound violation! \""+name+": "+str(viols[name]))
        return '\n'.join(ret)

    def printBoundsFeedback(self,*args,**kwargs):
        s = self.boundsFeedbackStr(*args,**kwargs)
        if s != '':
            print s

class Coll():
    collocationIsSetup = False
    def __init__(self, dae, nk=None, nicp=1, deg=4, collPoly='RADAU'):
        assert nk is not None
        assert isinstance(dae, Dae)

        self.dae = dae

        self.nk = nk
        self.nicp = nicp
        self.deg = deg
        self.collPoly = collPoly

        self._bounds = BoundsMap(self,"bounds")
        self._guess = collmaps.WriteableCollMap(self,"guess")

        self._constraints = Constraints()

        # setup NLP variables
        self._dvMap = collmaps.VectorizedReadOnlyCollMap(self,"design var map",CS.msym("V",self.getNV()))

        # quadratures
        self._quadratureManager = collmaps.QuadratureManager(self)


    def setQuadratureDdt(self,quadratureStateName,quadratureStateDotName):
        ''' Add a new quadrature state to the collocation problem by specifying the name of its derivative
        '''
        # run some checks and pass it to the less safe CollMapPlus.setQuadratureDdt
        if not self.collocationIsSetup:
            raise ValueError("Can't add quadratures until you call setupCollocation")

        # make sure this is a unique name (quadratureManager also checks)
        self.dae.assertUniqueName(quadratureStateName)

        # make sure nobody adds an output named the same as quadratureStateName
        self.dae._illegalNames.append(quadratureStateName)

        # setup the quadrature state
        self._quadratureManager.setQuadratureDdt(quadratureStateName,quadratureStateDotName,
                                                 self.lookup,self.lagrangePoly,self.h,self._dvMap.vectorize())

        """
        V = struct_msym([
           entry("x",repeat=[N,nicp,deg+1],struct=struct(dae.xNames())),
           entry("xN",struct=struct(dae.xNames())),
           entry("z",repeat=[N,nicp,deg],struct=struct(dae.zNames())),
           entry("u",repeat=N,struct=struct(dae.uNames())),
           entry("p",struct=struct(dae.pNames()))
        ])

        V.cat

        solution = V(solver.output())

        solution["x",:]

        lbx = V(0)
        ubx = V(0)

        lbx["x",0] = ubx["x",0] = 123
        """

    def setupCollocation(self,tf):
        if self.collocationIsSetup:
            raise ValueError("you can't setup collocation twice")
        self.collocationIsSetup = True

        ## -----------------------------------------------------------------------------
        ## Collocation setup
        ## -----------------------------------------------------------------------------
        # Size of the finite elements
        self.h = tf/float(self.nk*self.nicp)

        # make coefficients for collocation/continuity equations
        self.lagrangePoly = LagrangePoly(deg=self.deg,collPoly=self.collPoly)

        # function to get h out
        self.hfun = CS.MXFunction([self._dvMap.vectorize()],[self.h])
        self.hfun.init()

        # add collocation constraints
        ffcn = self._makeResidualFun()

        ndiff = self.xSize()
        nalg = self.zSize()

        self._xDot = np.resize(np.array([None]),(self.nk,self.nicp,self.deg+1))

        # For all finite elements
        for k in range(self.nk):
            for i in range(self.nicp):
                # For all collocation points
                for j in range(1,self.deg+1):
                    # Get an expression for the state derivative at the collocation point
                    xp_jk = 0
                    for j2 in range (self.deg+1):
                        # get the time derivative of the differential states (eq 10.19b)
                        xp_jk += self.lagrangePoly.lDotAtTauRoot[j,j2]*self.xVec(k,nicpIdx=i,degIdx=j2)
                    self._xDot[k,i,j] = xp_jk/self.h
                    # Add collocation equations to the NLP
                    [fk] = ffcn.call([self._xDot[k,i,j],
                                      self.xVec(k,nicpIdx=i,degIdx=j),
                                      self.zVec(k,nicpIdx=i,degIdx=j),
                                      self.uVec(k),
                                      self.pVec()])

                    # impose system dynamics (for the differential states (eq 10.19b))
                    self.constrain(fk,'==',0,tag=("implicit dynamic equation",(k,i,j)))

                # Get an expression for the state at the end of the finite element
                xf_k = 0
                for j in range(self.deg+1):
                    xf_k += self.lagrangePoly.lAtOne[j]*self.xVec(k,nicpIdx=i,degIdx=j)
#                    print "self.lagrangePoly.lAtOne["+str(j)+"]:" +str(self.lagrangePoly.lAtOne[j])

#                mxfun = CS.MXFunction([self._V],[xf_k])
#                mxfun.init()
#                sxfun = CS.SXFunction(mxfun)
#                sxfun.init()
#                print ""
#                print sxfun.outputSX()

                # Add continuity equation to NLP
                if i==self.nicp-1:
                    self.constrain(self.xVec(k+1,nicpIdx=0,degIdx=0), '==', xf_k, tag=("continuity",(k,i)))
                else:
                    self.constrain(self.xVec(k,nicpIdx=i+1,degIdx=0), '==', xf_k, tag=("continuity",(k,i)))

        # add outputs
        self._outputMapGenerator = collmaps.OutputMapGenerator(self, self._xDot)
        self._outputMap = collmaps.OutputMap(self._outputMapGenerator, self._dvMap.vectorize())

    def xVec(self,*args,**kwargs):
        return self._dvMap.xVec(*args,**kwargs)
    def zVec(self,*args,**kwargs):
        return self._dvMap.zVec(*args,**kwargs)
    def uVec(self,*args,**kwargs):
        return self._dvMap.uVec(*args,**kwargs)
    def pVec(self,*args,**kwargs):
        return self._dvMap.pVec(*args,**kwargs)

    def constrain(self,lhs,comparison,rhs,tag=('unnamed_constraint',None)):
        self._constraints.add(lhs,comparison,rhs,tag=tag)

    def constrainBnds(self,g,(lbg,ubg),tag=('unnamed_constraint',None)):
        self._constraints.addBnds(g,(lbg,ubg),tag=tag)

    def xSize(self):
        return len(self.dae.xNames())
    def zSize(self):
        return len(self.dae.zNames())
    def uSize(self):
        return len(self.dae.uNames())
    def pSize(self):
        return len(self.dae.pNames())

    # Total number of variables
    def getNV(self):
        NXD = self.nicp*self.nk*(self.deg+1)*self.xSize() # Collocated differential states
        NXA = self.nicp*self.nk*self.deg*self.zSize()     # Collocated algebraic states
        NU = self.nk*self.uSize()               # Parametrized controls
        NXF = self.xSize()                 # Final state (only the differential states)
        NP = self.pSize()
        return NXD+NXA+NU+NXF+NP

    def _makeResidualFun(self):
        residual = self.dae.getResidual()

        if (residual.size() != self.zSize()+self.xSize()):
            print 'WARNING: residual.size() != self.zSize() + self.xSize() ==> (%d != %d + %d)' % (residual.size(),self.zSize(), self.xSize())

        # residual function
        u = self.dae.uVec()
        xd = self.dae.xVec()
        xa = self.dae.zVec()
        xddot = CS.veccat([self.dae.ddt(name) for name in self.dae.xNames()])
        p  = self.dae.pVec()

        ffcn = CS.SXFunction([xddot,xd,xa,u,p],[residual])
        ffcn.init()

        return ffcn

    def interpolateInitialGuess(self,filename,force=False,quiet=False,numLoops=1):
        traj = trajectory.load_traj(filename)
        pps = trajectory.make_pps(traj)

        h = (traj.tgrid[-1,0,0] - traj.tgrid[0,0,0])/float(traj.dvMap._nk*traj.dvMap._nicp)
        h *= traj.dvMap._nk*traj.dvMap._nicp/float(self.nk*self.nicp)
        h *= numLoops

        sys.stdout.write('reticulating splines... ')
        missing = []

        t0 = 0.0
        x_interp_times = []
        for timestepIdx in range(self.nk):
            for nicpIdx in range(self.nicp):
                for degIdx in range(self.deg+1):
                    time = t0 + h*self.lagrangePoly.tau_root[degIdx]
                    if time > traj.tgrid[-1,0,0]:
                        time -= traj.tgrid[-1,0,0]
                    x_interp_times.append(time)
                t0 += h
                if t0 > traj.tgrid[-1,0,0]:
                    t0 -= traj.tgrid[-1,0,0]
        x_interp_times.append(t0)

        t0 = 0.0
        z_interp_times = []
        for timestepIdx in range(self.nk):
            for nicpIdx in range(self.nicp):
                for degIdx in range(1,self.deg+1):
                    time = t0 + h*self.lagrangePoly.tau_root[degIdx]
                    if time > traj.tgrid[-1,0,0]:
                        time -= traj.tgrid[-1,0,0]
                    z_interp_times.append(time)
                t0 += h

        # interpolate differential states
        for name in self.dae.xNames():
            sys.stdout.write(name+' '); sys.stdout.flush()
            if name not in pps:
                missing.append(name)
                continue
            # evaluate piecewise poly to set initial guess
            x_interp_values =  pps[name](x_interp_times)
            k = 0
            for timestepIdx in range(self.nk):
                for nicpIdx in range(self.nicp):
                    for degIdx in range(self.deg+1):
                        self.guess(name,x_interp_values[k],
                                   timestep=timestepIdx,nicpIdx=nicpIdx,
                                   degIdx=degIdx,force=force,quiet=quiet)
                        k += 1
            self.guess(name,x_interp_values[k],timestep=-1,nicpIdx=0,degIdx=0,force=force,quiet=quiet)

        # interpolate algebraic variables
        for name in self.dae.zNames():
            sys.stdout.write(name+' '); sys.stdout.flush()
            if name not in pps:
                missing.append(name)
                continue
            # evaluate piecewise poly to set initial guess
            z_interp_values =  pps[name](z_interp_times)
            k = 0
            for timestepIdx in range(self.nk):
                for nicpIdx in range(self.nicp):
                    for degIdx in range(1,self.deg+1):
                        self.guess(name,z_interp_values[k],
                                   timestep=timestepIdx,nicpIdx=nicpIdx,
                                   degIdx=degIdx,force=force,quiet=quiet)
                        k += 1

        # interpolate controls
        for name in self.dae.uNames():
            sys.stdout.write(name+' '); sys.stdout.flush()
            if name not in pps:
                missing.append(name)
                continue
            # evaluate piecewise poly to set initial guess
            u_interp_values = pps[name]([h*k for k in range(self.nk)])
            for timestepIdx in range(self.nk):
                self.guess(name,u_interp_values[timestepIdx],
                           timestep=timestepIdx,force=force,quiet=quiet)

        # set parameters
        for name in self.dae.pNames():
            sys.stdout.write(name+' '); sys.stdout.flush()
            if name not in traj.dvMap._pNames:
                missing.append(name)
                continue
            if name=='endTime':
                self.guess(name,traj.dvMap.lookup(name)*numLoops,force=force,quiet=quiet)
            else:
                self.guess(name,traj.dvMap.lookup(name),force=force,quiet=quiet)
        sys.stdout.write('\n')

        msg = "finished interpolating initial guess"
        if len(missing) > 0:
            msg += ", couldn't find fields: "+str(missing)
        else:
            msg += ", all fields found"
        print msg

    def setupSolver(self,solverOpts=[],constraintFunOpts=[],callback=None):
        if not self.collocationIsSetup:
            raise ValueError("you forgot to call setupCollocation")

        g =   self._constraints.getG()
        lbg = self._constraints.getLb()
        ubg = self._constraints.getUb()

        # Objective function/constraints of the NLP
        if not hasattr(self,'_objective'):
            raise ValueError('need to set objective function')
        nlp = CS.MXFunction(CS.nlpIn(x=self._dvMap.vectorize()),CS.nlpOut(f=self._objective, g=g))
        setFXOptions(nlp,constraintFunOpts)
        nlp.init()

        # solver callback (optional)
        if callback is not None:
            nd = self._dvMap.vectorize().size()
            nc = self._constraints.getG().size()

            c = CS.PyFunction( callback,
                               CS.nlpSolverOut(x = CS.sp_dense(nd,1),
                                               f = CS.sp_dense(1,1),
                                               lam_x = CS.sp_dense(nd,1),
                                               lam_g = CS.sp_dense(nc,1),
                                               lam_p = CS.sp_dense(0,1),
                                               g = CS.sp_dense(nc,1) ),
                               [CS.sp_dense(1,1)] )
            c.init()
            solverOpts.append( ("iteration_callback", c) )

        # Allocate an NLP solver
        self.solver = CS.IpoptSolver(nlp)
#        self.solver = CS.WorhpSolver(nlp)
#        self.solver = CS.SQPMethod(nlp)

        # Set options
        setFXOptions(self.solver, solverOpts)

        # initialize the solver
        self.solver.init()

        # Bounds on g
        self.solver.setInput(lbg,'lbg')
        self.solver.setInput(ubg,'ubg')

        ## Nonlinear constraint function, for debugging
        gfcn = CS.MXFunction([self._dvMap.vectorize()],[g])
        gfcn.init()
        setFXOptions(gfcn,constraintFunOpts)
        self._gfcn = gfcn

    def solve(self,xInit=None,warnZBounds=False,warnZGuess=False):
        if not hasattr(self, 'solver'):
            raise Exception("you need to call setupSolver before you can call solve")
        for name in self.dae.zNames():
            for k in range(self.nk):
                for j in range(self.nicp):
                    for d in range(1,self.deg+1):
                        # if algebraic states are missing from bounds, set to (-inf,inf)
                        if self._bounds.lookup(name,timestep=k,nicpIdx=j,degIdx=d) is None:
                            self._bounds.setVal(name,(-np.inf,np.inf),timestep=k,nicpIdx=j,degIdx=d)
                        # if algebraic states are missing from guess, set to 0
                        if self._guess.lookup(name,timestep=k,nicpIdx=j,degIdx=d) is None:
                            self._guess.setVal(name,0.0,timestep=k,nicpIdx=j,degIdx=d)

        # fill in missing differential states and make sure everything required is set
        def linearInterpMissing(tau,val0,val1):
            return val0*(1-tau) + val1*tau
        self._guess.fillInMissing("initial guess",linearInterpMissing)

        def ceilMissing(tau,val0,val1):
            (lb0,ub0) = val0
            (lb1,ub1) = val1
            lb = min(lb0,lb1)
            ub = max(ub0,ub1)
            return (lb,ub)
        self._bounds.fillInMissing("bounds",ceilMissing)

        vars_init = self._guess.vectorize()
        vars_lbub = self._bounds.vectorize()
        vars_lb = [lbub[0] for lbub in vars_lbub]
        vars_ub = [lbub[1] for lbub in vars_lbub]

        if xInit is not None:
            vars_init = xInit

        # Initial condition
        self.solver.setInput(vars_init, 'x0')

        # Bounds on x
        self.solver.setInput(vars_lb,'lbx')
        self.solver.setInput(vars_ub,'ubx')

        # Solve the problem
        self.solver.solve()
        traj = trajectory.TrajectoryPlotter(self,np.array(self.solver.output('x')))

        # print active bounds/constraints
        def printBoundViolation():
            # bounds feedback
            xOpt = traj.dvMap.vectorize()
            lbx = self.solver.input('lbx')
            ubx = self.solver.input('ubx')
            self._bounds.printBoundsFeedback(xOpt,lbx,ubx,reportThreshold=0)

        def printConstraintViolation():
            xOpt = traj.dvMap.vectorize()
            lbg = np.array(self.solver.input('lbg'))
            ubg = np.array(self.solver.input('ubg'))
            self._gfcn.setInput(xOpt,0)
            self._gfcn.evaluate()
            g = self._gfcn.output()
            self._constraints.printViolations(g,lbg,ubg,reportThreshold=0)
        printBoundViolation()
        printConstraintViolation()

        ret = self.solver.getStat('return_status')
        assert ret in ['Solve_Succeeded','Solved_To_Acceptable_Level'], 'Solver failed: '+ret

        # Print the optimal cost
        print "optimal cost: ", float(self.solver.output('f'))

        # Retrieve the solution
        return traj

    def bound(self,name,val,timestep=None,quiet=False,force=False):
        assert isinstance(name,str)
        assert isinstance(val,tuple)
        assert len(val)==2
        assert isinstance(val[0],numbers.Real)
        assert isinstance(val[1],numbers.Real)

        # handle timestep == None
        if timestep is None:
            if name in self._bounds._xMap:
                for k in range(self.nk+1):
                    self.bound(name,val,timestep=k,quiet=quiet)
                return
            if name in self._bounds._zMap:
                for k in range(self.nk):
                    self.bound(name,val,timestep=k,quiet=quiet)
                return
            if name in self._bounds._uMap:
                for k in range(self.nk):
                    self.bound(name,val,timestep=k,quiet=quiet)
                return
            if name in self._bounds._pMap:
                self._bounds.setVal(name,val,force=force)
                return
            raise ValueError("can't bound \""+name+"\" because it's not x/z/u/p")

        assert isinstance(timestep,int)
        if name in self._bounds._zMap:
            for nicpIdx in range(self.nicp):
                for degIdx in range(1,self.deg+1):
                    self._bounds.setVal(name,val,timestep=timestep,nicpIdx=nicpIdx,degIdx=degIdx,quiet=quiet)
        else:
            self._bounds.setVal(name,val,timestep=timestep,quiet=quiet)

    def guess(self,name,val,timestep=None,nicpIdx=None,degIdx=None,quiet=False,force=False):
        assert isinstance(name,str)
        if (type(val)==np.array or type(val)==np.ndarray) and val.size == 1:
            val = float(val)
        assert isinstance(val,numbers.Real), "your guess for \""+name+"\" must be a valid python number, you gave: "+repr(val)

        # handle timestep == None
        if timestep is None:
            assert (nicpIdx is None)
            assert (degIdx is None)
            if name in self._guess._xMap:
                for k in range(self.nk+1):
                    self.guess(name,val,timestep=k,quiet=quiet)
                return
            elif name in self._guess._zMap:
                for k in range(self.nk):
                    self.guess(name,val,timestep=k,quiet=quiet)
                return
            elif name in self._guess._uMap:
                for k in range(self.nk):
                    self.guess(name,val,timestep=k,quiet=quiet)
                return
            elif name in self._guess._pMap:
                self._guess.setVal(name,val,force=force)
                return
            else:
                raise Exception('"'+name+'" is not a design variable')

        assert isinstance(timestep,int), "timestep is not an int: " + str(timestep)
        if nicpIdx is not None:
            assert isinstance(nicpIdx,int)
        if degIdx is not None:
            assert isinstance(degIdx,int)
        self._guess.setVal(name,val,timestep=timestep,nicpIdx=nicpIdx,degIdx=degIdx,quiet=quiet)

    def _guessVec(self,val,names,**kwargs):
        if isinstance(val,list):
            length = len(val)
        elif hasattr(val,'size'):
            if hasattr(val.size,'__call__'):
                length = val.size()
            else:
                length = val.size
        else:
            raise ValueError("can't figure out how long "+str(val)+" is")
        assert len(names)==length,'guess{X,Z,U,P} got wrong length for guess'

        for k,name in enumerate(names):
            self.guess(name,float(val[k]),**kwargs)

    def guessX(self,val,**kwargs):
        self._guessVec(val,self.dae.xNames(),**kwargs)
    def guessZ(self,val,**kwargs):
        self._guessVec(val,self.dae.zNames(),**kwargs)
    def guessU(self,val,**kwargs):
        self._guessVec(val,self.dae.uNames(),**kwargs)
    def guessP(self,val,**kwargs):
        self._guessVec(val,self.dae.pNames(),**kwargs)

    def __call__(self,*args,**kwargs):
        return self.lookup(*args,**kwargs)

    def lookup(self,name,timestep=None,nicpIdx=None,degIdx=None):
        if name in self.dae.outputNames():
            if not hasattr(self,'_outputMap'):
                raise ValueError("Can't lookup outputs until you call setupOutputs")
            else:
                return self._outputMap.lookup(name,timestep=timestep,nicpIdx=nicpIdx,degIdx=degIdx)
        try:
            return self._dvMap.lookup(name,timestep=timestep,nicpIdx=nicpIdx,degIdx=degIdx)
        except NameError:
            pass
        try:
            return self._quadratureManager.lookup(name,timestep,nicpIdx,degIdx)
        except NameError:
            pass
        raise NameError("lookup fail, unrecognized name "+name)

    def setObjective(self,obj):
        if hasattr(self,'_objective'):
            raise ValueError("can't change objective function once it's already set")
        self._objective = obj
