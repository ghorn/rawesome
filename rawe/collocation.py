import casadi as CS
import numpy as np
import numbers
import pickle
from scipy.interpolate import PiecewisePolynomial

from ocputils import Constraints,setFXOptions
from collmap import CollMap
from collutils import mkCollocationPoints
from models import Dae
from trajectoryData import TrajectoryData

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
        tau = CS.ssym("tau")
          
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
                self.lDotAtTauRoot[j,k] = lfcn.fwdSens()

    def interp(self,tau,zs):
        ret = 0.0
        for j in range(deg+1):
            self.lfcns[j].setInput(tau)
            self.lfcns[j].evaluate()
            ret += self.lfcn.output()*zs[j]
            
        return ret

def boundsFeedback(x,lbx,ubx,bndtags,tolerance=0):
    violations = {}
    
    length = None
    if hasattr(x,'size'):
        if isinstance(x,list):
            length = len(x)
        elif hasattr(x,'size'):
            if hasattr(x.size,'__call__'):
                length = x.size()
            else:
                length = x.size
    if length is None:
        raise ValueError("couldn't determine length of x")

    for k in range(0,length):
        ubviol = x[k] - ubx[k]
        if ubviol > tolerance:
            (tagname,tagstep) = bndtags[k]
            err = (tagstep,'ub',float(ubviol))
            if tagname in violations:
                violations[tagname].append(err)
            else:
                violations[tagname] = [err]
        else:
            lbviol = lbx[k] - x[k]
            if lbviol > tolerance:
                (tagname,tagstep) = bndtags[k]
                err = (tagstep,'lb',float(lbviol))
                if tagname in violations:
                    violations[tagname].append(err)
                else:
                    violations[tagname] = [err]
    return violations

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

        self._bounds = CollMap(self,"bounds")
        self._guess = CollMap(self,"guess")

        self._constraints = Constraints()

        # set (V,XD,XA,U,P)
        self._setupDesignVars()


    def setupCollocation(self,tf):
        if self.collocationIsSetup:
            raise ValueError("you can't setup collocation twice")
        self.collocationIsSetup = True
        
        ## -----------------------------------------------------------------------------
        ## Collocation setup
        ## -----------------------------------------------------------------------------
        # Size of the finite elements
        self.h = tf/self.nk/self.nicp

        # make coefficients for collocation/continuity equations
        self.lagrangePoly = LagrangePoly(deg=self.deg,collPoly=self.collPoly)
        
        # function to get h out
        self.hfun = CS.MXFunction([self._V],[self.h])
        self.hfun.init()
        
        # add collocation constraints
        self._addCollocationConstraints()
        
    def _addCollocationConstraints(self):
        assert hasattr(self,'h')
        assert hasattr(self,'nicp')
        assert hasattr(self,'nk')
        assert hasattr(self,'deg')
        
        assert hasattr(self,'_XD')
        assert hasattr(self,'_XA')
        assert hasattr(self,'_U')
        assert hasattr(self,'_P')

        ffcn = self._makeResidualFun()
 
        ndiff = self.xSize()
        nalg = self.zSize()
        
        # For all finite elements
        for k in range(self.nk):
            for i in range(self.nicp):
                # For all collocation points
                for j in range(1,self.deg+1):
                    # Get an expression for the state derivative at the collocation point
                    xp_jk = 0
                    for j2 in range (self.deg+1):
                        xp_jk += self.lagrangePoly.lDotAtTauRoot[j2][j]*self._XD[k][i][j2] # get the time derivative of the differential states (eq 10.19b)
                    # Add collocation equations to the NLP
                    [fk] = ffcn.call([xp_jk/self.h, self._XD[k][i][j], self._XA[k][i][j-1], self._U[k], self._P])
                    
                    # impose system dynamics (for the differential states (eq 10.19b))
                    self.constrain(fk[:ndiff],'==',0,tag=
                                   "system dynamics, differential states, kIdx: %d,nicpIdx: %d, degIdx: %d" % (k,i,j))
                    
                    # impose system dynamics (for the algebraic states (eq 10.19b))
                    self.constrain(fk[ndiff:],'==',0,tag=
                                   "system dynamics, algebraic states, kIdx: %d,nicpIdx: %d, degIdx: %d" % (k,i,j))
                    
                # Get an expression for the state at the end of the finite element
                xf_k = 0
                for j in range(self.deg+1):
                    xf_k += self.lagrangePoly.lAtOne[j]*self._XD[k][i][j]

                # Add continuity equation to NLP
                if i==self.nicp-1:
                    self.constrain(self._XD[k+1][0][0], '==', xf_k, tag="continuity, kIdx: %d,nicpIdx: %d" % (k,i))
                else:
                    self.constrain(self._XD[k][i+1][0], '==', xf_k, tag="continuity, kIdx: %d,nicpIdx: %d" % (k,i))
                    

    def xVec(self,timestep=None,nicpIdx=0,degIdx=0):
        if not isinstance(timestep, int):
            raise ValueError("timestep needs to be an int")
        if timestep is self.nk:
            assert nicpIdx==0 and degIdx==0,"last timestep is only defined at nicpIdx=0,degIdx=0"
        return self._XD[timestep][nicpIdx][degIdx]
    
    def uVec(self,timestep=None):
        if not isinstance(timestep, int):
            raise ValueError("timestep needs to be an int")
        return self._U[timestep]
    
    def pVec(self):
        return self._P
    
    def zVec(self,timestep=None):
        raise ValueError("zVec not yet safe to use cause it's not defined at the beginning of the interval, need to implement moar inputs to zVec")
        if not isinstance(timestep, int):
            raise ValueError("timestep needs to be an int")
        if timestep==0:
            raise ValueError("there is no algebraic state at timestep 0")
        assert timestep>0
        return self._XA[timestep-1][-1][-1]

    def constrain(self,lhs,comparison,rhs,tag='unnamed_constraint'):
        self._constraints.add(lhs,comparison,rhs,tag)

    def constrainBnds(self,g,(lbg,ubg),tag='unnamed_constraint'):
        self._constraints.addBnds(g,(lbg,ubg),tag)

    def xSize(self):
        return len(self.dae.xNames())
        #return self.dae.xVec().size()
    def zSize(self):
        return len(self.dae.zNames())
#        return self.dae.zVec().size()
    def uSize(self):
        return len(self.dae.uNames())
#        return self.dae.uVec().size()
    def pSize(self):
        return len(self.dae.pNames())
#        return self.dae.pVec().size()

    # Total number of variables
    def getNV(self):
        NXD = self.nicp*self.nk*(self.deg+1)*self.xSize() # Collocated differential states
        NXA = self.nicp*self.nk*self.deg*self.zSize()     # Collocated algebraic states
        NU = self.nk*self.uSize()               # Parametrized controls
        NXF = self.xSize()                 # Final state (only the differential states)
        NP = self.pSize()
        return NXD+NXA+NU+NXF+NP

    def _makeResidualFun(self):
        if not hasattr(self.dae, 'stateDotDummy'):
            raise ValueError("need to set stateDotDummy")
        if not hasattr(self.dae, '_odeRes'):
            raise ValueError("need to set ode residual")
        
        residual = self.dae._odeRes
        
        if hasattr(self.dae,'_algRes'):
            residual = CS.veccat([residual, self.dae._algRes])
        else:
            if self.zSize()>0:
                raise ValueError("you've added algebraic states but haven't set the algebraic residual")
    
        # residual function
        u = self.dae.uVec()
        xd = self.dae.xVec()
        xa = self.dae.zVec()
        xddot = self.dae.stateDotDummy
        p  = self.dae.pVec()
        
        ffcn = CS.SXFunction([xddot,xd,xa,u,p],[residual])
        ffcn.init()
        return ffcn

    def _setupDesignVars(self):
        ## set these:
        # self._V  = V
        # self._XD = XD
        # self._XA = XA
        # self._U  = U
        # self._P  = P
        
        ndiff = self.xSize()
        nalg = self.zSize()
        nu = self.uSize()
        NP = self.pSize()
        nx = ndiff + nalg
        
        # Total number of variables
        NXD = self.nicp*self.nk*(self.deg+1)*ndiff # Collocated differential states
        NXA = self.nicp*self.nk*self.deg*nalg      # Collocated algebraic states
        NU = self.nk*nu                  # Parametrized controls
        NXF = ndiff                 # Final state (only the differential states)
        NV = NXD+NXA+NU+NXF+NP
    
        # NLP variable vector
        V = CS.msym("V",NV)

        # index map
        self._indexMap = CollMap(self,"index map")
        self._dvMap = CollMap(self,"design var map")
        
        offset = 0
        
        # Get the parameters
        P = V[offset:offset+NP]
        for k,name in enumerate(self.dae.pNames()):
            self._indexMap.setVal(name,k)
            self._dvMap.setVal(name,V[offset+k])
        offset += NP
        
        # Get collocated states and parametrized control
        XD = np.resize(np.array([],dtype=CS.MX),(self.nk+1,self.nicp,self.deg+1)) # NB: same name as above
        XA = np.resize(np.array([],dtype=CS.MX),(self.nk,self.nicp,self.deg)) # NB: same name as above
        U = np.resize(np.array([],dtype=CS.MX),self.nk)
        for k in range(self.nk):
            # Collocated states
            for i in range(self.nicp):
                for j in range(self.deg+1):
                    # Get the expression for the state vector
                    XD[k][i][j] = V[offset:offset+ndiff]
                    for m,name in enumerate(self.dae.xNames()):
                        self._indexMap.setVal(name,offset+m,timestep=k,nicpIdx=i,degIdx=j)
                        self._dvMap.setVal(name,V[offset+m],timestep=k,nicpIdx=i,degIdx=j)
                        
                    if j !=0:
                        XA[k][i][j-1] = V[offset+ndiff:offset+ndiff+nalg]
                        for m,name in enumerate(self.dae.zNames()):
                            self._indexMap.setVal(name,offset+ndiff+m,timestep=k,nicpIdx=i,degIdx=j)
                            self._dvMap.setVal(name,V[offset+ndiff+m],timestep=k,nicpIdx=i,degIdx=j)

                    index = (self.deg+1)*(self.nicp*k+i) + j
                    if k==0 and j==0 and i==0:
                        offset += ndiff
                    else:
                        if j!=0:
                            offset += nx
                        else:
                            offset += ndiff
            
            # Parametrized controls
            U[k] = V[offset:offset+nu]
            for m,name in enumerate(self.dae.uNames()):
                self._indexMap.setVal(name,offset+m,timestep=k)
                self._dvMap.setVal(name,V[offset+m],timestep=k)
            offset += nu
        
        # State at end time
        XD[self.nk][0][0] = V[offset:offset+ndiff]
        for m,name in enumerate(self.dae.xNames()):
            self._indexMap.setVal(name,offset+m,timestep=self.nk,degIdx=0,nicpIdx=0)
            self._dvMap.setVal(name,V[offset+m],timestep=self.nk,degIdx=0,nicpIdx=0)
        
        offset += ndiff
        assert offset==NV
    
        self._V  = V
        self._XD = XD
        self._XA = XA
        self._U  = U
        self._P  = P


    def _parseBoundsAndGuess(self,warnZBounds,warnZGuess):
        nk = self.nk
        nicp = self.nicp
        deg = self.deg
        
        # Differential state bounds and initial guess
        xD_min = np.array([[-CS.inf]*self.xSize()]*(nk+1))
        xD_max = np.array([[ CS.inf]*self.xSize()]*(nk+1))
        xD_init = np.array((nk*nicp*(deg+1)+1)*[[None]*self.xSize()]) # needs to be specified for every time interval
        for k,name in enumerate(self.dae.xNames()):
            for timestep in range(nk+1):
                xminmax=self._bounds.lookup(name,timestep=timestep)
                if xminmax is None:
                    raise ValueError('need to set bounds for \"'+name+'\" at timestep '+str(timestep))
                (xmin,xmax) = xminmax
                xD_min[timestep,k] = xmin
                xD_max[timestep,k] = xmax
    
            # linearly interpolate initial guess
            for timestep in range(nk):
                xd0=self._guess.lookup(name,timestep=timestep)
                xd1=self._guess.lookup(name,timestep=timestep+1)
                if xd0 is None:
                    raise ValueError("need to set initial guess for \""+name+ "\" at timestep "+str(timestep))
                if xd1 is None:
                    raise ValueError("need to set initial guess for \""+name+ "\" at timestep "+str(timestep+1))
    
                alpha = 0
                alphaIndex = 0
                for j in range(nicp):
                    for d in range(deg+1):
                        index = (deg+1)*(nicp*timestep+j) + d
                        alpha = alphaIndex/( (deg+1)*nicp - 1.0 )
                        xD_init[index,k] = xd0 + (xd1-xd0)*alpha
                        alphaIndex += 1
            if self._guess.lookup(name,timestep=-1) is None:
                raise ValueError("need to set initial guess for \""+name+ "\" at last timestep")
            xD_init[-1,k] = self._guess.lookup(name,timestep=-1)
        
        # Algebraic state bounds and initial guess
        xA_min = np.array([[-CS.inf]*self.zSize()]*nk)
        xA_max = np.array([[ CS.inf]*self.zSize()]*nk)
        xA_init = np.array((nk*nicp*(deg+1))*[[None]*self.zSize()])
        for k,name in enumerate(self.dae.zNames()):
            for timestep in range(nk):
                zminmax=self._bounds.lookup(name,timestep=timestep,degIdx=1)
                if zminmax is None:
                    if warnZBounds:
                        print "WARNING: didn't set bounds for algebraic state \""+name+"\" at timestep "+str(timestep)+", using (-inf,inf)"
                    zmin = -CS.inf
                    zmax =  CS.inf
                else:
                    (zmin,zmax) = zminmax
    
                xA_min[timestep,k] = zmin
                xA_max[timestep,k] = zmax
                
            # linearly interpolate initial guess
            for timestep in range(nk):
                xa0 = self._guess.lookup(name, timestep=timestep, degIdx=1)
                if timestep<nk-1:
                    xa1 = self._guess.lookup(name, timestep=timestep+1, degIdx=1)
                else:
                    xa1 = self._guess.lookup(name, timestep=timestep, degIdx=1)
    
                if xa0 is None:
                    if warnZGuess:
                        print "WARNING: initial guess for \""+name+ "\" not set at timestep "+str(timestep)+", using guess of 0.0"
                    xa0 = 0.0
                if xa1 is None:
                    # no print statement here to prevent two identical warnings
                    xa1 = 0.0
    
                alpha = 0
                alphaIndex = 0
                for j in range(nicp):
                    for d in range(deg):
                        index = (deg+1)*(nicp*timestep+j) + d
                        alpha = alphaIndex/( deg*nicp - 1.0 )
                        xA_init[index,k] = xa0 + (xa1-xa0)*alpha
                        alphaIndex += 1
        
        # Control bounds and initial guess
        u_min = np.array([[-CS.inf]*self.uSize()]*nk)
        u_max = np.array([[ CS.inf]*self.uSize()]*nk)
        u_init = np.array([[None]*self.uSize()]*nk)
        for k,name in enumerate(self.dae.uNames()):
            for timestep in range(nk):
                uminmax = self._bounds.lookup(name,timestep=timestep)
                if uminmax is None:
                    raise ValueError('need to set bounds for \"'+name+'\" at timestep '+str(timestep))
                (umin,umax) = uminmax
                u_min[timestep,k] = umin
                u_max[timestep,k] = umax
    
                # initial guess
                if self._guess.lookup(name,timestep=timestep) is None:
                    raise ValueError("need to set initial guess for \""+name+ "\" at timestep "+str(timestep))
                u_init[timestep,k] = self._guess.lookup(name,timestep=timestep)
        
        # Parameter bounds and initial guess
        p_min = np.array([-CS.inf]*self.pSize())
        p_max = np.array([ CS.inf]*self.pSize())
        p_init = np.array([None]*self.pSize())
        for k,name in enumerate(self.dae.pNames()):
            pminmax = self._bounds.lookup(name)
            if pminmax is None:
                raise ValueError('need to set bounds for \"'+name+'\"')
            (pmin,pmax) = pminmax
            p_min[k] = pmin
            p_max[k] = pmax
        for k,name in enumerate(self.dae.pNames()):
            if self._guess.lookup(name) is None:
                raise ValueError("need to set initial guess for \""+name+ "\"")
            p_init[k] = self._guess.lookup(name)
    
        return {'xD_init':xD_init, 'xD_min':xD_min, 'xD_max':xD_max,
                'xA_init':xA_init, 'xA_min':xA_min, 'xA_max':xA_max,
                 'u_init': u_init,  'u_min': u_min,  'u_max': u_max,
                 'p_init': p_init,  'p_min': p_min,  'p_max': p_max}

    def _vectorizeBoundsAndGuess(self, boundsAndGuess):
        p_init  = boundsAndGuess['p_init']
        p_min   = boundsAndGuess['p_min']
        p_max   = boundsAndGuess['p_max']
        xD_init = boundsAndGuess['xD_init']
        xD_min  = boundsAndGuess['xD_min']
        xD_max  = boundsAndGuess['xD_max']
        xA_init = boundsAndGuess['xA_init']
        xA_min  = boundsAndGuess['xA_min']
        xA_max  = boundsAndGuess['xA_max']
        u_init  = boundsAndGuess['u_init']
        u_min   = boundsAndGuess['u_min']
        u_max   = boundsAndGuess['u_max']
        
        nicp = self.nicp
        nk = self.nk
        deg = self.deg
        
        XD = self._XD
        XA = self._XA
        U = self._U
        P = self._P
        V = self._V
        
        # -----------------------------------------------------------------------------
        # NLP setup
        # -----------------------------------------------------------------------------
        # Dimensions of the problem
        ndiff = self.xSize() # number of differential states
        nalg = self.zSize()  # number of algebraic states
        nu   = self.uSize()  # number of controls
        NP   = self.pSize()  # number of parameters
        nx   = ndiff + nalg  # total number of states
        NV = self.getNV()

        # tags for the bounds
        bndtags = []
        
        # All variables with bounds and initial guess
        vars_lb = np.zeros(NV)
        vars_ub = np.zeros(NV)
        vars_init = np.zeros(NV)
        offset = 0
        
        # Get the parameters
        vars_init[offset:offset+NP] = p_init
        vars_lb[offset:offset+NP] = p_min
        vars_ub[offset:offset+NP] = p_max
        for name in self.dae.pNames():
            bndtags.append((name,None))
        offset += NP

        # Get collocated states and parametrized control
        for k in range(nk):  
            # Collocated states
            for i in range(nicp):
                for j in range(deg+1):
                    # Add the initial condition
                    index = (deg+1)*(nicp*k+i) + j
                    if k==0 and j==0 and i==0:
                        vars_init[offset:offset+ndiff] = xD_init[index,:]
                        
                        vars_lb[offset:offset+ndiff] = xD_min[k,:]
                        vars_ub[offset:offset+ndiff] = xD_max[k,:]

                        for name in self.dae.xNames():
                            bndtags.append((name,(k,i,j)))
                        
                        offset += ndiff
                    else:
                        if j!=0:
                            vars_init[offset:offset+nx] = np.concatenate((xD_init[index,:],xA_init[index-1,:]))
    
                            vars_lb[offset:offset+nx] = np.concatenate((np.minimum(xD_min[k,:],xD_min[k+1,:]),xA_min[k,:]))
                            vars_ub[offset:offset+nx] = np.concatenate((np.maximum(xD_max[k,:],xD_max[k+1,:]),xA_max[k,:]))
                            for name in self.dae.xNames()+self.dae.zNames():
                                bndtags.append((name,(k,i,j)))
                            
                            offset += nx
                        else:
                            vars_init[offset:offset+ndiff] = xD_init[index,:]
    
                            for name in self.dae.xNames():
                                bndtags.append((name,(k,i,j)))
                            
                            if i==0:
                                vars_lb[offset:offset+ndiff] = xD_min[k,:]
                                vars_ub[offset:offset+ndiff] = xD_max[k,:]
                            else:
                                vars_lb[offset:offset+ndiff] = np.minimum(xD_min[k,:],xD_min[k+1,:])
                                vars_ub[offset:offset+ndiff] = np.maximum(xD_max[k,:],xD_max[k+1,:])
                                
                            offset += ndiff
            
            # Parametrized controls
            vars_lb[offset:offset+nu] = u_min[k]
            vars_ub[offset:offset+nu] = u_max[k]
            vars_init[offset:offset+nu] = u_init[k]
            for name in self.dae.uNames():
                bndtags.append((name,(k,i,j)))
            offset += nu
        
        # State at end time
        vars_lb[offset:offset+ndiff] = xD_min[nk,:]
        vars_ub[offset:offset+ndiff] = xD_max[nk,:]
        vars_init[offset:offset+ndiff] = xD_init[-1,:]
        for name in self.dae.xNames():
            bndtags.append((name,(nk,0,0)))
        offset += ndiff
        assert offset==NV
        self.bndtags = np.array(bndtags)
        return (vars_init,vars_lb,vars_ub)

    def interpolateInitialGuess(self,filename,force=False,quiet=False):
        print "interpolating initial guess..."
        f=open(filename,'r')
        traj = pickle.load(f)
        f.close()

        assert isinstance(traj,TrajectoryData), "the file \""+filename+"\" doean't have a pickled TrajectoryData"

        if traj.collMap._nicp == 1:
            h = traj.tgrid[1][0][0] - traj.tgrid[0][0][0]
        else:
            h = traj.tgrid[0][1][0] - traj.tgrid[0][0][0]
        h *= (traj.collMap._nk-1)*traj.collMap._nicp/float((self.nk-1)*self.nicp)
            
        missing = set(self.dae.xNames()+self.dae.zNames()+self.dae.uNames()+self.dae.pNames())

        # interpolate differential states
        for name in traj.collMap._xNames:
            if name not in missing:
                continue
            missing.remove(name)

            # make piecewise poly
            pp = None
            for timestepIdx in range(traj.collMap._nk):
                for nicpIdx in range(traj.collMap._nicp):
                    ts = []
                    ys = []
                    for degIdx in range(traj.collMap._deg+1):
                        ts.append(traj.tgrid[timestepIdx,nicpIdx,degIdx])
                        ys.append([traj.collMap.lookup(name,timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx)])
                    if pp is None:
                        pp = PiecewisePolynomial(ts,ys)
                    else:
                        pp.extend(ts,ys)
            pp.extend([traj.tgrid[-1,0,0]],[[traj.collMap.lookup(name,timestep=-1,nicpIdx=0,degIdx=0)]])

            # evaluate piecewise poly to set initial guess
            t0 = 0.0
            for timestepIdx in range(self.nk):
                for nicpIdx in range(self.nicp):
                    for degIdx in range(self.deg+1):
                        time = t0 + h*self.lagrangePoly.tau_root[degIdx]
                        self.guess(name,pp(time),timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx,force=force,quiet=quiet)
                    t0 += h
            self.guess(name,pp(t0),timestep=-1,nicpIdx=0,degIdx=0,force=force,quiet=quiet)

        # interpolate algebraic variables
        for name in traj.collMap._zNames:
            if name not in missing:
                continue
            missing.remove(name)

            # make piecewise poly
            pp = None
            for timestepIdx in range(traj.collMap._nk):
                for nicpIdx in range(traj.collMap._nicp):
                    ts = []
                    ys = []
                    for degIdx in range(1,traj.collMap._deg+1):
                        ts.append(traj.tgrid[timestepIdx,nicpIdx,degIdx])
                        ys.append([traj.collMap.lookup(name,timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx)])
                    if pp is None:
                        pp = PiecewisePolynomial(ts,ys)
                    else:
                        pp.extend(ts,ys)

            # evaluate piecewise poly to set initial guess
            t0 = 0.0
            for timestepIdx in range(self.nk):
                for nicpIdx in range(self.nicp):
                    for degIdx in range(1,self.deg+1):
                        time = t0 + h*self.lagrangePoly.tau_root[degIdx]
                        self.guess(name,pp(time),timestep=timestepIdx,nicpIdx=nicpIdx,degIdx=degIdx,force=force,quiet=quiet)
                    t0 += h
                
        # interpolate controls
        for name in traj.collMap._uNames:
            if name not in missing:
                continue
            missing.remove(name)

            # make piecewise poly
            ts = []
            ys = []
            for timestepIdx in range(traj.collMap._nk):
                ts.append(traj.tgrid[timestepIdx,0,0])
                ys.append([traj.collMap.lookup(name,timestep=timestepIdx)])
            pp = PiecewisePolynomial(ts,ys)

            # evaluate piecewise poly to set initial guess
            t0 = 0.0
            for timestepIdx in range(self.nk):
                self.guess(name,pp(t0),timestep=timestepIdx,force=force,quiet=quiet)
                t0 += h
        
        # set parameters
        for name in traj.collMap._pNames:
            if name not in missing:
                continue
            missing.remove(name)
            self.guess(name,traj.collMap.lookup(name),force=force,quiet=quiet)

        msg = "finished interpolating initial guess"
        if len(missing) > 0:
            msg += ", couldn't find fields: "+str(list(missing))
        else:
            msg += ", all fields found"
        print msg

    def setupSolver(self,solverOpts=[],constraintFunOpts=[],callback=None):
        if not self.collocationIsSetup:
            raise ValueError("you forgot to call setupCollocation")
        
        g =   self._constraints.getG()
        lbg = self._constraints.getLb()
        ubg = self._constraints.getUb()

        # Nonlinear constraint function
        gfcn = CS.MXFunction([self._V],[g])
        setFXOptions(gfcn,constraintFunOpts)
        
        # Objective function of the NLP
        if not hasattr(self,'_objective'):
            raise ValueError('need to set objective function')
        ofcn = CS.MXFunction([self._V],[self._objective])

        # solver callback (optional)
        if callback is not None:
            nd = self._V.size()
            nc = self._constraints.getG().size()
    
            c = CS.PyFunction( callback, CS.nlpsolverOut(x_opt=CS.sp_dense(nd,1), cost=CS.sp_dense(1,1), lambda_x=CS.sp_dense(nd,1), lambda_g = CS.sp_dense(nc,1), g = CS.sp_dense(nc,1) ), [CS.sp_dense(1,1)] )
            c.init()
            solverOpts.append( ("iteration_callback", c) )

        # Allocate an NLP solver
        self.solver = CS.IpoptSolver(ofcn,gfcn)
        
        # Set options
        setFXOptions(self.solver, solverOpts)
        
        # initialize the solver
        self.solver.init()
        
        # Bounds on g
        self.solver.setInput(lbg,CS.NLP_LBG)
        self.solver.setInput(ubg,CS.NLP_UBG)

    def solve(self,xInit=None,warnZBounds=False,warnZGuess=False):
        (vars_init,vars_lb,vars_ub) = self._vectorizeBoundsAndGuess( self._parseBoundsAndGuess(warnZBounds,warnZGuess) )
        if xInit is not None:
            vars_init = xInit

        # Initial condition
        self.solver.setInput(vars_init, CS.NLP_X_INIT)
        
        # Bounds on x
        self.solver.setInput(vars_lb,CS.NLP_LBX)
        self.solver.setInput(vars_ub,CS.NLP_UBX)
        
        # Solve the problem
        self.solver.solve()
        
        # Print the optimal cost
        print "optimal cost: ", float(self.solver.output(CS.NLP_COST))
        
        # Retrieve the solution
        v_opt = np.array(self.solver.output(CS.NLP_X_OPT))
        
        return self.devectorize(v_opt)
    

    def devectorize(self,v_opt):
        return CollMap(self,'devectorized design vars',devectorize=v_opt)

    def mkTimeGrid(self,v_opt):
        self.hfun.setInput(v_opt)
        self.hfun.evaluate()
        h = float(self.hfun.output())

        tgrid = np.resize([],(self.nk+1,self.nicp,self.deg+1))
        tf = 0.0
        for k in range(self.nk):
            for i in range(self.nicp):
                tgrid[k,i,:] = tf + h*np.array(self.lagrangePoly.tau_root)
                tf += h
        tgrid[self.nk,0,0] = tf
        return tgrid

    def mkTimeGridVec(self,v_opt):
        self.hfun.setInput(v_opt)
        self.hfun.evaluate()
        h = float(self.hfun.output())
        tg = np.array(self.lagrangePoly.tau_root)*h
        for k in range(self.nk*self.nicp):
            if k == 0:
                tgrid = tg
            else:
                tgrid = np.append(tgrid,tgrid[-1]+tg)
        tgrid = np.append(tgrid,tgrid[-1])
        return tgrid

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
        
        assert isinstance(timestep,int)
        self._bounds.setVal(name,val,timestep=timestep,quiet=quiet)

    def guess(self,name,val,timestep=None,nicpIdx=None,degIdx=None,quiet=False,force=False):
        assert isinstance(name,str)
        assert isinstance(val,numbers.Real)

        # handle timestep == None
        if timestep is None:
            assert (nicpIdx is None)
            assert (degIdx is None)
            if name in self._guess._xMap:
                for k in range(self.nk+1):
                    self.guess(name,val,timestep=k,quiet=quiet)
                return
            if name in self._guess._zMap:
                for k in range(self.nk):
                    self.guess(name,val,timestep=k,quiet=quiet)
                return
            if name in self._guess._uMap:
                for k in range(self.nk):
                    self.guess(name,val,timestep=k,quiet=quiet)
                return
            if name in self._guess._pMap:
                self._guess.setVal(name,val,force=force)
                return
        
        assert isinstance(timestep,int)
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
        

    def lookup(self,name,timestep=None,nicpIdx=0,degIdx=0):
        assert isinstance(name,str)
        assert nicpIdx<self.nicp,'nicpIdx must be < nicp'
        assert degIdx<self.deg+1,'degIdx must be < deg+1'
        assert degIdx>=0,'degIdx must be positive'
        
        if name in self.dae.pNames():
            if timestep is not None:
                raise ValueError("can't lookup a parameter ("+name+")at a certain timestep")
            k = self.dae.pNames().index(name)
            return self._P[k]

        if timestep is None:
            if name in self.dae.xNames()+self.dae.zNames():
                return CS.veccat([self.lookup(name,k,nicpIdx,degIdx) for k in range(self.nk+1)])
            elif name in self.dae.uNames():
                return CS.veccat([self.lookup(name,k,nicpIdx,degIdx) for k in range(self.nk)])
            else:
                raise ValueError("unrecognized key \""+name+"\"")

        assert isinstance(timestep,int)
        if name in self.dae.xNames():
            k = self.dae.xNames().index(name)
            return self._XD[timestep][nicpIdx][degIdx][k]
        if name in self.dae.zNames():
            assert degIdx!=0,'there is no algebraic state at degIdx=0'
            k = self.dae.zNames().index(name)
            return self._XA[timestep][nicpIdx][degIdx-1][k]
        if name in self.dae.uNames():
            k = self.dae.uNames().index(name)
            return self._U[timestep][k]

        raise ValueError("unrecognized variable \""+name+"\"")

    def setObjective(self,obj):
        if hasattr(self,'_objective'):
            raise ValueError("can't change objective function once it's already set")
        self._objective = obj


if __name__=='__main__':
    dae = Dae()
    # -----------------------------------------------------------------------------
    # Model setup
    # -----------------------------------------------------------------------------
    # Declare variables (use scalar graph)
    t  = CS.ssym("t")          # time
    u  = dae.addU("u")          # control
    xd = CS.veccat(dae.addX(["x0","x1","x2"]))  # differential state
#    xa = CS.ssym("xa",0,1)     # algebraic state
#    xddot = CS.ssym("xdot",3)  # differential state time derivative
#    p     = CS.ssym("p",0,1)      # parameters
    dae.addOutput('u^2',u*u)
    dae.stateDotDummy = CS.veccat( [CS.ssym(name+"DotDummy") for name in dae._xNames] )

    # ODE right hand side function
    rhs = CS.vertcat([(1 - xd[1]*xd[1])*xd[0] - xd[1] + u, \
           xd[0], \
           xd[0]*xd[0] + xd[1]*xd[1] + u*u])
    dae.setOdeRes( rhs - dae.stateDotDummy )

    nicp_ = 1        # Number of (intermediate) collocation points per control interval
    nk_ = 10         # Control discretization
    deg_ = 4         # Degree of interpolating polynomial
    collPoly_ = 'RADAU'  # use Radau collocation points
    #collPoly = 'LEGENDRE'    # use LEGENDRE collocation points

    coll = Coll(dae, nk_, nicp_, deg_, collPoly_)
    coll.bound('x0',(-CS.inf,CS.inf))
    coll.bound('x1',(-CS.inf,CS.inf))
    coll.bound('x2',(-CS.inf,CS.inf))

    coll.bound('x0',(0,0),timestep=0)
    coll.bound('x1',(1,1),timestep=0)
    coll.bound('x2',(0,0),timestep=0)
    
    coll.bound('x0',(0,0),timestep=-1)
    coll.bound('x1',(0,0),timestep=-1)

    coll.bound('u',(-0.75,1))

    coll.guess('x0',0)
    coll.guess('x1',0)
    coll.guess('x2',0)
    coll.guess('u',0)

    coll.setObjective( coll.lookup('x2',timestep=-1) )

    coll.setupCollocation(2.5)
    coll.setupSolver()
    opt = coll.solve()
