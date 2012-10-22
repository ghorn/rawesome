import casadi as CS
import numpy as np
import matplotlib.pyplot as plt
import numbers
from ocputils import Constraints,setFXOptions

from dae import Dae

def mkCollocationPoints():
    # Legendre collocation points
    legendre_points1 = [0,0.500000]
    legendre_points2 = [0,0.211325,0.788675]
    legendre_points3 = [0,0.112702,0.500000,0.887298]
    legendre_points4 = [0,0.069432,0.330009,0.669991,0.930568]
    legendre_points5 = [0,0.046910,0.230765,0.500000,0.769235,0.953090]
    legendre_points = [0,legendre_points1,legendre_points2,legendre_points3,legendre_points4,legendre_points5]
    
    # Radau collocation points
    radau_points1 = [0,1.000000]
    radau_points2 = [0,0.333333,1.000000]
    radau_points3 = [0,0.155051,0.644949,1.000000]
    radau_points4 = [0,0.088588,0.409467,0.787659,1.000000]
    radau_points5 = [0,0.057104,0.276843,0.583590,0.860240,1.000000]
    radau_points = [0,radau_points1,radau_points2,radau_points3,radau_points4,radau_points5]
    
    # Type of collocation points
    return {'LEGENDRE':legendre_points, 'RADAU':radau_points}

def setupCoeffs(deg=None,collPoly=None,nk=None,h=None):
    assert(deg is not None)
    assert(collPoly is not None)
    assert(nk is not None)
    assert(h is not None)
    
    collocation_points = mkCollocationPoints()
    assert( collPoly in collocation_points )
  
    # Coefficients of the collocation equation
    C = np.zeros((deg+1,deg+1))
    # Coefficients of the continuity equation
    D = np.zeros(deg+1)
    
    # Collocation point
    tau = CS.ssym("__collocation_tau__")
      
    # All collocation time points
    tau_root = collocation_points[collPoly][deg]
#    T = np.zeros((nk,deg+1))
#    for i in range(nk):
#      for j in range(deg+1):
#    	T[i][j] = h*(i + tau_root[j])
    
    # For all collocation points: eq 10.4 or 10.17 in Biegler's book
    # Construct Lagrange polynomials to get the polynomial basis at the collocation point
    for j in range(deg+1):
        L = 1
        for j2 in range(deg+1):
            if j2 != j:
                L *= (tau-tau_root[j2])/(tau_root[j]-tau_root[j2])
        lfcn = CS.SXFunction([tau],[L])
        lfcn.init()
        # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
        lfcn.setInput(1.0)
        lfcn.evaluate()
        D[j] = lfcn.output()
        # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
        for j2 in range(deg+1):
            lfcn.setInput(tau_root[j2])
            lfcn.setFwdSeed(1.0)
            lfcn.evaluate(1,0)
            C[j][j2] = lfcn.fwdSens()
  
    return (C,D,tau,tau_root)


def setupCollocationConstraints(coll,C,ffcn,h,D):
    nicp = coll.nicp
    nk = coll.nk
    deg = coll.deg
    
    XD = coll._XD
    XA = coll._XA
    U = coll._U
    P = coll._P

    ndiff = coll.xSize()
    nalg = coll.zSize()
    
    # Constraint function for the NLP
    g = []
    lbg = []
    ubg = []
    
    # For all finite elements
    for k in range(nk):
        for i in range(nicp):
            # For all collocation points
            for j in range(1,deg+1):   		
                # Get an expression for the state derivative at the collocation point
                xp_jk = 0
                for j2 in range (deg+1):
                    xp_jk += C[j2][j]*XD[k][i][j2] # get the time derivative of the differential states (eq 10.19b)
                
                # Add collocation equations to the NLP
                [fk] = ffcn.call([xp_jk/h, XD[k][i][j], XA[k][i][j-1], U[k], P])
                g += [fk[:ndiff]]           # impose system dynamics (for the differential states (eq 10.19b))
                lbg.append(np.zeros(ndiff)) # equality constraints
                ubg.append(np.zeros(ndiff)) # equality constraints
                
                g += [fk[ndiff:]]          # impose system dynamics (for the algebraic states (eq 10.19b))
                lbg.append(np.zeros(nalg)) # equality constraints
                ubg.append(np.zeros(nalg)) # equality constraints
                
            # Get an expression for the state at the end of the finite element
            xf_k = 0
            for j in range(deg+1):
                xf_k += D[j]*XD[k][i][j]
                
            # Add continuity equation to NLP
            if i==nicp-1:
    #            print "a ", k, i
                g += [XD[k+1][0][0] - xf_k]
            else:
    #            print "b ", k, i
                g += [XD[k][i+1][0] - xf_k]
            
            lbg.append(np.zeros(ndiff))
            ubg.append(np.zeros(ndiff))
    
    return (g,lbg,ubg)

class Coll():
    collocationIsSetup = False
    def __init__(self, dae, nk=None, nicp=1, deg=4, collPoly='RADAU'):
        assert(nk is not None)
        assert(isinstance(dae, Dae))
        
        print "collocation init, yay"
        self.dae = dae
        
        self.nk = nk
        self.nicp = nicp
        self.deg = deg
        self.collPoly = collPoly

        self._xBounds = dict( [(name,[None]*(nk+1)) for name in dae.xNames()] )
        self._zBounds = dict( [(name,[None]*nk)     for name in dae.zNames()] )
        self._uBounds = dict( [(name,[None]*nk)     for name in dae.uNames()] )
        self._pBounds = dict( [(name,None)          for name in dae.pNames()] )

        self._xGuess = dict( [(name,[None]*(nk+1)) for name in dae.xNames()] )
        self._zGuess = dict( [(name,[None]*nk)     for name in dae.zNames()] )
        self._uGuess = dict( [(name,[None]*nk)     for name in dae.uNames()] )
        self._pGuess = dict( [(name,None)          for name in dae.pNames()] )

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
        C,D,tau,tau_root = setupCoeffs(deg=self.deg,collPoly=self.collPoly,nk=self.nk,h=self.h)
        self.tau_root = tau_root
        self.tau = tau
        
        ffcn = self.makeResidualFun()

        # function to get h out
        self.hfun = CS.MXFunction([self._V],[self.h])
        self.hfun.init()
        
        # add collocaiton constraints
        (g_coll,lbg_coll,ubg_coll) = setupCollocationConstraints(self,C,ffcn,self.h,D)

        assert(len(g_coll)==len(lbg_coll) and len(g_coll)==len(ubg_coll))
        for k in range(len(g_coll)):
            assert(lbg_coll[k].size == ubg_coll[k].size)
            if lbg_coll[k].size>0:
                assert(g_coll[k].size()==lbg_coll[k].size)
                self.constrainBnds(g_coll[k],(lbg_coll[k],ubg_coll[k]))
        

    def xVec(self,timestep=None):
        if not isinstance(timestep, int):
            raise ValueError("timestep needs to be an int")
        return self._XD[timestep][0][0]
    
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
        assert(timestep>0)
        return self._XA[timestep-1][-1][-1]

    def constrain(self,lhs,comparison,rhs):
        self._constraints.add(lhs,comparison,rhs)

    def constrainBnds(self,g,(lbg,ubg)):
        self._constraints.addBnds(g,(lbg,ubg))

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

    def makeResidualFun(self):
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
        
        ndiff = len(self.dae.xNames())
        nalg = len(self.dae.zNames())
        nu = len(self.dae.uNames())
        NP = len(self.dae.pNames())
        nx = ndiff + nalg
        
        # Total number of variables
        NXD = self.nicp*self.nk*(self.deg+1)*ndiff # Collocated differential states
        NXA = self.nicp*self.nk*self.deg*nalg      # Collocated algebraic states
        NU = self.nk*nu                  # Parametrized controls
        NXF = ndiff                 # Final state (only the differential states)
        NV = NXD+NXA+NU+NXF+NP
    
        # NLP variable vector
        V = CS.msym("V",NV)
        
        offset = 0
        
        # Get the parameters
        P = V[offset:offset+NP]
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
                    if j !=0:
                        XA[k][i][j-1] = V[offset+ndiff:offset+ndiff+nalg]
                    # Add the initial condition
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
            offset += nu
        
        # State at end time
        XD[self.nk][0][0] = V[offset:offset+ndiff]
        
        offset += ndiff
        assert(offset==NV)
    
        self._V  = V
        self._XD = XD
        self._XA = XA
        self._U  = U
        self._P  = P


    def _parseBoundsAndGuess(self):
        nk = self.nk
        nicp = self.nicp
        deg = self.deg
        
        # Differential state bounds and initial guess
        xD_min = np.array([[-CS.inf]*self.xSize()]*(nk+1))
        xD_max = np.array([[ CS.inf]*self.xSize()]*(nk+1))
        xD_init = np.array((nk*nicp*(deg+1)+1)*[[None]*self.xSize()]) # needs to be specified for every time interval
        for k,name in enumerate(self.dae.xNames()):
            for timestep in range(nk+1):
                xminmax = self._xBounds[name][timestep]
                if xminmax is None:
                    raise ValueError('need to set bounds for \"'+name+'\" at timestep '+str(timestep))
                (xmin,xmax) = xminmax
                xD_min[timestep,k] = xmin
                xD_max[timestep,k] = xmax
    
            # linearly interpolate initial guess
            for timestep in range(nk):
                xd0 = self._xGuess[name][timestep]
                xd1 = self._xGuess[name][timestep+1]
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
            if self._xGuess[name][-1] is None:
                raise ValueError("need to set initial guess for \""+name+ "\" at last timestep")
            xD_init[-1,k] = self._xGuess[name][-1]
        
        # Algebraic state bounds and initial guess
        xA_min = np.array([[-CS.inf]*self.zSize()]*nk)
        xA_max = np.array([[ CS.inf]*self.zSize()]*nk)
        xA_init = np.array((nk*nicp*(deg+1))*[[None]*self.zSize()])
        for k,name in enumerate(self.dae.zNames()):
            for timestep in range(nk):
                zminmax = self._zBounds[name][timestep]
                if zminmax is None:
                    print "WARNING: didn't set bounds for algebraic state \""+name+"\" at timestep "+str(timestep)+", using (-inf,inf)"
                    zmin = -CS.inf
                    zmax =  CS.inf
                else:
                    (zmin,zmax) = zminmax
    
                xA_min[timestep,k] = zmin
                xA_max[timestep,k] = zmax
                
            # linearly interpolate initial guess
            for timestep in range(nk):
                xa0 = self._zGuess[name][timestep]
                if timestep<nk-1:
                    xa1 = self._zGuess[name][k+1]
                else:
                    xa1 = self._zGuess[name][k]
    
                if xa0 is None:
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
        u_init = np.array([[None]*self.uSize()]*nk) # needs to be specified for every time interval (even though it stays constant)
        for k,name in enumerate(self.dae.uNames()):
            for timestep in range(nk):
                uminmax = self._uBounds[name][timestep]
                if uminmax is None:
                    raise ValueError('need to set bounds for \"'+name+'\" at timestep '+str(timestep))
                (umin,umax) = uminmax
                u_min[timestep,k] = umin
                u_max[timestep,k] = umax
    
                # initial guess
                if self._uGuess[name][timestep] is None:
                    raise ValueError("need to set initial guess for \""+name+ "\" at timestep "+str(timestep))
                u_init[timestep,k] = self._uGuess[name][timestep]
        
        # Parameter bounds and initial guess
        p_min = np.array([-CS.inf]*self.pSize())
        p_max = np.array([ CS.inf]*self.pSize())
        p_init = np.array([None]*self.pSize())
        for k,name in enumerate(self.dae.pNames()):
            pminmax = self._pBounds[name]
            if pminmax is None:
                raise ValueError('need to set bounds for \"'+name+'\"')
            (pmin,pmax) = pminmax
            p_min[k] = pmin
            p_max[k] = pmax
        for k,name in enumerate(self.dae.pNames()):
            if self._pGuess[name] is None:
                raise ValueError("need to set initial guess for \""+name+ "\"")
            p_init[k] = self._pGuess[name]
    
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
        
        # All variables with bounds and initial guess
        vars_lb = np.zeros(NV)
        vars_ub = np.zeros(NV)
        vars_init = np.zeros(NV)
        offset = 0
        
        # Get the parameters
        vars_init[offset:offset+NP] = p_init
        vars_lb[offset:offset+NP] = p_min
        vars_ub[offset:offset+NP] = p_max
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
                        offset += ndiff
                    else:
                        if j!=0:
                            vars_init[offset:offset+nx] = np.concatenate((xD_init[index,:],xA_init[index-1,:]))
    
                            vars_lb[offset:offset+nx] = np.concatenate((np.minimum(xD_min[k,:],xD_min[k+1,:]),xA_min[k,:]))
                            vars_ub[offset:offset+nx] = np.concatenate((np.maximum(xD_max[k,:],xD_max[k+1,:]),xA_max[k,:]))
                            offset += nx
                        else:
                            vars_init[offset:offset+ndiff] = xD_init[index,:]
    
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
            offset += nu
        
        # State at end time
        vars_lb[offset:offset+ndiff] = xD_min[nk,:]
        vars_ub[offset:offset+ndiff] = xD_max[nk,:]
        vars_init[offset:offset+ndiff] = xD_init[-1,:]
        offset += ndiff
        assert(offset==NV)
        return (vars_init,vars_lb,vars_ub)

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

    def solve(self,xInit=None):
        (vars_init,vars_lb,vars_ub) = self._vectorizeBoundsAndGuess( self._parseBoundsAndGuess() )
            
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
        deg = self.deg
        nicp = self.nicp
        nk = self.nk

        ndiff = len(self.dae.xNames())
        nalg = len(self.dae.zNames())
        nu = len(self.dae.uNames())
        NP = len(self.dae.pNames())
        
        ## ----
        ## RETRIEVE THE SOLUTION
        ## ----
        xD_opt = np.resize(np.array([]),(ndiff,(deg+1)*nicp*(nk)+1))
        xA_opt = np.resize(np.array([]),(nalg,(deg)*nicp*(nk)))
        u_opt = np.resize(np.array([]),(nu,(deg+1)*nicp*(nk)+1))
        offset = 0
        offset2 = 0
        offset3 = 0
        offset4 = 0
    
        p_opt = v_opt[offset:offset+NP]
        offset += NP
        
        for k in range(nk):  
            for i in range(nicp):
                for j in range(deg+1):
                    xD_opt[:,offset2] = v_opt[offset:offset+ndiff][:,0]
                    offset2 += 1
                    offset += ndiff
                    if j!=0:
                        xA_opt[:,offset4] = v_opt[offset:offset+nalg][:,0]
                        offset4 += 1
                        offset += nalg
            utemp = v_opt[offset:offset+nu][:,0]
            for i in range(nicp):
                for j in range(deg+1):
                    u_opt[:,offset3] = utemp
                    offset3 += 1
            #    u_opt += v_opt[offset:offset+nu]
            offset += nu
            
        xD_opt[:,-1] = v_opt[offset:offset+ndiff][:,0]
    
    
        # The algebraic states are not defined at the first collocation point of the finite elements:
        # with the polynomials we compute them at that point
        Da = np.zeros(deg)
        for j in range(1,deg+1):
            # Lagrange polynomials for the algebraic states: exclude the first point
            La = 1
            for j2 in range(1,deg+1):
                if j2 != j:
                    La *= (self.tau-self.tau_root[j2])/(self.tau_root[j]-self.tau_root[j2])
            lafcn = CS.SXFunction([self.tau],[La])
            lafcn.init()
            lafcn.setInput(self.tau_root[0])
            lafcn.evaluate()
            Da[j-1] = lafcn.output()
        
        xA_plt = np.resize(np.array([]),(nalg,(deg+1)*nicp*(nk)+1))
        offset4=0
        offset5=0
        for k in range(nk):  
            for i in range(nicp):
                for j in range(deg+1):
                    if j!=0:         
                        xA_plt[:,offset5] = xA_opt[:,offset4]
                        offset4 += 1
                        offset5 += 1
                    else:
                        xa0 = 0
                        for j in range(deg):
                            xa0 += Da[j]*xA_opt[:,offset4+j]
                        xA_plt[:,offset5] = xa0
                        #xA_plt[:,offset5] = xA_opt[:,offset4]
                        offset5 += 1
        
        xA_plt[:,-1] = xA_plt[:,-2]    

        self.hfun.setInput(v_opt)
        self.hfun.evaluate()
        h = float(self.hfun.output())
        tg = np.array(self.tau_root)*h
        for k in range(self.nk*self.nicp):
            if k == 0:
                tgrid = tg
            else:
                tgrid = np.append(tgrid,tgrid[-1]+tg)
        tgrid = np.append(tgrid,tgrid[-1])

        # make nice dict of outputs
        vardict = {}
        for k,name in enumerate(self.dae.pNames()):
            vardict[name] = float(v_opt[k])
        for k,name in enumerate(self.dae.xNames()):
            vardict[name] = xD_opt[k,:]
        for k,name in enumerate(self.dae.zNames()):
            vardict[name] = xA_plt[k,:]
        for k,name in enumerate(self.dae.uNames()):
            vardict[name] = u_opt[k,:]
            
        return {'vardict':vardict, 'tgrid':tgrid, 'x':xD_opt, 'u':u_opt, 'z':xA_opt, 'zPlt':xA_plt, 'p':v_opt[:NP], 'X_OPT':v_opt}


    def bound(self,name,val,timestep=None,quiet=False,force=False):
        assert(isinstance(name,str))
        assert(isinstance(val,tuple))
        assert(len(val)==2)
        assert(isinstance(val[0],numbers.Real))
        assert(isinstance(val[1],numbers.Real))

        if not any([name in d for d in [self._xBounds, self._zBounds, self._uBounds, self._pBounds]]):
            raise KeyError("Unrecognized variable \""+name+"\"")

        # handle timestep == None
        if timestep is None:
            if name in self._xBounds:
                for k in range(self.nk+1):
                    self.bound(name,val,timestep=k,quiet=quiet)
                return
            if name in self._zBounds:
                for k in range(self.nk):
                    self.bound(name,val,timestep=k,quiet=quiet)
                return
            if name in self._uBounds:
                for k in range(self.nk):
                    self.bound(name,val,timestep=k,quiet=quiet)
                return
            if name in self._pBounds:
                oldbound = self._pBounds[name]
                if ((self._pBounds[name] is not None) and (not force)):
                    raise ValueError("can't change bound on \""+name+"\" once it's set unless you use force=True (tried to change bound "+str(oldbound)+" to "+str(val))
                self._pBounds[name] = val
                return
        
        assert(isinstance(timestep,int))
        if name in self._pBounds:
            raise ValueError("can't set bounds on a parameter on a specific timestep")
        if name in self._xBounds:
            if timestep>self.nk+1:
                raise ValueError("differential state \""+name+"\" bound timestep too high ("+str(timestep)+">"+str(self.nk+1))
            oldbound = self._xBounds[name][timestep]
            if (oldbound is not None) and (quiet is False):
                print "WARNING: changing \"%s\" bounds at timestep %d from %s to %s" % (name,timestep,str(oldbound),str(val))
            self._xBounds[name][timestep] = val
            return
        if name in self._zBounds:
            if timestep>self.nk:
                raise ValueError("algebraic state \""+name+"\" bound timestep too high ("+str(timestep)+">"+str(self.nk))
            oldbound = self._zBounds[name][timestep]
            if (oldbound is not None) and (quiet is False):
                print "WARNING: changing \"%s\" bounds at timestep %d from %s to %s" % (name,timestep,str(oldbound),str(val))
            self._zBounds[name][timestep] = val
            return
        if name in self._uBounds:
            if timestep>self.nk:
                raise ValueError("action \""+name+"\" bound timestep too high ("+str(timestep)+">"+str(self.nk))
            oldbound = self._uBounds[name][timestep]
            if (oldbound is not None) and (quiet is False):
                print "WARNING: changing \"%s\" bounds at timestep %d from %s to %s" % (name,timestep,str(oldbound),str(val))
            self._uBounds[name][timestep] = val
            return

    def guess(self,name,val,timestep=None,quiet=False):
        assert(isinstance(name,str))
        assert(isinstance(val,numbers.Real))

        if not any([name in d for d in [self._xBounds, self._zBounds, self._uBounds, self._pBounds]]):
            raise KeyError("Unrecognized variable \""+name+"\"")

        # handle timestep == None
        if timestep is None:
            if name in self._xGuess:
                for k in range(self.nk+1):
                    self.guess(name,val,timestep=k,quiet=quiet)
                return
            if name in self._zGuess:
                for k in range(self.nk):
                    self.guess(name,val,timestep=k,quiet=quiet)
                return
            if name in self._uGuess:
                for k in range(self.nk):
                    self.guess(name,val,timestep=k,quiet=quiet)
                return
            if name in self._pGuess:
                oldguess = self._pGuess[name]
                if (self._pGuess[name] is not None):
                    raise ValueError("can't change guess on \""+name+"\" once it's set (tried to change guess "+str(oldguess)+" to "+str(val))
                self._pGuess[name] = val
                return
        
        assert(isinstance(timestep,int))
        if name in self._pGuess:
            raise ValueError("can't set guess on a parameter on a specific timestep")
        if name in self._xGuess:
            if timestep>self.nk+1:
                raise ValueError("differential state \""+name+"\" guess timestep too high ("+str(timestep)+">"+str(self.nk+1))
            oldguess = self._xGuess[name][timestep]
            if (oldguess is not None) and (quiet is False):
                print "WARNING: changing \"%s\" guess at timestep %d from %s to %s" % (name,timestep,str(oldguess),str(val))
            self._xGuess[name][timestep] = val
            return
        if name in self._zGuess:
            if timestep>self.nk:
                raise ValueError("algebraic state \""+name+"\" guess timestep too high ("+str(timestep)+">"+str(self.nk))
            oldguess = self._zGuess[name][timestep]
            if (oldguess is not None) and (quiet is False):
                print "WARNING: changing \"%s\" guess at timestep %d from %s to %s" % (name,timestep,str(oldguess),str(val))
            self._zGuess[name][timestep] = val
            return
        if name in self._uGuess:
            if timestep>self.nk:
                raise ValueError("action \""+name+"\" guess timestep too high ("+str(timestep)+">"+str(self.nk))
            oldguess = self._uGuess[name][timestep]
            if (oldguess is not None) and (quiet is False):
                print "WARNING: changing \"%s\" guess at timestep %d from %s to %s" % (name,timestep,str(oldguess),str(val))
            self._uGuess[name][timestep] = val
            return

    def _guessVec(self,val,names,**kwargs):
        if isinstance(val,list):
            length = len(val)
        elif hasattr(val,'size'):
            length = val.size()
        else:
            raise ValueError("can't figure out how long "+str(val)+" is")
        assert(len(names)==length)

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
        

    def lookup(self,name,timestep=None):
        assert(isinstance(name,str))
        
        if name in self.dae.pNames():
            if timestep is not None:
                raise ValueError("can't lookup a parameter ("+name+")at a certain timestep")
            k = self.dae.pNames().index(name)
            return self._P[k]

        if timestep is None:
            return CS.veccat([self.lookup(name,k) for k in range(self.nk+1)])

        assert(isinstance(timestep,int))
        if name in self.dae.xNames():
            k = self.dae.xNames().index(name)
            return self._XD[timestep][0][0][k]
        if name in self.dae.uNames():
            k = self.dae.uNames().index(name)
            return self._U[timestep][k]
        if name in self.dae.zNames():
            raise ValueError("can't lookup algebraic states at this time, bug Greg")

    def setObjective(self,obj):
        if hasattr(self,'_objective'):
            raise ValueError("can't change objective function once it's already set")
        self._objective = obj

    def plot(self,names,opt,title=None):
        if isinstance(names,str):
            names = [names]
        assert( isinstance(names,list) )
        plt.figure()
        plt.clf()
        legend = []
        for name in names:
            assert(isinstance(name,str))
            legend.append(name)

            # if it's a simple state or action
            if name in self.dae.xNames()+self.dae.uNames():
                plt.plot(opt['tgrid'],opt['vardict'][name])

            # if it's a Dae output
            elif name in self.dae._outputDict:
                y = []
                f = CS.SXFunction( [ self.dae.xVec(),self.dae.uVec(),self.dae.pVec()], [self.dae.output(name)] )
                f.init()
                for k,t in enumerate(opt['tgrid']):
                    f.setInput(opt['x'][:,k],0)
                    f.setInput(opt['u'][:,k],1)
                    f.setInput(opt['p'],2)
                    f.evaluate()
                    y.append(CS.DMatrix(f.output(0)))
                plt.plot(opt['tgrid'],y)
            else:
                raise ValueError("unrecognized name: \""+name+"\"")
        if isinstance(title,str):
            plt.title(title)
        plt.xlabel('time')
        plt.legend(legend)
        plt.grid()

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

    # Plot the results
    coll.plot(['x0','x1'],opt)
    coll.plot(['u','u^2'],opt)
    plt.show()
