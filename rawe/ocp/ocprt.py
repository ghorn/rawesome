import ctypes
import numpy
import copy
import matplotlib.pyplot as plt
import casadi as C
import scipy

def dlqr(A, B, Q, R, N):
    
    if N == None:
        N = numpy.zeros((Q.shape[0],R.shape[0]))
    
    P = scipy.linalg.solve_discrete_are(A, B, Q, R)
    
    k1 = numpy.dot(numpy.dot(B.T, P), B) + R
    k2 = numpy.dot(numpy.dot(B.T, P), A)
    K = numpy.linalg.solve(k1,k2)
    
    return K, P

#class Logger(object):
#    def __init__(self,ocprt,dae):
#
#        if ocprt._lib.py_get_ACADO_INITIAL_STATE_FIXED():
#            self._canonicalNames = ocprt._canonicalNames
#        else:
#            self._canonicalNames = list(set(ocprt._canonicalNames)-set(['x0']))
#        
#        self.xNames = dae.xNames()
#        self.uNames = dae.uNames()
#        self.Ts = ocprt.getTs()
#
#        self._ocprt = ocprt
#        self._dae = dae
#        self._log = {}
#        for field in self._canonicalNames:
#            self._log[field] = []
#        self._log['kkt'] = []
#        self._log['objective'] = []
#        self._log['timing'] = []
#        self.log(ocprt)
#        
#    def log(self, ocprt):
#        for field in self._canonicalNames:
#            self._log[field].append(copy.deepcopy(getattr(ocprt, field)))
#        self._log['kkt'].append(ocprt.getKKT())
#        self._log['objective'].append(ocprt.getObjective())
#
#    def subplot(self,names,title=None,style='',when=0,showLegend=True):
#        assert isinstance(names,list)
#
#        fig = plt.figure()
#        if title is None:
#            if isinstance(names,str):
#                title = names
#            else:
#                assert isinstance(names,list)
#                if len(names) == 1:
#                    title = names[0]
#                else:
#                    title = str(names)
#        fig.canvas.set_window_title(str(title))
#
#        plt.clf()
#        n = len(names)
#        if style is '':
#            style = ['']*n
#        for k,name in enumerate(names):
#            plt.subplot(n,1,k+1)
#            if k==0:
#                self._plot(name,title,style[k],when=when,showLegend=showLegend)
#            else:
#                self._plot(name,None,style[k],when=when,showLegend=showLegend)
#
#    def plot(self,names,title=None,style='',when=0,showLegend=True):
#
#        fig = plt.figure()
#        if title is None:
#            if isinstance(names,str):
#                title = names
#            else:
#                assert isinstance(names,list)
#                if len(names) == 1:
#                    title = names[0]
#                else:
#                    title = str(names)
#        fig.canvas.set_window_title(str(title))
#
#        plt.clf()
#        self._plot(names,title,style,when=when,showLegend=showLegend)
#
#
#    def _plot(self,names,title,style,when=0,showLegend=True):
#        if isinstance(names,str):
#            names = [names]
#        assert isinstance(names,list)
#        
##        if style == None:
##            style = ''
#        
#        legend = []
#        for name in names:
#            assert isinstance(name,str)
#            legend.append(name)
#
#            # if it's a differential state
#            if name in self.xNames:
#                index = self.xNames.index(name)
#                if when == 'all':
#                    for k in range(numpy.array(self._log['x']).shape[0]-1):
#                        ys = numpy.array(self._log['x'])[1+k,:,index]
#                        ts = numpy.arange(len(ys))*self.Ts + self.Ts*k
#                        plt.plot(ts,ys,style)
#                else:
#                    ys = numpy.array(self._log['x'])[1:,when,index]
#                    ts = numpy.arange(len(ys))*self.Ts
#                    plt.plot(ts,ys,style)
#
#            # if it's a control
#            if name in self.uNames:
#                index = self.uNames.index(name)
#                ys = numpy.array(self._log['u'])[1:,when,index]
#                ts = numpy.arange(len(ys))*self.Ts
#                plt.step(ts,ys,style)
#
#        if title is not None:
#            assert isinstance(title,str), "title must be a string"
#            plt.title(title)
#        plt.xlabel('time [s]')
#        if showLegend is True:
#            plt.legend(legend)
#        plt.grid()




class OcpRT(object):
    _canonicalNames = ['x','u','z','y','yN','x0','S','SN']
    def __init__(self,libpath, ts, dae):
        self._dae = dae
        self._ts = ts
        print 'loading "'+libpath+'"'
        self._lib = ctypes.cdll.LoadLibrary(libpath)

        # set return types of KKT,objective,etc
        self._lib.getKKT.restype = ctypes.c_double
        self._lib.getObjective.restype = ctypes.c_double
        self._lib.preparationStepTimed.restype = ctypes.c_double
        self._lib.feedbackStepTimed.restype = ctypes.c_double

        self.preparationTime = 0.0
        self.feedbackTime = 0.0

        print 'initializing solver'
        self._lib.py_initialize()
        self._libpath = libpath

        self.x  = numpy.zeros((self._lib.py_get_ACADO_N()+1,
                               self._lib.py_get_ACADO_NX()))
        self.u  = numpy.zeros((self._lib.py_get_ACADO_N(),
                               self._lib.py_get_ACADO_NU()))
        self.y  = numpy.zeros((self._lib.py_get_ACADO_N(),
                               self._lib.py_get_ACADO_NY()))
        if self._lib.py_get_ACADO_NXA() > 0:
            self.z  = numpy.zeros((self._lib.py_get_ACADO_N(),
                                   self._lib.py_get_ACADO_NXA()))

        self.yN = numpy.zeros(self._lib.py_get_ACADO_NYN())
        wmt = self._lib.py_get_ACADO_WEIGHTING_MATRICES_TYPE()
        if wmt == 1:
            self.S  = numpy.zeros((self._lib.py_get_ACADO_NY(),
                                   self._lib.py_get_ACADO_NY()))
            self.SN = numpy.zeros((self._lib.py_get_ACADO_NYN(),
                                   self._lib.py_get_ACADO_NYN()))
        elif wmt == 2:
            self.S  = numpy.zeros((self._lib.py_get_ACADO_N()*self._lib.py_get_ACADO_NY(),
                                   self._lib.py_get_ACADO_NY()))
            self.SN = numpy.zeros((self._lib.py_get_ACADO_NYN(),
                                   self._lib.py_get_ACADO_NYN()))
        else:
            raise Exception('unrecognized ACADO_WEIGHING_MATRICES_TYPE '+str(wmt))

        if self._lib.py_get_ACADO_INITIAL_STATE_FIXED():
            self.x0 = numpy.zeros(self._lib.py_get_ACADO_NX())

        self._lib.py_initialize()
        self._getAll()
        
        if not self._lib.py_get_ACADO_INITIAL_STATE_FIXED():
            self._canonicalNames = list(set(ocprt._canonicalNames)-set(['x0']))
        
        self.xNames = self.dae.xNames()
        self.uNames = self.dae.uNames()
        self.Ts = self.getTs()

        self._log = {}
        for field in self._canonicalNames:
            self._log[field] = []
        self._log['kkt'] = []
        self._log['objective'] = []
        self._log['timing'] = []
        self.log()

    def __setattr__(self, name, value):
        if name in self._canonicalNames:
            if type(value)==C.DMatrix:
                value = numpy.array(value)
            if hasattr(self, name):
                assert value.shape == getattr(self, name).shape, \
                    name+' has dimension '+str(getattr(self,name).shape)+' but you tried to '+\
                    'assign it something with dimension '+str(value.shape)
            object.__setattr__(self, name, numpy.ascontiguousarray(value, dtype=numpy.double))
        else:
            object.__setattr__(self, name, value)

    def _callMat(self,call,mat):
        sh = mat.shape
        # if it's 1 dimensional with size n, treat is as shape (n,1)
        if len(sh) == 1:
            sh = (sh[0], 1)
        (nr,nc) = sh
        ret = call(ctypes.c_void_p(mat.ctypes.data), nr, nc)
        assert 0 == ret, "dimension mismatch in "+str(call)
        return call(ctypes.c_void_p(mat.ctypes.data), nr, nc)

    def _setAll(self):
        self._callMat(self._lib.py_set_x,  self.x)
        self._callMat(self._lib.py_set_u,  self.u)
        self._callMat(self._lib.py_set_y,  self.y)
        self._callMat(self._lib.py_set_yN, self.yN)
        self._callMat(self._lib.py_set_S,  self.S)
        self._callMat(self._lib.py_set_SN, self.SN)
        if self._lib.py_get_ACADO_INITIAL_STATE_FIXED():
            self._callMat(self._lib.py_set_x0, self.x0)
        if hasattr(self, 'z'):
            self._callMat(self._lib.py_set_z, self.z)

    def _getAll(self):
        self._callMat(self._lib.py_get_x,  self.x)
        self._callMat(self._lib.py_get_u,  self.u)
        self._callMat(self._lib.py_get_y,  self.y)
        self._callMat(self._lib.py_get_yN, self.yN)
        self._callMat(self._lib.py_get_S,  self.S)
        self._callMat(self._lib.py_get_SN, self.SN)
        if self._lib.py_get_ACADO_INITIAL_STATE_FIXED():
            self._callMat(self._lib.py_get_x0, self.x0)
        if hasattr(self, 'z'):
            self._callMat(self._lib.py_get_z, self.z)

    def preparationStep(self):
        self._setAll()
        self.preparationTime = self._lib.preparationStepTimed()
        self._getAll()

    def feedbackStep(self):
        self._setAll()
        ret = ctypes.c_int(0)
        self.feedbackTime = self._lib.feedbackStepTimed(ctypes.byref(ret))
        self._getAll()
        if ret.value != 0:
            raise Exception("feedbackStep returned error code "+str(ret.value))


    def initializeNodesByForwardSimulation(self):
        self._setAll()
        self._lib.initializeNodesByForwardSimulation()
        self._getAll()

    def shift(self,new_x=None,new_u=None,sim=None,new_y=None,new_yN=None,new_S=None,new_SN=None):
        # Shift weighting matrices
        if new_S != None:
            wmt = self._lib.py_get_ACADO_WEIGHTING_MATRICES_TYPE()
            if wmt == 1:    # Constant weighting matrices
                self.S =  new_S

            elif wmt == 2:  # Varying weighting matrices
                weight_S = self.S[1:(self._lib.py_get_ACADO_N())*self._lib.py_get_ACADO_NY(),:]
                self.S = numpy.ascontiguousarray(numpy.append(weight_S,new_S,axis=0), dtype=numpy.double)

            else:
                raise Exception('unrecognized ACADO_WEIGHING_MATRICES_TYPE '+str(wmt))

        if new_SN != None:
            self.SN = new_SN

        # Shift states and controls
        guess_x = self.x[1:,:]
        guess_u = self.u[1:,:]

        if new_x != None and new_u == None:
            raise Exception('if you provide new_x you must also provide new_u')
        if new_x != None and sim == None:
            raise Exception('you cannot provide both new_x and sim')
        if new_u != None and new_x != None and sim != None:
            raise Exception('your sim will never be used')
        if new_u != None:   # If the new control is provided, either integrate to get the new state or use the provided one
            if new_x != None: # If the new state is provided, use it
                self.x = numpy.append(guess_x,new_x,axis=0)
                self.u = numpy.append(guess_u,new_u,axis=0)
            elif sim != None: # Integrate the system forward with the provided integrator
                # Integrate the system forward using the last control
                new_x = sim.step(self.x[-1,:],new_u,{})

                self.u = numpy.append(guess_u,new_u,axis=0)
                self.x = numpy.append(guess_x,new_x.T,axis=0)
            else:
                raise Exception('If a new control is provided as an input to the shift function either a state or an integrator must also be provided')

        elif sim != None: # If an integrator is provided, use it
            if new_x != None:
                raise Exception('your new_x will never be used')
            new_u = self.u[-1,:]
            # Integrate the system forward using the last control
            new_x = sim.step(self.x[-1,:],new_u,{})

            self.u = numpy.append(guess_u,[new_u],axis=0)
            self.x = numpy.append(guess_x,new_x.T,axis=0)

        else:
            self.x = numpy.append(guess_x,[self.x[-1,:]],axis=0)
            self.u = numpy.append(guess_u,[self.u[-1,:]],axis=0)

        # Shift the reference (if provided, else keep the previous one)
        if new_y != None:
            self.y = numpy.append(self.y[1:,:],[new_y],axis=0)

        if new_yN != None:
            self.yN = new_yN
            
    def log(self):
        for field in self._canonicalNames:
            self._log[field].append(copy.deepcopy(getattr(self, field)))
        self._log['kkt'].append(self.getKKT())
        self._log['objective'].append(self.getObjective())


#     def shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd ):
#         void shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd );

#     def shiftControls( real_t* const uEnd ):
#         void shiftControls( real_t* const uEnd );

#    def getTiming(self):
#        blabla
        
    def getKKT(self):
        self._setAll()
        return self._lib.getKKT()

    def getObjective(self):
        self._setAll()
        return self._lib.getObjective()

    def getTs(self):
        return self._ts
        
    def subplot(self,names,title=None,style='',when=0,showLegend=True):
        assert isinstance(names,list)

        fig = plt.figure()
        if title is None:
            if isinstance(names,str):
                title = names
            else:
                assert isinstance(names,list)
                if len(names) == 1:
                    title = names[0]
                else:
                    title = str(names)
        fig.canvas.set_window_title(str(title))

        plt.clf()
        n = len(names)
        if style is '':
            style = ['']*n
        for k,name in enumerate(names):
            plt.subplot(n,1,k+1)
            if k==0:
                self._plot(name,title,style[k],when=when,showLegend=showLegend)
            else:
                self._plot(name,None,style[k],when=when,showLegend=showLegend)

    def plot(self,names,title=None,style='',when=0,showLegend=True):

        fig = plt.figure()
        if title is None:
            if isinstance(names,str):
                title = names
            else:
                assert isinstance(names,list)
                if len(names) == 1:
                    title = names[0]
                else:
                    title = str(names)
        fig.canvas.set_window_title(str(title))

        plt.clf()
        self._plot(names,title,style,when=when,showLegend=showLegend)


    def _plot(self,names,title,style,when=0,showLegend=True):
        if isinstance(names,str):
            names = [names]
        assert isinstance(names,list)
        
#        if style == None:
#            style = ''
        
        legend = []
        for name in names:
            assert isinstance(name,str)
            legend.append(name)

            # if it's a differential state
            if name in self.xNames:
                index = self.xNames.index(name)
                if when == 'all':
                    for k in range(numpy.array(self._log['x']).shape[0]-1):
                        ys = numpy.array(self._log['x'])[1+k,:,index]
                        ts = numpy.arange(len(ys))*self.Ts + self.Ts*k
                        plt.plot(ts,ys,style)
                else:
                    ys = numpy.array(self._log['x'])[1:,when,index]
                    ts = numpy.arange(len(ys))*self.Ts
                    plt.plot(ts,ys,style)

            # if it's a control
            if name in self.uNames:
                index = self.uNames.index(name)
                ys = numpy.array(self._log['u'])[1:,when,index]
                ts = numpy.arange(len(ys))*self.Ts
                plt.step(ts,ys,style)

        if title is not None:
            assert isinstance(title,str), "title must be a string"
            plt.title(title)
        plt.xlabel('time [s]')
        if showLegend is True:
            plt.legend(legend)
        plt.grid()


class MpcRT(OcpRT):
    def __init__(self, dae, lqrDae, libpath, ts, referenceStr, reference):
        OcpRT.__init__(self, dae, libpath, ts)
        self.lqrDae = lqrDae
        self.referenceStr = referenceStr
        self.reference = reference

    def computeLqr(self):
        nx = self.x.shape[1]
#        nu = self.u.shape[1]
        
        self.integrator.x = self.y[-1,:nx]
        self.integrator.u = self.y[-1,nx:]
        self.integrator.step()
        A = self.integrator.dx1_dx0
        B = self.integrator.dx1_du
    #    integrator.getOutputs()
        
        K, P = dlqr(A, B, self.Q, self.R, self.N)
        
        self.K = K
        self.SN = P

class MheRT(OcpRT):
    def UpdateArrivalCost(self): 
        ''' Arrival cost implementation.
            Approximate the solution of:
            min_{xL_,uL_,xL1_} ||  pL ( xL_-xL )         ||^2
                               ||  vL ( yL-h(xL_,uL_) )  ||
                               ||  wL wx                 ||_2
                         s.t.  wx = xL1_ - f(xL_,uL_)
            where:
                    PL = pL^T pL is the last kalman prediction covariance matrix
                    VL = vL^T vL is the measurement noise covariance matrix
                    WL = wL^T wL is the state noise covariance matrix
            
            Linearization (at the last MHE estimate x,u which is different from xL,uL):
            f(xL_,uL_) ~= f(x,u) + df(x,u)/dx (xL_-x) + df(x,u)/du (uL_-u)
                       ~= f(x,u) +         Xx (xL_-x) +         Xu (uL_-u)
                       ~= f(x,u) - Xx x - Xu u + Xx xL_ + Xu uL_
                       ~= x_tilde              + Xx xL_ + Xu uL_
            h(xL_,uL_) ~= h(x,u) + dh(x,u)/dx (xL_-x) + dh(x,u)/du (uL_-u)
                       ~= f(x,u) +         Hx (xL_-x) +         Hu (uL_-u)
                       ~= h(x,u) - Hx x - Hu u + Hx xL_ + Hu uL_
                       ~= h_tilde              + Hx xL_ + Hu uL_
                       
            Linearized problem:
            min_{xL_,uL_,xL1_} ||  pL ( xL_ - xL )                          ||^2
                               ||  vL ( yL - h_tilde - Hx xL_ - Hu uL_ )    ||
                               ||  wL ( xL1_ - x_tilde - Xx xL_ - Xu uL_ )  ||_2
            
            Rewrite as:
            min_{xL_,uL_,xL1_} ||  M ( xL_, uL_, xL1_ ) + res  ||^2_2
            
            After QR factorization of M:
            min_{xL_,uL_,xL1_} ||  R ( xL_, uL_, xL1_ ) + rho  ||^2_2
                
            '''        
        pL = self.pL
        vL = self.vL
        wL = self.wL
        
        xL = self.xL        # Last kalman update state prediction for the initial state
        yL = self.y[0,:]    # Initial measurement
        
        x = self.x[0,:]     # Last MHE state prediction
        u = self.u[0,:]     # Last MHE control prediction
        
        nx = x.shape[0]
        nu = u.shape[0]
        nV = vL.shape[0]
        
        self.integrator.x = x
        self.integrator.u = u
        h = self.integrator.y
        x1 = self.integrator.step()
        Xx = self.integrator.dx1_dx0
        Xu = self.integrator.dx1_du
        
        Hx = self.integrator.dy_dx0
        Hu = self.integrator.dy_du
        
        x_tilde = x1 - numpy.dot(Xx,x) - numpy.dot(Xu,u)
        h_tilde =  h - numpy.dot(Hx,x) - numpy.dot(Hu,u)
        
        res = numpy.bmat([ -numpy.dot(pL, xL),
                            numpy.dot(vL, yL - h_tilde),
                           -numpy.dot(wL, x_tilde) ])
        res = numpy.squeeze(numpy.array(res))
        
        M = numpy.bmat([[                pL,  numpy.zeros((nx,nu)), numpy.zeros((nx,nx)) ],
                        [ -numpy.dot(vL,Hx),     -numpy.dot(vL,Hu), numpy.zeros((nV,nx)) ],
                        [ -numpy.dot(wL,Xx),     -numpy.dot(wL,Xu),                   wL ]])
        
        Q, R = numpy.linalg.qr(M)
        
    #    R1  = R[:nx+nu,:nx+nu]
    #    R12 = R[:nx+nu,nx+nu:]
        R2  = R[nx+nu:,nx+nu:]
        
    #    rho = numpy.linalg.solve(Q,res)
        rho = numpy.squeeze(numpy.array(numpy.dot(Q.T,res)))
        rho2 = rho[nx+nu:]
        
        pL1 = R2
        xL1 = -numpy.linalg.solve(R2,rho2)
        
        self.pL = numpy.array( pL1 )
        self.AC = numpy.dot( pL1.T, pL1 )
        
        self.xL = numpy.array( xL1 )
