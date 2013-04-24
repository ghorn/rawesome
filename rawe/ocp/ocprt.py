import ctypes
import numpy
import casadi as C

class OcpRT(object):
    _canonicalNames = ['x','u','y','yN','x0','S','SN']
    def __init__(self,libpath):
        print 'loading "'+libpath+'"'
        self._lib = ctypes.cdll.LoadLibrary(libpath)

        # set return types of KKT and objective
        self._lib.getKKT.restype = ctypes.c_double
        self._lib.getObjective.restype = ctypes.c_double


        print 'initializing solver'
        self._lib.py_initialize()
        self._libpath = libpath

        self.x  = numpy.zeros((self._lib.py_get_ACADO_N()+1,
                               self._lib.py_get_ACADO_NX()))
        self.u  = numpy.zeros((self._lib.py_get_ACADO_N(),
                               self._lib.py_get_ACADO_NU()))
        self.y  = numpy.zeros((self._lib.py_get_ACADO_N(),
                               self._lib.py_get_ACADO_NY()))
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
        self.getAll()

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

    def setAll(self):
        self._callMat(self._lib.py_set_x,  self.x)
        self._callMat(self._lib.py_set_u,  self.u)
        self._callMat(self._lib.py_set_y,  self.y)
        self._callMat(self._lib.py_set_yN, self.yN)
        self._callMat(self._lib.py_set_S,  self.S)
        self._callMat(self._lib.py_set_SN, self.SN)
        if self._lib.py_get_ACADO_INITIAL_STATE_FIXED():
            self._callMat(self._lib.py_set_x0, self.x0)

    def getAll(self):
        self._callMat(self._lib.py_get_x,  self.x)
        self._callMat(self._lib.py_get_u,  self.u)
        self._callMat(self._lib.py_get_y,  self.y)
        self._callMat(self._lib.py_get_yN, self.yN)
        self._callMat(self._lib.py_get_S,  self.S)
        self._callMat(self._lib.py_get_SN, self.SN)
        if self._lib.py_get_ACADO_INITIAL_STATE_FIXED():
            self._callMat(self._lib.py_get_x0, self.x0)

    def preparationStep(self):
        self.setAll()
        ret = self._lib.preparationStep()
        self.getAll()
        return ret

    def feedbackStep(self):
        self.setAll()
        ret = self._lib.feedbackStep()
        self.getAll()
        return ret

    def initializeNodesByForwardSimulation(self):
        self.setAll()
        self._lib.initializeNodesByForwardSimulation()
        self.getAll()
                
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
            self.y = numpy.append(self.y[1:,:],new_y,axis=0)
            
        if new_yN != None:
            self.yN = new_yN
           
        
        
#     def shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd ):
#         void shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd );

#     def shiftControls( real_t* const uEnd ):
#         void shiftControls( real_t* const uEnd );

    def getKKT(self):
        self.setAll()
        return self._lib.getKKT()

    def getObjective(self):
        self.setAll()
        return self._lib.getObjective()
