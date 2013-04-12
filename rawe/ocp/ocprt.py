import ctypes
import numpy

class OcpRT(object):
    def __init__(self,libpath):
        print 'loading "'+libpath+'"'
        self._lib = ctypes.cdll.LoadLibrary(libpath)

        # set return types of KKT and objective
        self._lib.getKKT.restype = ctypes.c_double
        self._lib.getObjective.restype = ctypes.c_double


        print 'initializing solver'
        self._lib.py_initialize()
        self._libpath = libpath

        self.x  = self._myZeros(self._lib.py_get_ACADO_N()+1,
                                self._lib.py_get_ACADO_NX())
        self.u  = self._myZeros(self._lib.py_get_ACADO_N(),
                                self._lib.py_get_ACADO_NU())
        self.y  = self._myZeros(self._lib.py_get_ACADO_N(),
                                self._lib.py_get_ACADO_NY())
        self.yN = self._myZeros(self._lib.py_get_ACADO_NYN(), 1)
        wmt = self._lib.py_get_ACADO_WEIGHTING_MATRICES_TYPE()
        if wmt == 1:
            self.S  = self._myZeros(self._lib.py_get_ACADO_NY(),
                                    self._lib.py_get_ACADO_NY())
            self.SN = self._myZeros(self._lib.py_get_ACADO_NYN(),
                                    self._lib.py_get_ACADO_NYN())
        elif wmt == 2:
            self.S  = self._myZeros(self._lib.py_get_ACADO_N()*self._lib.py_get_ACADO_NY(),
                                    self._lib.py_get_ACADO_NY())
            self.SN = self._myZeros(self._lib.py_get_ACADO_NYN(),
                                    self._lib.py_get_ACADO_NYN())
        else:
            raise Exception('unrecognized ACADO_WEIGHING_MATRICES_TYPE '+str(wmt))

        if self._lib.py_get_ACADO_INITIAL_STATE_FIXED():
            self.x0 = self._myZeros(self._lib.py_get_ACADO_NX(), 1)

        self._lib.py_initialize()
        self.getAll()

    def _myZeros(self, nr, nc):
        z = numpy.zeros( (nr, nc), dtype=numpy.double )
        return numpy.ascontiguousarray(z, dtype=numpy.double)

    def _callMat(self,call,mat):
        (nr,nc) = mat.shape
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
