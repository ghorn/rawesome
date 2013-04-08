import subprocess
import ctypes
import os
import numpy

import phase1
import ocg_interface
from rawe.utils import codegen,pkgconfig

def exportOcp(ocp, cgOptions, acadoOptions, qpSolver):
    assert isinstance(cgOptions, dict), "codegen options must be a dictionary"
    if 'CXX' not in cgOptions:
        cgOptions['CXX'] = 'g++'
    if 'CC' not in cgOptions:
        cgOptions['CC'] = 'gcc'

    # validate acado options (sort of)
    acadoOpsMsg = "acadoOptions must be a list of (string,string) tuples"
    assert isinstance(acadoOptions,list), acadoOpsMsg
    try:
        for key,val in acadoOptions:
            assert type(key) is str, acadoOpsMsg
            assert type(val) is str, acadoOpsMsg
    except Exception:
        raise Exception(acadoOpsMsg)

    # write the OCP exporter and run it, returning an exported OCP
    files = phase1.runPhase1(ocp, cgOptions, acadoOptions, qpSolver)

    # add model for rien integrator
    files['model.c'] = '''\
#include "qpoases/solver.hpp"

%(rhsAndJacString)s
''' % ocp._dae.makeRienModel(ocp._ts)

    # add python_interface.c
    files['python_interface.c'] = ocg_interface.ocg_interface

    if qpSolver is 'QP_OASES':
        ocpret = exportQpOases(cgOptions, files)
    else:
        raise Exception('the impossible happened, unsupported qp solver: "'+str(qpSolver)+'"')

    return ocpret


def qpoasesMakefile(cgOptions, qposrc):
    qposrc = ' \\\n'.join(['\t'+os.path.join('qpoases', q.split('qpoases'+os.sep)[1]) for q in qposrc])
    makefile = """\
CXX      = %(CXX)s
CC       = %(CC)s
CXXFLAGS = -O3 -fPIC -finline-functions
CFLAGS   = -O3 -fPIC -finline-functions

#CFLAGS   += -Wall -Wextra
#CXXFLAGS += -Wall -Wextra

LDFLAGS = -lm -lrt

CXX_SRC = \\
%(qpo_src)s \\
\tqpoases/solver.cpp

C_SRC = \\
\tworkspace.c \\
\tpython_interface.c \\
\tmodel.c \\
\tacado_integrator.c \\
\tacado_solver.c \\
\tacado_auxiliary_functions.c

QPO_INC = \\
\t-I. \\
\t-Iqpoases \\
\t-Iqpoases/INCLUDE \\
\t-Iqpoases/SRC

CXXFLAGS += $(QPO_INC)

CXX_OBJ = $(CXX_SRC:%%.cpp=%%.o)
C_OBJ = $(C_SRC:%%.c=%%.o)

HEADERS = \\
\tqpoases/solver.hpp \\
\tacado_auxiliary_functions.h \\
\tacado_common.h

.PHONY: clean all ocp.a ocp.so
all : $(CXX_OBJ) $(C_OBJ) ocp.a ocp.so

$(CXX_OBJ) : %%.o : %%.cpp $(HEADERS)
\t@echo CXX $@: $(CXX) $(CXXFLAGS) -c $< -o $@
\t@$(CXX) $(CXXFLAGS) -c $< -o $@

model.o :: CFLAGS += -Wall -Wextra -Werror -Wno-unused-variable
workspace.o :: CFLAGS += -Wall -Wextra -Werror
python_interface.o :: CFLAGS += -Wall -Wextra -Werror
$(C_OBJ) : %%.o : %%.c $(HEADERS)
\t@echo CC $@: $(CC) $(CFLAGS) -c $< -o $@
\t@$(CC) $(CFLAGS) -c $< -o $@

ocp.so : $(CXX_OBJ) $(C_OBJ)
\t@echo LD $@: $(CXX) -shared -o $@ $(CXX_OBJ) $(C_OBJ) $(LDFLAGS)
\t@$(CXX) -shared -o $@ $(CXX_OBJ) $(C_OBJ) $(LDFLAGS)

ocp.a : $(CXX_OBJ) $(C_OBJ)
\t@echo AR $@ : ar r $@ $?
\t@ar r $@ $?

clean :
\t@echo rm -f ocp.so ocp.a $(CXX_OBJ) $(C_OBJ) #*.o *.a ./qpoases/SRC/*.o ./qpoases/SRC/*.a test
\t@rm -f ocp.so ocp.a $(CXX_OBJ) $(C_OBJ)
""" % {'CXX':cgOptions['CXX'], 'CC':cgOptions['CC'], 'qpo_src':qposrc}
    return makefile

def exportQpOases(cgOptions, phase1src):
    # call pkg-config to get qpoases source and includes
    qpoStuff = {}
    for name in ['qpOASESsrc', 'qpOASESinc']:
        qpoStuff[name] = pkgconfig.call(['--variable',name,'acado'])
    qpoStuff['qpOASESsrc'] = qpoStuff['qpOASESsrc'].split(' ')

    # get qpoases source as file dictionary
    qpoSrcPath = os.path.join(qpoStuff['qpOASESsrc'][0].split('qpoases')[0], 'qpoases')
    phase2src = codegen.directoryToDict(qpoSrcPath)

    # merge qpoases source with phase 1 source
    def mergeAll(srcdict,destdict):
        for name,src in srcdict.items():
            if isinstance(src,dict):
                if name not in destdict:
                    destdict[name] = {}
                assert isinstance(destdict[name], dict), "dictionary merge failed, source was directory but destination was a file"
                destdict[name] = mergeAll(src,destdict[name])
            else:
                destdict[name] = src
        return destdict
    genfiles = mergeAll(phase1src, {'qpoases':phase2src})

    # add makefile
    genfiles['Makefile'] = qpoasesMakefile(cgOptions, qpoStuff['qpOASESsrc'])
    genfiles['workspace.c'] ='''\
#include "acado_common.h"
ACADOworkspace acadoWorkspace;
ACADOvariables acadoVariables;
'''

    # write all this
    exportpath = codegen.memoizeFiles(genfiles)
    print exportpath

    # compile!
    p = subprocess.Popen(['make',codegen.makeJobs()], stderr=subprocess.PIPE, cwd=exportpath)
    if p.wait() != 0:
        raise Exception("ocp compilation failed:\n"+p.stderr.read())

    # load the result
    libpath = os.path.join(exportpath, 'ocp.so')

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

            self.x  = numpy.zeros( (self._lib.py_get_ACADO_N()+1,
                                    self._lib.py_get_ACADO_NX()), dtype=numpy.double)
            self.u  = numpy.zeros( (self._lib.py_get_ACADO_N(),
                                    self._lib.py_get_ACADO_NU()), dtype=numpy.double)
            self.y  = numpy.zeros( (self._lib.py_get_ACADO_N(),
                                    self._lib.py_get_ACADO_NY()), dtype=numpy.double)
            self.yN = numpy.zeros( (self._lib.py_get_ACADO_NYN(), 1), dtype=numpy.double)
            wmt = self._lib.py_get_ACADO_WEIGHTING_MATRICES_TYPE()
            if wmt == 1:
                self.S  = numpy.zeros( (self._lib.py_get_ACADO_NY(),
                                        self._lib.py_get_ACADO_NY()),
                                       dtype=numpy.double)
                self.SN = numpy.zeros( (self._lib.py_get_ACADO_NYN(),
                                        self._lib.py_get_ACADO_NYN()),
                                       dtype=numpy.double)
            elif wmt == 2:
                self.S  = numpy.zeros( (self._lib.py_get_ACADO_N()*self._lib.py_get_ACADO_NY(),
                                        self._lib.py_get_ACADO_NY()),
                                       dtype=numpy.double)
                self.SN = numpy.zeros( (self._lib.py_get_ACADO_NYN(),
                                        self._lib.py_get_ACADO_NYN()),
                                       dtype=numpy.double)
            else:
                raise Exception('unrecognized ACADO_WEIGHING_MATRICES_TYPE '+str(wmt))

            if self._lib.py_get_ACADO_INITIAL_STATE_FIXED():
                self.x0 = numpy.zeros( (self._lib.py_get_ACADO_NX(), 1), dtype=numpy.double)

            self._lib.py_initialize()
            self.getAll()

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

#        def shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd ):
#            void shiftStates( int strategy, real_t* const xEnd, real_t* const uEnd );

#        def shiftControls( real_t* const uEnd ):
#            void shiftControls( real_t* const uEnd );

        def getKKT(self):
            return self._lib.getKKT()

        def getObjective(self):
            return self._lib.getObjective()

    return OcpRT(libpath)
