import ctypes
import os
import numpy
import rienIntegratorInterface

from ..utils import codegen, subprocess_tee

def loadIntegratorInterface():
    # write the interface file
    files = {'rienIntegratorInterface.cpp':rienIntegratorInterface.phase1src,
             'Makefile':rienIntegratorInterface.phase1makefile}
    interfaceDir = codegen.memoizeFiles(files)

    # call make to make sure shared lib is build
    (ret, msgs) = subprocess_tee.call(['make',codegen.makeJobs()], cwd=interfaceDir)
    if ret != 0:
        raise Exception("integrator compilation failed:\n"+msgs)

    # load the shared object
    return ctypes.cdll.LoadLibrary(os.path.join(interfaceDir, 'rienIntegratorInterface.so'))

def makeMakefile(cfiles):
    return """\
CC      = gcc
CFLAGS  = -O3 -fPIC -finline-functions -I.
LDFLAGS = -lm

C_SRC = %(cfiles)s
OBJ = $(C_SRC:%%.c=%%.o)

.PHONY: clean all
all : $(OBJ) model.so integrator.so

%%.o : %%.c acado.h
\t@echo CC $@: $(CC) $(CFLAGS) -c $< -o $@
\t@$(CC) $(CFLAGS) -c $< -o $@

%%.so : $(OBJ)
\t@echo LD $@: $(CC) -shared -o $@ $(OBJ) $(LDFLAGS)
\t@$(CC) -shared -o $@ $(OBJ) $(LDFLAGS)

clean :
\trm -f *.o *.so
""" % {'cfiles':' '.join(cfiles)}


def writeRienIntegrator(dae, options):
    nx = len(dae.xNames())
    nz = len(dae.zNames())
    nup = len(dae.uNames()) + len(dae.pNames())

    # call makeRienIntegrator
    lib = loadIntegratorInterface()
    numIntervals = 1 # should maybe be hard-coded
    def call(path):
        ret = lib.makeRienIntegrator(ctypes.c_char_p(path),
                                     numIntervals,
                                     ctypes.c_double(1.0),
                                     ctypes.c_char_p(options['integratorType']),
                                     options['integratorGrid'],
                                     options['numIntegratorSteps'],
                                     nx, nz, nup)
        if ret != 0:
            raise Exception("Rien integrator creater failed")
    return codegen.withTempdir(call)

def exportIntegrator(dae, options):
    # get the exported integrator files
    exportedFiles = writeRienIntegrator(dae, options)

    # model file
    rienModelGen = dae.makeRienModel(options['timestep'])
    rienModelGen['modelFile'] = '#include "acado.h"\n\n'+rienModelGen['modelFile']
    modelFile = rienModelGen['modelFile'] 

    # write the makefile
    makefile = makeMakefile(['workspace.c', 'model.c', 'integrator.c'])

    # write the static workspace file (temporary)
    workspace = """\
#include <acado.h>
ACADOworkspace acadoWorkspace;
ACADOvariables acadoVariables;
"""

    genfiles = {'integrator.c': exportedFiles['integrator.c'],
                'acado.h': exportedFiles['acado.h'],
                'model.c': modelFile,
                'workspace.c': workspace,
                'Makefile': makefile}
    exportpath = codegen.memoizeFiles(genfiles)

    # compile the code
    (ret, msgs) = subprocess_tee.call(['make',codegen.makeJobs()], cwd=exportpath)
    if ret != 0:
        raise Exception("integrator compilation failed:\n"+msgs)

    print 'loading '+exportpath+'/integrator.so'
    integratorLib = ctypes.cdll.LoadLibrary(exportpath+'/integrator.so')
    print 'loading '+exportpath+'/model.so'
    modelLib = ctypes.cdll.LoadLibrary(exportpath+'/model.so')
    return (integratorLib, modelLib, rienModelGen)


class RienIntegrator(object):
    def __init__(self, dae, ts, numIntegratorSteps=10, integratorType='INT_IRK_GL4'):
        self._dae = dae

        # set some options
        options = {}
        options['timestep'] = ts # because we scale xdot
        options['numIntegratorSteps'] = numIntegratorSteps
        options['integratorType'] = integratorType
        options['integratorGrid'] = None

        (integratorLib, modelLib, rienModelGen) = exportIntegrator(self._dae, options)
        self._integratorLib = integratorLib
        self._modelLib = modelLib
        self._rienModelGen = rienModelGen
        
        self._initIntegrator = 1

        nx = len( self._dae.xNames() )
        nz = len( self._dae.zNames() )
        nu = len( self._dae.uNames() )
        np = len( self._dae.pNames() )
        N = nx + nz + nu + np + (nx+nz)*nx + (nu+np)*nx
        
        self._data = numpy.zeros( N, dtype=numpy.double )
        self._xvec = self._data[0:nx]
        self._zvec = self._data[nx:nx+nz]
        self._pvec = self._data[(N-np):N]
        self._uvec = self._data[(N-np-nu):(N-np)]

        assert self._xvec.size == nx, str(self._xvec.size)+'!='+str(nx)
        assert self._zvec.size == nz, str(self._zvec.size)+'!='+str(nz)
        assert self._pvec.size == np, str(self._pvec.size)+'!='+str(np)
        assert self._uvec.size == nu, str(self._uvec.size)+'!='+str(nu)

    def rhs(self,xdot,x,z,u,p, compareWithSX=False):
        xdot = numpy.array([xdot[n] for n in self._dae.xNames()],dtype=numpy.double)
        x    = numpy.array([x[n]    for n in self._dae.xNames()],dtype=numpy.double)
        z    = numpy.array([z[n]    for n in self._dae.zNames()],dtype=numpy.double)
        u    = numpy.array([u[n]    for n in self._dae.uNames()],dtype=numpy.double)
        p    = numpy.array([p[n]    for n in self._dae.pNames()],dtype=numpy.double)
        dataIn = numpy.concatenate((x,z,u,p,xdot))
        dataOut = numpy.zeros(x.size + z.size, dtype=numpy.double)
        
        self._modelLib.rhs(ctypes.c_void_p(dataIn.ctypes.data),
                           ctypes.c_void_p(dataOut.ctypes.data),
                           )

        if compareWithSX:
            f = self._rienModelGen['rhs']
            f.setInput(dataIn)
            f.evaluate()
            print f.output() - dataOut

        return dataOut

    def rhsJac(self,xdot,x,z,u,p, compareWithSX=False):
        xdot = numpy.array([xdot[n] for n in self._dae.xNames()],dtype=numpy.double)
        x    = numpy.array([x[n]    for n in self._dae.xNames()],dtype=numpy.double)
        z    = numpy.array([z[n]    for n in self._dae.zNames()],dtype=numpy.double)
        u    = numpy.array([u[n]    for n in self._dae.uNames()],dtype=numpy.double)
        p    = numpy.array([p[n]    for n in self._dae.pNames()],dtype=numpy.double)
        dataIn = numpy.concatenate((x,z,u,p,xdot))
        dataOut = numpy.zeros((x.size + z.size)*(2*x.size+z.size+u.size+p.size), dtype=numpy.double)
        
        self._modelLib.rhs_jac(ctypes.c_void_p(dataIn.ctypes.data),
                               ctypes.c_void_p(dataOut.ctypes.data),
                               )
        if compareWithSX:
            f = self._rienModelGen['rhsJacob']
            f.setInput(dataIn)
            f.evaluate()
            print (f.output() - dataOut)
        return dataOut

    def run(self,*args,**kwargs):
        raise Exception("to step a rien integrator, you now have to call .step(x,u,p) instead of .run(x,u,p)")

    def step(self,x,u,p):
        # x,u,p can be dicts or array-like
        # if x is a dict, the return value is a dict, otherwise it's a numpy array

        # vectorize inputs
        if type(x) == dict:
            for k,name in enumerate(self._dae.xNames()):
                self._xvec[k] = x[name]
        else:
            for k in range(len(self._dae.xNames())):
                self._xvec[k] = x[k]

        if type(u) == dict:
            for k,name in enumerate(self._dae.uNames()):
                self._uvec[k] = u[name]
        else:
            for k in range(len(self._dae.uNames())):
                self._uvec[k] = u[k]

        if type(p) == dict:
            for k,name in enumerate(self._dae.pNames()):
                self._pvec[k] = p[name]
        else:
            for k in range(len(self._dae.pNames())):
                self._pvec[k] = p[k]

        # call integrator
        ret = self._integratorLib.integrate(ctypes.c_void_p(self._data.ctypes.data), self._initIntegrator)
        self._initIntegrator = 0

        # devectorize outputs
        if type(x) == dict:
            xret = {}
            for k,name in enumerate(self._dae.xNames()):
                xret[name] = self._xvec[k]
        else:
            xret = numpy.copy(self._xvec)
        return xret
