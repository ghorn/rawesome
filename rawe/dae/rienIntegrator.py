import ctypes
import os
import subprocess
import numpy

from ..utils.codegen import memoizeFiles, withTempdir

def loadIntegratorInterface():
    # get the integrator interface directory
    filename = __file__.rstrip('.pyc')
    filename = filename.rstrip('.py')
    filename = filename.rstrip('rienIntegrator')
    interfaceDir = os.path.join(filename , 'rienIntegratorInterface')

    # call make to make sure shared lib is build
    p = subprocess.Popen(['make'], stdout=subprocess.PIPE, cwd=interfaceDir)
    ret = p.wait()

    if ret != 0:
        print p.stdout.read()
        raise Exception("integrator compilation failed, return code "+str(ret))

    # load the shared object
    return ctypes.cdll.LoadLibrary(os.path.join(interfaceDir, 'rienIntegratorInterface.so'))

def makeMakefile(cfiles):
    return """\
CC      = gcc
CFLAGS  = -O3 -fPIC -finline-functions -I.
LDFLAGS = -lm

C_SRC = %(cfiles)s
OBJ = $(C_SRC:%%.c=%%.o)

.PHONY: clean libs obj

all : model.so integrator.so

%%.o : %%.c acado.h
\t@echo CC $@: $(CC) $(CFLAGS) -c $< -o $@
\t@$(CC) $(CFLAGS) -c $< -o $@

%%.so : $(OBJ)
\t@echo LD $@: $(CC) -shared -Wl,-soname,$@ -o $@ $(OBJ) $(LDFLAGS)
\t@$(CC)   -shared -Wl,-soname,$@ -o $@ $(OBJ) $(LDFLAGS)

# aliases
obj : $(OBJ)
clean :
\trm -f *.o *.so
""" % {'cfiles':' '.join(cfiles)}


def writeRienIntegrator(dae, path, options):
    nx = len(dae.xNames())
    nz = len(dae.zNames())
    nup = len(dae.uNames()) + len(dae.pNames())

    # call makeRienIntegrator
    lib = loadIntegratorInterface()
    ret = lib.makeRienIntegrator(ctypes.c_char_p(path),
                                 options['numIntervals'],
                                 ctypes.c_double(1.0),
                                 ctypes.c_char_p(options['integratorType']),
                                 options['integratorGrid'],
                                 options['numIntegratorSteps'],
                                 nx, nz, nup)
    if ret != 0:
        raise Exception("rien integrator creater, what goon set bad options?")


def exportIntegrator(dae, options):
    # get the exported integrator files
    exportedFiles = withTempdir(lambda tmppath : writeRienIntegrator(dae, tmppath, options))

    # model file
    rienModelGen = dae.makeRienModel(options['timestep'])
    modelFile = rienModelGen['modelFile']

    # write the makefile
    makefile = makeMakefile(['workspace.c', 'model.c', 'integrator.c'])

    # write the static workspace file (temporary)
    workspace = """\
#include <acado.h>
ACADOworkspace acadoWorkspace;
ACADOvariables acadoVariables;
"""

    genfiles = [('integrator.c', exportedFiles['integrator.c']),
                ('acado.h', exportedFiles['acado.h']),
                ('model.c', modelFile),
                ('workspace.c', workspace),
                ('Makefile', makefile)]
    exportpath = memoizeFiles(genfiles)

    # compile the code
    p = subprocess.Popen(['make'], cwd=exportpath)
    ret = p.wait()
    if ret != 0:
        raise Exception("integrator compilation failed, return code "+str(ret))

    #p = subprocess.Popen(['make'], stdout=subprocess.PIPE, cwd=exportpath)
    #p.wait()
    #ret = p.stdout.read()
    #print "ret: " +str(ret)

    print 'loading '+exportpath+'/integrator.so'
    integratorLib = ctypes.cdll.LoadLibrary(exportpath+'/integrator.so')
    print 'loading '+exportpath+'/model.so'
    modelLib = ctypes.cdll.LoadLibrary(exportpath+'/model.so')
    return (integratorLib, modelLib, rienModelGen)


class RienIntegrator(object):
    def __init__(self, dae, ts, numIntervals=1, numIntegratorSteps=10, integratorType='INT_IRK_RIIA3'):
        self._dae = dae

        # set some options
        options = {}
        options['numIntervals'] = numIntervals
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

    def run(self,x,u,p):
        # vectorize inputs
        for k,name in enumerate(self._dae.xNames()):
            self._xvec[k] = x[name]
        for k,name in enumerate(self._dae.uNames()):
            self._uvec[k] = u[name]
        for k,name in enumerate(self._dae.pNames()):
            self._pvec[k] = p[name]

        # call integrator
        ret = self._integratorLib.integrate(ctypes.c_void_p(self._data.ctypes.data), self._initIntegrator)
        self._initIntegrator = 0

        # devectorize outputs
        xret = {}
        for k,name in enumerate(self._dae.xNames()):
            xret[name] = self._xvec[k]
        return xret
