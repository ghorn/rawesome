import ctypes
import os
import numpy

import casadi as C

import rienModelExport
import rienIntegratorInterface

from ..utils import codegen, subprocess_tee

def loadIntegratorInterface():
    # write the interface file
    files = {'rienIntegratorInterface.cpp':rienIntegratorInterface.phase1src(),
             'Makefile':rienIntegratorInterface.phase1makefile()}
    interfaceDir = codegen.memoizeFiles(files)

    # call make to make sure shared lib is build
    (ret, msgs) = subprocess_tee.call(['make',codegen.makeJobs()], cwd=interfaceDir)
    if ret != 0:
        raise Exception("integrator compilation failed:\n"+msgs)

    # load the shared object
    return ctypes.cdll.LoadLibrary(os.path.join(interfaceDir, 'rienIntegratorInterface.so'))

def makeMakefile(cfiles, cxxfiles):
    return """\
CC      = gcc
CFLAGS  = -O3 -fPIC -finline-functions -I.
CXX     = g++
CXXFLAGS = -O3 -fPIC -finline-functions -I.
LDFLAGS = -lm

C_SRC = %(cfiles)s
CXX_SRC = %(cxxfiles)s
OBJ = $(C_SRC:%%.c=%%.o)
OBJ += $(CXX_SRC:%%.cpp=%%.o)

.PHONY: clean all
all : $(OBJ) model.so integrator.so

%%.o : %%.c acado.h
\t@echo CC $@: $(CC) $(CFLAGS) -c $< -o $@
\t@$(CC) $(CFLAGS) -c $< -o $@

%%.o : %%.cpp acado.h
\t@echo CXX $@: $(CXX) $(CXXFLAGS) -c $< -o $@
\t@$(CXX) $(CXXFLAGS) -c $< -o $@

%%.so : $(OBJ)
\t@echo LD $@: $(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)
\t@$(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)

clean :
\trm -f *.o *.so
""" % {'cfiles':' '.join(cfiles), 'cxxfiles':' '.join(cxxfiles)}


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
    rienModelGen = rienModelExport.generateCModel(dae,options['timestep'])
    modelFile = '''\
#include "acado.h"
#include "rhs.h"
#include "rhsJacob.h"
'''
    # write the makefile
    makefile = makeMakefile(['workspace.c', 'model.c', 'integrator.c'],
                            ['rhs.cpp','rhsJacob.cpp'])

    # write the static workspace file (temporary)
    workspace = """\
#include <acado.h>
ACADOworkspace acadoWorkspace;
ACADOvariables acadoVariables;
"""
    genfiles = {'integrator.c': exportedFiles['integrator.c'],
                'acado.h': exportedFiles['acado.h'],
                'model.c': modelFile,
                'rhs.cpp': '#include "rhs.h"\n'+rienModelGen['rhsFile'][0],
                'rhs.h': rienModelGen['rhsFile'][1],
                'rhsJacob.cpp': '#include "rhsJacob.h"\n'+rienModelGen['rhsJacobFile'][0],
                'rhsJacob.h': rienModelGen['rhsJacobFile'][1],
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
    _canonicalNames = ['x','z',
                       'dx1_dx0','dz0_dx0','dx1_du','dx1_dp','dz0_du','dz0_dp',
                       'u','p',
                       '_dx1z0_dx0','_dx1z0_dup',
                       '_data']
    def __setattr__(self, name, value):
        if name in self._canonicalNames:
            if type(value)==C.DMatrix:
                value = numpy.array(value)
            if type(value)==dict:
                if name == 'x':
                    value = numpy.array([value[n] for n in self._dae.xNames()])
                elif name == 'z':
                    value = numpy.array([value[n] for n in self._dae.zNames()])
                elif name == 'u':
                    value = numpy.array([value[n] for n in self._dae.uNames()])
                elif name == 'p':
                    value = numpy.array([value[n] for n in self._dae.pNames()])
                else:
                    raise Exception('you can only pass a dict for [x,z,u,p], not for '+name)

            if hasattr(self, name):
                assert value.shape == getattr(self, name).shape, \
                    name+' has dimension '+str(getattr(self,name).shape)+' but you tried to '+\
                    'assign it something with dimension '+str(value.shape)
            object.__setattr__(self, name, numpy.ascontiguousarray(value, dtype=numpy.double))
        else:
            object.__setattr__(self, name, value)

    def _setData(self):
        self._dx1z0_dx0 = numpy.vstack( (self.dx1_dx0,
                                         self.dz0_dx0) )
        self._dx1z0_dup = numpy.vstack((numpy.hstack( (self.dx1_du, self.dx1_dp) ),
                                        numpy.hstack( (self.dz0_du, self.dz0_dp) )))
        self._data = numpy.concatenate((self.x,
                                        self.z,
                                        self._dx1z0_dx0.flatten(),
                                        self._dx1z0_dup.flatten(),
                                        self.u,
                                        self.p))

    def _getData(self):
        i0 = 0
        i1 = 0
        for field in ['x','z','_dx1z0_dx0','_dx1z0_dup','u','p']:
            i0  = i1
            i1 += getattr(self,field).size
            shape = getattr(self,field).shape
            self.__setattr__(field, self._data[i0:i1].reshape(shape))
        assert i1 == self._data.size
        # unpack dx1z0_dx0, dx1z0_dup
        nx = self.x.size
        nu = self.u.size
        self.dx1_dx0 = self._dx1z0_dx0[:nx,:]
        self.dz0_dx0 = self._dx1z0_dx0[nx:,:]
        self.dx1_du  = self._dx1z0_dup[:nx,:nu]
        self.dx1_dp  = self._dx1z0_dup[:nx,nu:]
        self.dz0_du  = self._dx1z0_dup[nx:,:nu]
        self.dz0_dp  = self._dx1z0_dup[nx:,nu:]

    def __init__(self, dae, ts, numIntegratorSteps=10, integratorType='INT_IRK_GL4'):
        self._dae = dae
        self._ts = ts

        # set some options
        options = {}
        options['timestep'] = ts # because we scale xdot
        options['numIntegratorSteps'] = numIntegratorSteps
        options['integratorType'] = integratorType
        options['integratorGrid'] = None

        # setup outputs function
        self._outputsFun = self._dae.outputsFunWithSolve()

        (integratorLib, modelLib, rienModelGen) = exportIntegrator(self._dae, options)
        self._integratorLib = integratorLib
        self._modelLib = modelLib
        self._rienModelGen = rienModelGen
        
        self._initIntegrator = 1

        nx = len( self._dae.xNames() )
        nz = len( self._dae.zNames() )
        nu = len( self._dae.uNames() )
        np = len( self._dae.pNames() )

#        [ x z d(x,z)/dx d(x,z)/d(u,p) u p]
        self.x = numpy.zeros( nx )
        self.z = numpy.zeros( nz )
        self.u = numpy.zeros( nu )
        self.p = numpy.zeros( np )
        self.dx1_dx0  = numpy.zeros( (nx, nx) )
        self.dz0_dx0  = numpy.zeros( (nz, nx) )
        self.dx1_du = numpy.zeros( (nx, nu) )
        self.dx1_dp = numpy.zeros( (nx, np) )
        self.dz0_du = numpy.zeros( (nz, nu) )
        self.dz0_dp = numpy.zeros( (nz, np) )

        self._dx1z0_dx0 = numpy.zeros( (nx+nz, nx) )
        self._dx1z0_dup = numpy.zeros( (nx+nz, nu+np) )

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

    def step(self,x=None,u=None,p=None):
        # x,u,p can be dicts or array-like
        # if x is a dict, the return value is a dict, otherwise it's a numpy array

        # vectorize inputs
        if x != None:
            self.x = x
        if u != None:
            self.u = u
        if p != None:
            self.p = p

        # call integrator
        self._setData()
        ret = self._integratorLib.integrate(ctypes.c_void_p(self._data.ctypes.data),
                                            self._initIntegrator)
        self._getData()
        self._initIntegrator = 0

        # devectorize outputs
        if x != None and type(x) == dict:
            xret = {}
            for k,name in enumerate(self._dae.xNames()):
                xret[name] = self.x[k]
        else:
            xret = numpy.copy(self.x)
        return xret

    def getOutputs(self, x=None, u=None, p=None):
        # vectorize inputs
        if x != None:
            self.x = x
        if u != None:
            self.u = u
        if p != None:
            self.p = p
        self._outputsFun.setInput(self.x, 0)
        self._outputsFun.setInput(self.u, 1)
        self._outputsFun.setInput(self.p, 2)
        self._outputsFun.evaluate()
        ret = {}
        for k,name in enumerate(self._dae.outputNames()):
            ret[name] = numpy.array(self._outputsFun.output(k))
        return ret
