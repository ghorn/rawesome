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

import ctypes
import numpy

import casadi as C

from rtIntegratorExport import exportIntegrator

from ..utils import codegen, subprocess_tee
from ..utils.options import Options, OptStr, OptInt, OptBool

class RtIntegratorOptions(Options):
    def __init__(self):
        Options.__init__(self, 'RtIntegrator')
        self.add(OptStr('LINEAR_ALGEBRA_SOLVER',['GAUSS_LU','HOUSEHOLDER_QR'],default='GAUSS_LU'))
        self.add(OptInt('NUM_INTEGRATOR_STEPS',default=1))
        integratorTypes = \
            ['INT_EX_EULER',
             'INT_RK2','INT_RK3','INT_RK4',
             'INT_IRK_GL2','INT_IRK_GL4','INT_IRK_GL6','INT_IRK_GL8',
             'INT_IRK_RIIA1','INT_IRK_RIIA3','INT_IRK_RIIA5',
             'INT_DIRK3','INT_DIRK4','INT_DIRK5',
             'INT_DT','INT_NARX']
        self.add(OptStr('INTEGRATOR_TYPE',integratorTypes,default='INT_IRK_GL4'))
        self.add(OptStr('IMPLICIT_INTEGRATOR_MODE',['IFTR','IFT'],default='IFTR'))
        self.add(OptInt('IMPLICIT_INTEGRATOR_NUM_ITS',default=3))
        self.add(OptInt('IMPLICIT_INTEGRATOR_NUM_ITS_INIT',default=0))
        self.add(OptBool('UNROLL_LINEAR_SOLVER',default=False))
#        self.add(OptStr('MEASUREMENT_GRID',['EQUIDISTANT_SUBGRID','EQUIDISTANT_GRID','ONLINE_GRID']))

class RtIntegrator(object):
    _canonicalNames = ['x','z',
                       'dx1_dx0','dz0_dx0','dx1_du','dx1_dp','dz0_du','dz0_dp',
                       'u','p',
                       '_dx1z0_dx0','_dx1z0_dup',
                       'h','dh_dx0','dh_du','dh_dp','_dh_dup', '_measData',# measurements
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
        if self._measurements is not None:
            self._dh_dup = numpy.hstack( (self.dh_du, self.dh_dp) )
            self._measData = numpy.concatenate((self.h,
                                                self.dh_dx0.flatten(),
                                                self._dh_dup.flatten()))

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

        if self._measurements is not None:
            i0 = 0
            i1 = 0
            for field in ['h','dh_dx0','_dh_dup']:
                i0  = i1
                i1 += getattr(self,field).size
                shape = getattr(self,field).shape
                self.__setattr__(field, self._measData[i0:i1].reshape(shape))
            assert i1 == self._measData.size
            # unpack dh_dup
            self.dh_du = self._dh_dup[:,:nu]
            self.dh_dp = self._dh_dup[:,nu:]

    def __init__(self, dae, ts, measurements=None, options=RtIntegratorOptions()):
        self._dae = dae
        self._ts = ts
        if measurements is None:
            self._measurements = measurements
        else:
            if isinstance(measurements,list):
                measurements = C.veccat(measurements)
            self._measurements = measurements

        # setup outputs function
        self._outputsFun = self._dae.outputsFunWithSolve()

        (integratorLib, modelLib, rtModelGen) = exportIntegrator(self._dae, ts, options, self._measurements)
        self._integratorLib = integratorLib
        self._modelLib = modelLib
        self._rtModelGen = rtModelGen
        
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

        if self._measurements is not None:
            nh = self._measurements.size()
            self.h = numpy.zeros( nh )
            self.dh_dx0 = numpy.zeros( (nh, nx) )
            self.dh_du = numpy.zeros( (nh, nu) )
            self.dh_dp = numpy.zeros( (nh, np) )

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
            f = self._rtModelGen['rhs']
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
            f = self._rtModelGen['rhsJacob']
            f.setInput(dataIn)
            f.evaluate()
            print (f.output() - dataOut)
        return dataOut

    def run(self,*args,**kwargs):
        raise Exception("to step an rt integrator, you now have to call .step(x,u,p) instead of .run(x,u,p)")

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
        if self._measurements is None:
            ret = self._integratorLib.integrate(ctypes.c_void_p(self._data.ctypes.data),
                                                self._initIntegrator)
        else:
            ret = self._integratorLib.integrate(ctypes.c_void_p(self._data.ctypes.data),
                                                ctypes.c_void_p(self._measData.ctypes.data),
                                                self._initIntegrator)
        assert ret==0, "integrator returned error: "+str(ret)
        self._getData()
        self._initIntegrator = 0
#        print "h:"
#        print self.h
#        print "dh_dx0:"
#        print self.dh_dx0
#        print "dh_du:"
#        print self.dh_du
#        print "dh_dp:"
#        print self.dh_dp
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
