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

import rawe
import casadi as C

intOpts = rawe.RtIntegratorOptions()
intOpts['INTEGRATOR_TYPE'] = 'INT_IRK_GL2'
intOpts['NUM_INTEGRATOR_STEPS'] = 40
intOpts['IMPLICIT_INTEGRATOR_NUM_ITS'] = 3
intOpts['IMPLICIT_INTEGRATOR_NUM_ITS_INIT'] = 0
intOpts['LINEAR_ALGEBRA_SOLVER'] = 'HOUSEHOLDER_QR'
intOpts['UNROLL_LINEAR_SOLVER'] = False
intOpts['IMPLICIT_INTEGRATOR_MODE'] = 'IFTR'

def makeMhe(dae,N,dt):
    from rawe.ocp import Ocp
    mhe = Ocp(dae, N=N, ts=dt)
    
    ocpOpts = rawe.OcpExportOptions()
    ocpOpts['HESSIAN_APPROXIMATION'] = 'GAUSS_NEWTON'
    ocpOpts['DISCRETIZATION_TYPE'] = 'MULTIPLE_SHOOTING'
    ocpOpts['QP_SOLVER'] = 'QP_QPOASES'
    ocpOpts['HOTSTART_QP'] = False
    ocpOpts['SPARSE_QP_SOLUTION'] = 'CONDENSING'
#   ocpOpts['SPARSE_QP_SOLUTION'] = 'FULL_CONDENSING_U2'
#   ocpOpts['AX_NUM_QP_ITERATIONS'] = '30'
    ocpOpts['FIX_INITIAL_STATE'] = False
               
#    mhe.minimizeLsq(C.veccat([mhe['x'],mhe['u']]))
#    mhe.minimizeLsqEndTerm(C.veccat([mhe['x']]))
    mhe.minimizeLsq(mhe['measurements'])
    mhe.minimizeLsqEndTerm(mhe['measurementsN'])

    cgOpts = {'CXX':'g++', 'CC':'gcc'}
    print "makeMhe calling rawe.OcpRT"
    return rawe.OcpRT(mhe, ocpOptions=ocpOpts, integratorOptions=intOpts,
                       codegenOptions=cgOpts)
