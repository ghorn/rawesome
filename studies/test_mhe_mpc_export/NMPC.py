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
    
def makeNmpc(dae,N,dt):
    mpc = rawe.Ocp(dae, N=N, ts=dt)
    
    ocpOpts = rawe.OcpExportOptions()
    ocpOpts['HESSIAN_APPROXIMATION'] = 'GAUSS_NEWTON'
    ocpOpts['DISCRETIZATION_TYPE'] = 'MULTIPLE_SHOOTING'
    ocpOpts['QP_SOLVER'] = 'QP_QPOASES'
    ocpOpts['HOTSTART_QP'] = False
    ocpOpts['SPARSE_QP_SOLUTION'] = 'CONDENSING'
#   ocpOpts['SPARSE_QP_SOLUTION'] = 'FULL_CONDENSING_U2'
#   ocpOpts['AX_NUM_QP_ITERATIONS'] = '30'
    ocpOpts['FIX_INITIAL_STATE'] = True
               
    mpc.minimizeLsq(C.veccat([mpc['x'],mpc['v'],mpc['u']]))
    mpc.minimizeLsqEndTerm(C.veccat([mpc['x'],mpc['v']]))

    cgOpts = {'CXX':'g++', 'CC':'gcc'}
    return rawe.OcpRT(mpc, ocpOptions=ocpOpts, integratorOptions=intOpts,
                       codegenOptions=cgOpts)
