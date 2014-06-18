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

import casadi as C

from ..utils import codegen

def generateCModel(dae,timeScaling,measurements):
    xdot = C.veccat([dae.ddt(name) for name in dae.xNames()])
    inputs = C.veccat([dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), xdot])
    jacobian_inputs = C.veccat([dae.xVec(), dae.zVec(), dae.uVec(), xdot])
    f = dae.getResidual()

    # dae residual
    rhs = C.SXFunction( [inputs], [f] )
    rhs.init()
    # handle time scaling
    [f] = rhs.eval([C.veccat([dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), xdot/timeScaling])])
    rhs = C.SXFunction( [inputs], [C.densify(f)] )
    rhs.init()
    rhsString = codegen.writeCCode(rhs, 'rhs')

    # dae residual jacobian
    jf = C.veccat( [ C.jacobian(f,jacobian_inputs).T ] )
    rhsJacob = C.SXFunction( [inputs], [C.densify(jf)] )
    rhsJacob.init()
    rhsJacobString = codegen.writeCCode(rhsJacob, 'rhsJacob')

    ret = {'rhs':rhs,
           'rhsJacob':rhsJacob,
           'rhsFile':rhsString,
           'rhsJacobFile':rhsJacobString}

    if measurements is not None:
        # measurements
        measurementsFun = C.SXFunction( [inputs], [measurements] )
        measurementsFun.init()
        [measurements] = measurementsFun.eval([C.veccat([dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), xdot/timeScaling])])
        measurementsFun = C.SXFunction( [inputs], [C.densify(measurements)] )
        measurementsFun.init()
        measurementsString = codegen.writeCCode(measurementsFun, 'measurements')
        ret['measurements'] = measurementsFun
        ret['measurementsFile'] = measurementsString

        # measurements jacobian
        jo = C.veccat( [ C.jacobian(measurements,jacobian_inputs).T ] )
        measurementsJacobFun = C.SXFunction( [inputs], [C.densify(jo)] )
        measurementsJacobFun.init()
        measurementsJacobString = codegen.writeCCode(measurementsJacobFun, 'measurementsJacob')
        ret['measurementsJacob'] = measurementsJacobFun
        ret['measurementsJacobFile'] = measurementsJacobString

    return ret
