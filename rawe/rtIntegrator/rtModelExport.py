import casadi as C

from ..utils import codegen

def generateCModel(dae,timeScaling,measurements):
    xdot = C.veccat([dae.ddt(name) for name in dae.xNames()])
    inputs = C.veccat([dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), xdot])
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
    jf = C.veccat( [ C.jacobian(f,inputs).T ] )
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
        jo = C.veccat( [ C.jacobian(measurements,inputs).T ] )
        measurementsJacobFun = C.SXFunction( [inputs], [C.densify(jo)] )
        measurementsJacobFun.init()
        measurementsJacobString = codegen.writeCCode(measurementsJacobFun, 'measurementsJacob')
        ret['measurementsJacob'] = measurementsJacobFun
        ret['measurementsJacobFile'] = measurementsJacobString

    return ret
