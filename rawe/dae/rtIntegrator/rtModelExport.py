import casadi as C

from ...utils import codegen

def generateCModel(dae,timeScaling):
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

    # outputs
    o = C.veccat( [dae[outname] for outname in dae.outputNames()] )
    outputs = C.SXFunction( [inputs], [C.densify(o)] )
    outputs.init()
    outputsString = codegen.writeCCode(outputs, 'outputs')

    # outputs jacobian
    jo = C.veccat( [ C.jacobian(o,inputs).T ] )
    outputsJacob = C.SXFunction( [inputs], [C.densify(jo)] )
    outputsJacob.init()
    outputsJacobString = codegen.writeCCode(outputsJacob, 'outputsJacob')

    return {'rhs':rhs,
            'rhsJacob':rhsJacob,
            'oututs':outputs,
            'outputsJacob':outputs,
            'rhsFile':rhsString,
            'rhsJacobFile':rhsJacobString,
            'oututsFile':outputsString,
            'outputsJacobFile':outputsJacob}
