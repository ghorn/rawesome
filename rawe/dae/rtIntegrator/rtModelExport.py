import casadi as C

from ...utils import codegen

def generateCModel(dae,timeScaling,outputs):
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

    if outputs is not None:
        # outputs
        outputsFun = C.SXFunction( [inputs], [outputs] )
        outputsFun.init()
        [outputs] = outputsFun.eval([C.veccat([dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec(), xdot/timeScaling])])
        outputsFun = C.SXFunction( [inputs], [C.densify(outputs)] )
        outputsFun.init()
        outputsString = codegen.writeCCode(outputsFun, 'outputs')
        ret['outputs'] = outputsFun
        ret['outputsFile'] = outputsString

        # outputs jacobian
        jo = C.veccat( [ C.jacobian(outputs,inputs).T ] )
        outputsJacobFun = C.SXFunction( [inputs], [C.densify(jo)] )
        outputsJacobFun.init()
        outputsJacobString = codegen.writeCCode(outputsJacobFun, 'outputsJacob')
        ret['outputsJacob'] = outputsJacobFun
        ret['outputsJacobFile'] = outputsJacobString

    return ret
