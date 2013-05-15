import casadi as C

import phase1
import qpoases
from ocprt import OcpRT
import ocg_interface
from ..dae.rtIntegrator import rtModelExport
from ..utils import codegen

def validateOptions(defaultOpts, userOpts, optName):
    assert isinstance(userOpts,dict), optName+" options must be a dictionary"
    # fill in missing default options
    for name in defaultOpts:
        if name not in userOpts:
            userOpts[name] = defaultOpts[name]

    # throw error on extra, unrecognized options
    for name in userOpts:
        if name not in defaultOpts:
            raise Exception(optName+' option "'+name+'" unrecognized, valid options: '+\
                                str(defaultOpts.keys()))

def writeObjective(ocp, out0, exportName):
    dae = ocp._dae

    # first make out not a function of xDot or z
    xdot = C.veccat([dae.ddt(name) for name in dae.xNames()])
    inputs0 = [xdot, dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec()]
    outputFun0 = C.SXFunction(inputs0, [out0])

    (xDotDict, zDict) = dae.solveForXDotAndZ()
    xDot = C.veccat([xDotDict[name] for name in dae.xNames()])
    z    = C.veccat([zDict[name] for name in dae.zNames()])

    # plug in xdot, z solution to outputs fun
    outputFun0.init()
    [out] = outputFun0.eval([xDot, dae.xVec(), z, dae.uVec(), dae.pVec()])

    # make new SXFunction that is only fcn of [x, u, p]
    assert len(dae.pNames()) == 0, "parameters not supported right now in ocp export, sorry"
    if exportName == 'lsqExtern':
        inputs = C.veccat([dae.xVec(), dae.uVec(), dae.pVec()])
        outs = C.veccat( [ out, C.jacobian(out,dae.xVec()).T, C.jacobian(out,dae.uVec()).T ] )
        outputFun = C.SXFunction([inputs], [C.densify(outs)])
        outputFun.init()
        assert len(outputFun.getFree()) == 0, 'the "impossible" happened >_<'
    elif exportName == 'lsqEndTermExtern':
        inputs = dae.xVec()
        outs = C.veccat( [ out, C.jacobian(out,dae.xVec()).T ] )
        outputFun = C.SXFunction([inputs], [C.densify(outs)])
        outputFun.init()
        assert len(outputFun.getFree()) == 0, 'lsqEndTermExtern cannot be a function of controls u, saw: '+str(outputFun.getFree())
    else:
        raise Exception('unrecognized name "'+exportName+'"')

    return codegen.writeCCode(outputFun,exportName)


def exportOcp(ocp, cgOptions, acadoOptions, phase1Options):
    defaultCgOptions = {'CXX':'g++', 'CC':'gcc'}
    defaultPhase1Options = {'CXX':'g++'}
    validateOptions(defaultCgOptions, cgOptions, "codegen")
    validateOptions(defaultPhase1Options, phase1Options, "phase 1")

    # validate acado options (sort of)
    acadoOpsMsg = "acadoOptions must be a list of (string,string) tuples"
    assert isinstance(acadoOptions,list), acadoOpsMsg
    try:
        for key,val in acadoOptions:
            assert type(key) is str, acadoOpsMsg
            assert type(val) is str, acadoOpsMsg
    except Exception:
        raise Exception(acadoOpsMsg)
    try:
        qpSolver = dict(acadoOptions)['QP_SOLVER']
    except Exception:
        raise Exception('you must specify the QP_SOLVER in acado options')
    supportedQps = ['QP_QPOASES','QP_QPDUNES']
    assert qpSolver in supportedQps, "qp solver must be one of " + str(supportedQps)

    # write the OCP exporter and run it, returning an exported OCP
    files = phase1.runPhase1(ocp, phase1Options, acadoOptions, qpSolver)

    # add model for rt integrator
    files['model.c'] = '''\
#include "qpoases/solver.hpp"
#include "rhs.h"
#include "rhsJacob.h"
'''
    rtModelGen = rtModelExport.generateCModel(ocp._dae, ocp._ts)
    files['rhs.cpp'] = '#include "rhs.h"\n'+rtModelGen['rhsFile'][0]
    files['rhsJacob.cpp'] = '#include "rhsJacob.h"\n'+rtModelGen['rhsJacobFile'][0]
    files['rhs.h'] = rtModelGen['rhsFile'][1]
    files['rhsJacob.h'] = rtModelGen['rhsJacobFile'][1]

    # add objective and jacobian
    externObj    = writeObjective(ocp, ocp._minLsq, 'lsqExtern')
    externObjEnd = writeObjective(ocp, ocp._minLsqEndTerm, 'lsqEndTermExtern')
    files['lsqExtern.cpp'] = '#include "lsqExtern.h"\n'+externObj[0]
    files['lsqExtern.h']   = externObj[1]
    files['lsqEndTermExtern.cpp'] = '#include "lsqEndTermExtern.h"\n'+externObjEnd[0]
    files['lsqEndTermExtern.h']   = externObjEnd[1]

    # #include objective/jacobian in acado_solver.c
    files['acado_solver.c'] = '#include "lsqExtern.h"\n#include "lsqEndTermExtern.h"\n'+\
                              files['acado_solver.c']

    # add python_interface.c
    files['python_interface.c'] = ocg_interface.ocg_interface

    if qpSolver == 'QP_QPOASES':
        ocpSoPath = qpoases.exportPhase2(cgOptions, files)
    else:
        raise Exception('the impossible happened, unsupported qp solver: "'+str(qpSolver)+'"')

    return (ocpSoPath, ocp._ts, ocp._dae)
