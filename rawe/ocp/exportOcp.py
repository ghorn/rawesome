import phase1
import qpoases
from ocprt import OcpRT
import ocg_interface

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

    # add model for rien integrator
    files['model.c'] = '''\
#include "qpoases/solver.hpp"

%(rhsAndJacString)s
''' % ocp._dae.makeRienModel(ocp._ts)

    # add python_interface.c
    files['python_interface.c'] = ocg_interface.ocg_interface

    if qpSolver is 'QP_QPOASES':
        ocpSoPath = qpoases.exportPhase2(cgOptions, files)
    else:
        raise Exception('the impossible happened, unsupported qp solver: "'+str(qpSolver)+'"')

    return OcpRT(ocpSoPath)
