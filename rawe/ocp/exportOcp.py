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

import phase1
import qpoases
import ocg_interface
from ..rtIntegrator import rtModelExport
from ..utils import codegen

def validateOptions(defaultOpts, userOpts, optName):
    '''
    Fill in any missing options which have defaults.
    Throw an error if any given option is unrecognized (has no default).
    '''
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

# make sure each element in the output is only a function of x or u, not both
def testSeparation(dae,out,exportName):
    badOnes = {}
    msgs = []
    for k in range(out.size()):
        fx = C.SXFunction([dae.xVec(),dae.pVec()],[out[k]])
        fu = C.SXFunction([dae.uVec(),dae.pVec()],[out[k]])
        fx.init()
        fu.init()
        us = fx.getFree()
        xs = fu.getFree()
        if len(us) > 0 and len(xs) > 0:
            msgs.append('output '+str(k)+', xs: '+str(xs)+', us: '+str(us))
    if len(msgs) > 0:
        msg = str(len(msgs))+ ' compenents of '+exportName+' are functions '+\
              'of both x and u:\n'+'\n'.join(msgs)
        raise Exception(msg)

def writeObjective(ocp, out0, exportName):
    dae = ocp.dae

    # first make out not a function of xDot or z
    inputs0 = [dae.xDotVec(), dae.xVec(), dae.zVec(), dae.uVec(), dae.pVec()]
    outputFun0 = C.SXFunction(inputs0, [out0])

    (xDotDict, zDict) = dae.solveForXDotAndZ()
    xDot = C.veccat([xDotDict[name] for name in dae.xNames()])
    z    = C.veccat([zDict[name] for name in dae.zNames()])

    # plug in xdot, z solution to outputs fun
    outputFun0.init()
    [out] = outputFun0.eval([xDot, dae.xVec(), z, dae.uVec(), dae.pVec()])

    assert len(dae.pNames()) == 0, "parameters not supported right now in ocp export, sorry"

    # make sure each element in the output is only a function of x or u, not both
    testSeparation(dae,out,exportName)

    # make new SXFunction that is only fcn of [x, u, p]
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


def exportOcp(ocp, ocpOptions, integratorOptions, cgOptions, phase1Options):
    defaultCgOptions = {'CXX':'g++', 'CC':'gcc',
                        'CXXFLAGS':'-O3 -fPIC -finline-functions',
                        'CFLAGS':'-O3 -fPIC -finline-functions',
                        'hideSymbols':False,
                        'export_without_build_path':None}
    defaultPhase1Options = {'CXX':'g++'}
    validateOptions(defaultCgOptions, cgOptions, "codegen")
    validateOptions(defaultPhase1Options, phase1Options, "phase 1")
    cgOptions['hashPrefix'] = ocp.hashPrefix
    phase1Options['hashPrefix'] = ocp.hashPrefix

    # write the OCP exporter and run it, returning an exported OCP
    files = phase1.runPhase1(ocp, phase1Options, integratorOptions, ocpOptions)

    # add model for rt integrator
    files['model.c'] = '''\
#include "qpoases/solver.hpp"
#include "rhs.h"
#include "rhsJacob.h"
'''
    rtModelGen = rtModelExport.generateCModel(ocp.dae, ocp.ts, None)
    files['rhs.cpp'] = '#include "rhs.h"\n'+rtModelGen['rhsFile'][0]
    files['rhsJacob.cpp'] = '#include "rhsJacob.h"\n'+rtModelGen['rhsJacobFile'][0]
    files['rhs.h'] = rtModelGen['rhsFile'][1]
    files['rhsJacob.h'] = rtModelGen['rhsJacobFile'][1]

    # add objective and jacobian
    externObj    = writeObjective(ocp, ocp._minLsq, 'lsqExtern')
    externObjEnd = writeObjective(ocp, ocp._minLsqEndTerm, 'lsqEndTermExtern')
    externFile  = '''\
#include "acado_external_functions.h"
#include <math.h>
#include <stdio.h>
'''
    externFile += externObj[0] + '\n'
    externFile += externObjEnd[0]
    files['acado_external_functions.cpp'] = externFile
    externHeader  = externObj[1] + '\n'
    externHeader += externObjEnd[1]
    files['acado_external_functions.h'] = externHeader

    # #include objective/jacobian in acado_solver.c
    files['acado_solver.c'] = '#include "acado_external_functions.h"\n'+\
                              files['acado_solver.c']

    # add python_interface.c
    files['python_interface.c'] = ocg_interface.ocg_interface

    if ocpOptions['QP_SOLVER'] == 'QP_QPOASES':
        exportPath = qpoases.exportPhase2(cgOptions, files)
    else:
        raise Exception('the impossible happened, unsupported qp solver: "'+str(ocpOptions['QP_SOLVER'])+'"')

    return exportPath
