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
import os
from multiprocessing import Process, Queue

from ..utils import codegen, subprocess_tee
import rtModelExport
import rtIntegratorInterface

def makeMakefile(cfiles, cxxfiles):
    return """\
CC      = gcc
CFLAGS  = -O3 -fPIC -finline-functions -I.
CXX     = g++
CXXFLAGS = -O3 -fPIC -finline-functions -I.
LDFLAGS = -lm

C_SRC = %(cfiles)s
CXX_SRC = %(cxxfiles)s
OBJ = $(C_SRC:%%.c=%%.o)
OBJ += $(CXX_SRC:%%.cpp=%%.o)

.PHONY: clean all
all : $(OBJ) model.so integrator.so

%%.o : %%.c acado.h
\t@echo CC $@: $(CC) $(CFLAGS) -c $< -o $@
\t@$(CC) $(CFLAGS) -c $< -o $@

%%.o : %%.cpp acado.h
\t@echo CXX $@: $(CXX) $(CXXFLAGS) -c $< -o $@
\t@$(CXX) $(CXXFLAGS) -c $< -o $@

%%.so : $(OBJ)
\t@echo LD $@: $(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)
\t@$(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)

clean :
\trm -f *.o *.so
""" % {'cfiles':' '.join(cfiles), 'cxxfiles':' '.join(cxxfiles)}


def writeRtIntegrator(dae, options, measurements):
    # write the exporter file
    files = {'export_integrator.cpp':rtIntegratorInterface.phase1src(dae, options, measurements),
             'Makefile':rtIntegratorInterface.phase1makefile()}
    interfaceDir = codegen.memoizeFiles(files,prefix='rt_integrator_phase1__')

    # call make to make sure shared lib is build
    (ret, msgs) = subprocess_tee.call(['make',codegen.makeJobs()], cwd=interfaceDir)
    if ret != 0:
        raise Exception("integrator compilation failed:\n"+msgs)

    # call makeRtIntegrator
    def call(path):
        # load the shared object
        lib = ctypes.cdll.LoadLibrary(os.path.join(interfaceDir, 'export_integrator.so'))
        ret = lib.export_integrator(ctypes.c_char_p(path))
        if ret != 0:
            raise Exception("Rt integrator creater failed")
    def callInProcess(q):
        try:
            q.put(codegen.withTempdir(call))
        finally:
            q.put(None)

    q = Queue()
    p = Process(target=callInProcess,args=(q,))
    p.start()
    ret = q.get()
    p.join()
    assert (0 == p.exitcode) and (ret is not None), \
        "error exporting integrator, see stdout/stderr above"
    return ret

def exportIntegrator(dae, timestep, options, measurements):
    # get the exported integrator files
    exportedFiles = writeRtIntegrator(dae, options, measurements)

    # model file
    rtModelGen = rtModelExport.generateCModel(dae,timestep, measurements)
    modelFile = '''\
#include "acado.h"
#include "rhs.h"
#include "rhsJacob.h"
'''
    if measurements is not None:
        modelFile += '''\
#include "measurements.h"
#include "measurementsJacob.h"
'''

    # write the makefile
    symbolicsFiles = ['rhs.cpp','rhsJacob.cpp']
    if measurements is not None:
        symbolicsFiles += ['measurements.cpp', 'measurementsJacob.cpp']
    makefile = makeMakefile(['workspace.c', 'model.c', 'integrator.c'],
                            symbolicsFiles)

    # write the static workspace file (temporary)
    workspace = """\
#include <acado.h>
ACADOworkspace acadoWorkspace;
ACADOvariables acadoVariables;
"""
    genfiles = {'integrator.c': exportedFiles['integrator.c'],
                'acado.h': exportedFiles['acado.h'],
                'model.c': modelFile,
                'rhs.cpp': '#include "rhs.h"\n'+rtModelGen['rhsFile'][0],
                'rhs.h': rtModelGen['rhsFile'][1],
                'rhsJacob.cpp': '#include "rhsJacob.h"\n'+rtModelGen['rhsJacobFile'][0],
                'rhsJacob.h': rtModelGen['rhsJacobFile'][1],
                'workspace.c': workspace,
                'Makefile': makefile}
    if measurements is not None:
        genfiles['measurements.cpp'] = '#include "measurements.h"\n'+rtModelGen['measurementsFile'][0]
        genfiles['measurements.h'] = rtModelGen['measurementsFile'][1]
        genfiles['measurementsJacob.cpp'] = '#include "measurementsJacob.h"\n'+rtModelGen['measurementsJacobFile'][0]
        genfiles['measurementsJacob.h'] = rtModelGen['measurementsJacobFile'][1]
    exportpath = codegen.memoizeFiles(genfiles,prefix='rt_integrator__')

    # compile the code
    (ret, msgs) = subprocess_tee.call(['make',codegen.makeJobs()], cwd=exportpath)
    if ret != 0:
        raise Exception("integrator compilation failed:\n"+msgs)

    print 'loading '+exportpath+'/integrator.so'
    integratorLib = ctypes.cdll.LoadLibrary(exportpath+'/integrator.so')
    print 'loading '+exportpath+'/model.so'
    modelLib = ctypes.cdll.LoadLibrary(exportpath+'/model.so')
    return (integratorLib, modelLib, rtModelGen)
