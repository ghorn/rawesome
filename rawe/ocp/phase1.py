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

from rawe.utils import codegen,pkgconfig,subprocess_tee
import writeAcadoOcpExport

def makeExportMakefile(phase1Options):
    rpathAcado = pkgconfig.getRpath('acado')
    rpathOcg2 = pkgconfig.getRpath('ocg2')

    makefile = """\
CXX       = %(CXX)s
CXXFLAGS  = -O2 -fPIC -finline-functions -I. `pkg-config --cflags acado` `pkg-config --cflags ocg2`
LDFLAGS = -lm `pkg-config --libs acado` `pkg-config --libs ocg2`

CXX_SRC = export_ocp.cpp
OBJ = $(CXX_SRC:%%.cpp=%%.o)

.PHONY: clean all export_ocp.so run_export
all : $(OBJ) export_ocp.so run_export

%%.o : %%.cpp
\t@echo CXX $@: $(CXX) $(CXXFLAGS) -c $< -o $@
\t@$(CXX) $(CXXFLAGS) -c $< -o $@

export_ocp.so::LDFLAGS+=-Wl,-rpath,%(rpathAcado)s
export_ocp.so::LDFLAGS+=-Wl,-rpath,%(rpathOcg2)s

export_ocp.so : $(OBJ)
\t@echo LD $@ : $(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)
\t@$(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)

run_export : export_ocp.so
	@echo $(CXX) export_ocp.so -Wl,-rpath,. -o $@
	@$(CXX) run_export.cpp export_ocp.so -Wl,-rpath,. -o $@

clean :
\trm -f *.o *.so
""" % {'CXX':phase1Options['CXX'],'rpathAcado':rpathAcado,'rpathOcg2':rpathOcg2}
    return makefile

# This writes and runs the ocp exporter, returning an exported OCP as a
# dictionary of files.
def runPhase1(ocp, phase1Options, integratorOptions, ocpOptions):
    # write the ocp exporter cpp file
    genfiles = {'export_ocp.cpp':writeAcadoOcpExport.generateAcadoOcp(ocp, integratorOptions, ocpOptions),
                'Makefile':makeExportMakefile(phase1Options)}
    # add a file which just runs the export in the current directory
    genfiles['run_export.cpp'] = '''\
extern "C" int exportOcp(const char * exportDir);

int main(void){
    return exportOcp(".");
}
'''
    exportpath = codegen.memoizeFiles(genfiles,prefix=phase1Options['hashPrefix']+'_phase1__')

    # compile the ocp exporter
    (ret, msgs) = subprocess_tee.call(['make',codegen.makeJobs()], cwd=exportpath)
    if ret != 0:
        raise Exception("exportOcp phase 1 compilation failed:\n\n"+msgs)

    # run the ocp exporter
    def runOcpExporter(path):
        # load the ocp exporter
        lib = ctypes.cdll.LoadLibrary(os.path.join(exportpath, 'export_ocp.so'))

        if ocpOptions['QP_SOLVER'] == 'QP_QPOASES':
            os.mkdir(os.path.join(path,'qpoases'))

        ret = lib.exportOcp(ctypes.c_char_p(path))
        if ret != 0:
            raise Exception("call to export_ocp.so failed")

    def callExporterInProcess(q):
        try:
            q.put(codegen.withTempdir(runOcpExporter))
        finally:
            q.put(None)

    q = Queue()
    p = Process(target=callExporterInProcess,args=(q,))
    p.start()
    ret = q.get()
    p.join()
    assert (0 == p.exitcode) and (ret is not None), \
        "error exporting ocp, see stdout/stderr above"
    return ret
