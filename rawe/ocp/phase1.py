import ctypes
import os

from rawe.utils import codegen,pkgconfig,subprocess_tee
import writeAcadoOcpExport

def makeExportMakefile(cgOptions):
    rpath = None
    for blah in pkgconfig.call(['--libs','ocg2']).split(' '):
        blah = blah.strip()
        if len(blah) > 1 and blah[:2] == '-L':
            rpath = blah[2:]
            break
    assert rpath is not None, "couldn't detect the library path of ocg2 :("

    makefile = """\
CXX       = %(CXX)s
CXXFLAGS  = -O2 -fPIC -finline-functions -I. `pkg-config --cflags acado` `pkg-config --cflags ocg2`
LDFLAGS = -lm `pkg-config --libs acado` `pkg-config --libs ocg2`

CXX_SRC = export_ocp.cpp
OBJ = $(CXX_SRC:%%.cpp=%%.o)

.PHONY: clean all export_ocp.so
all : $(OBJ) export_ocp.so

%%.o : %%.cpp #acado.h
\t@echo CXX $@: $(CXX) $(CXXFLAGS) -c $< -o $@
\t@$(CXX) $(CXXFLAGS) -c $< -o $@

export_ocp.so::LDFLAGS+=-Wl,-rpath,%(rpath)s

export_ocp.so : $(OBJ)
\t@echo LD $@ : $(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)
\t@$(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)

clean :
\trm -f *.o *.so
""" % {'CXX':cgOptions['CXX'],'rpath':rpath}
    return makefile

# This writes and runs the ocp exporter, returning an exported OCP as a
# dictionary of files.
def runPhase1(ocp, cgOptions, acadoOptions, qpSolver):
    supportedQps = ['QP_OASES']
    assert qpSolver in supportedQps, "qp solver must be one of " + str(supportedQps)

    # write the ocp exporter cpp file
    genfiles = {'export_ocp.cpp':writeAcadoOcpExport.generateAcadoOcp(ocp, acadoOptions),
                'Makefile':makeExportMakefile(cgOptions)}
    exportpath = codegen.memoizeFiles(genfiles)

    # compile the ocp exporter
    (ret, msgs) = subprocess_tee.call(['make',codegen.makeJobs()], cwd=exportpath)
    if ret != 0:
        raise Exception("exportOcp phase 1 compilation failed:\n\n"+msgs)

    # load the ocp exporter
    lib = ctypes.cdll.LoadLibrary(os.path.join(exportpath, 'export_ocp.so'))
    lib.exportOcp.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_char_p]

    # run the ocp exporter
    def runOcpExporter(path):
        if qpSolver == 'QP_OASES':
            os.mkdir(os.path.join(path,'qpoases'))

        ret = lib.exportOcp(ocp._nk,
                            ctypes.c_double(ocp._ts),
                            ctypes.c_char_p(path))
        if ret != 0:
            raise Exception("call to export_ocp.so failed")
    return codegen.withTempdir(runOcpExporter)
