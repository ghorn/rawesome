import subprocess
import ctypes
import os

from ..utils.codegen import memoizeFiles,withTempdir
import writeAcadoOcpExport

def makeMakefile(CXX):
    makefile = """\
CXX       = %(CXX)s
CXXFLAGS  = -O2 -fPIC -finline-functions -I. `pkg-config --cflags acado` -I/home/ghorn/acado/experimental/mvukov/ocg2
LDFLAGS = -lm `pkg-config --libs acado` -L/home/ghorn/acado/build/experimental/mvukov/ocg2/lib -locg2

CXX_SRC = export_ocp.cpp
OBJ = $(CXX_SRC:%%.cpp=%%.o)

.PHONY: clean libs obj

all : export_ocp.so

%%.o : %%.cpp #acado.h
\t@echo CXX $@: $(CXX) $(CXXFLAGS) -c $< -o $@
\t@$(CXX) $(CXXFLAGS) -c $< -o $@

%%.so : $(OBJ)
\t@echo LD $@: $(CXX) -shared -Wl,-soname,$@ -o $@ $(OBJ) $(LDFLAGS)
\t@$(CXX)   -shared -Wl,-soname,$@ -o $@ $(OBJ) $(LDFLAGS)

# aliases
obj : $(OBJ)
clean :
\trm -f *.o *.so
""" % {'CXX':CXX}
    return makefile

def exportOcp(ocp, CXX):
    # write the ocp exporter cpp file
    genfiles = [('export_ocp.cpp', writeAcadoOcpExport.generateAcadoOcp(ocp)),
                ('Makefile',makeMakefile(CXX))]
    exportpath = memoizeFiles(genfiles)

    # compile the ocp exporter
    p = subprocess.Popen(['make'], stdout=subprocess.PIPE, cwd=exportpath)
    ret = p.wait()
    if ret != 0:
        print "stdout: "+p.stdout.read()
        raise Exception("integrator compilation failed, return code "+str(ret))
    
    # load the ocp exporter
    Ni = 5
    print os.path.join(exportpath, 'export_ocp.so')
    lib = ctypes.cdll.LoadLibrary(os.path.join(exportpath, 'export_ocp.so'))
    
    # run the ocp exporter
    def runOcpExporter(path):
        ret = lib.exportOcp(ocp._nk,
                            Ni,
                            ctypes.c_double(ocp._ts),
                            ctypes.c_char_p(path))
        if ret != 0:
            raise Exception("call to export_ocp.so failed")
    files = withTempdir(runOcpExporter)
    print "exported files: "+str([name for (name,_) in files])
