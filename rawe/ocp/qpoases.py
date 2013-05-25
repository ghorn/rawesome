import os
from rawe.utils import pkgconfig, codegen, subprocess_tee

def mkMakefile(cgOptions, qposrc):
    qposrc = ' \\\n'.join(['\t'+os.path.join('qpoases', q.split('qpoases'+os.sep)[1]) for q in qposrc])
    if cgOptions['hideSymbols']:
        c_visibility = ' -fvisibility=hidden'
        cxx_visibility = ' -fvisibility=hidden -fvisibility-inlines-hidden'
    else:
        c_visibility = ''
        cxx_visibility = ''

    makefile = """\
CXX      = %(CXX)s
CC       = %(CC)s
CXXFLAGS = -O3 -fPIC -finline-functions%(cxx_visibility)s
CFLAGS   = -O3 -fPIC -finline-functions%(c_visibility)s

#CFLAGS   += -Wall -Wextra
#CXXFLAGS += -Wall -Wextra
#CXXFLAGS += -DPC_DEBUG # make qpoases print out a bunch of debugging info

LDFLAGS = -lm

UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
#	LDFLAGS += -L/opt/local/lib
#	INCLUDES += -I/opt/local/include
#	INCLUDES += -isystem /usr/local/include
else
	LDFLAGS += -lrt
endif


CXX_SRC = \\
\trhs.cpp \\
\trhsJacob.cpp \\
\tacado_external_functions.cpp \\
%(qpo_src)s \\
\tqpoases/solver.cpp

C_SRC = \\
\tworkspace.c \\
\tpython_interface.c \\
\tmodel.c \\
\tacado_integrator.c \\
\tacado_solver.c \\
\tacado_auxiliary_functions.c

QPO_INC = \\
\t-I. \\
\t-Iqpoases \\
\t-Iqpoases/INCLUDE \\
\t-Iqpoases/SRC

CXXFLAGS += $(QPO_INC)

CXX_OBJ = $(CXX_SRC:%%.cpp=%%.o)
C_OBJ = $(C_SRC:%%.c=%%.o)

HEADERS = \\
\tqpoases/solver.hpp \\
\tacado_auxiliary_functions.h \\
\tacado_common.h

.PHONY: clean all ocp.a ocp.so
all : $(CXX_OBJ) $(C_OBJ) ocp.a ocp.so ocp.o

$(CXX_OBJ) : %%.o : %%.cpp $(HEADERS)
\t@echo CXX $@: $(CXX) $(CXXFLAGS) -c $< -o $@
\t@$(CXX) $(CXXFLAGS) -c $< -o $@

model.o :: CFLAGS += -Wall -Wextra -Werror -Wno-unused-variable
workspace.o :: CFLAGS += -Wall -Wextra -Werror
python_interface.o :: CFLAGS += -Wall -Wextra -Werror
$(C_OBJ) : %%.o : %%.c $(HEADERS)
\t@echo CC $@: $(CC) $(CFLAGS) -c $< -o $@
\t@$(CC) $(CFLAGS) -c $< -o $@

ocp.so : $(CXX_OBJ) $(C_OBJ)
\t@echo LD $@: $(CXX) -shared -o $@ $(CXX_OBJ) $(C_OBJ) $(LDFLAGS)
\t@$(CXX) -shared -o $@ $(CXX_OBJ) $(C_OBJ) $(LDFLAGS)

ocp.a : $(CXX_OBJ) $(C_OBJ)
\t@echo AR $@ : ar r $@ $?
\t@ar r $@ $?

ocp.o : $(CXX_OBJ) $(C_OBJ)
\t@echo ld $@ : ld -r $? -o $@
\t@ld -r $? -o $@

clean :
\t@echo rm -f ocp.a $(CXX_OBJ) $(C_OBJ) ocp.so
\t@rm -f ocp.a ocp.so ocp.o $(CXX_OBJ) $(C_OBJ)
""" % {'CXX':cgOptions['CXX'], 'CC':cgOptions['CC'],
       'c_visibility':c_visibility,
       'cxx_visibility':cxx_visibility,
       'qpo_src':qposrc}
    return makefile


def exportPhase2(cgOptions, phase1src):
    # call pkg-config to get qpoases source and includes
    qpoStuff = {}
    for name in ['qpOASESsrc', 'qpOASESinc']:
        qpoStuff[name] = pkgconfig.call(['--variable',name,'acado'])
    qpoStuff['qpOASESsrc'] = qpoStuff['qpOASESsrc'].split(' ')

    # get qpoases source as file dictionary
    qpoSrcPath = os.path.join(qpoStuff['qpOASESsrc'][0].split('qpoases')[0], 'qpoases')
    phase2src = codegen.directoryToDict(qpoSrcPath)

    # merge qpoases source with phase 1 output
    def mergeAll(srcdict,destdict):
        for name,src in srcdict.items():
            if isinstance(src,dict):
                if name not in destdict:
                    destdict[name] = {}
                assert isinstance(destdict[name], dict), "dictionary merge failed, source was directory but destination was a file"
                destdict[name] = mergeAll(src,destdict[name])
            else:
                destdict[name] = src
        return destdict
    genfiles = mergeAll(phase1src, {'qpoases':phase2src})

    # add makefile
    genfiles['Makefile'] = mkMakefile(cgOptions, qpoStuff['qpOASESsrc'])
    genfiles['workspace.c'] ='''\
#include "acado_common.h"
ACADOworkspace acadoWorkspace;
ACADOvariables acadoVariables;
'''

    # write all this
    exportpath = codegen.memoizeFiles(genfiles,prefix=cgOptions['hashPrefix']+'__')

    # compile!
    (ret, msgs) = subprocess_tee.call(['make',codegen.makeJobs()], cwd=exportpath)
    if ret != 0:
        raise Exception("ocp compilation failed:\n\n"+msgs)

    # return shared object
    return exportpath
