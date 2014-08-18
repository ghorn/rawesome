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

from rawe.utils import codegen, subprocess_tee

def mkMakefile( cgOptions ):
    if cgOptions['hideSymbols']:
        c_visibility = ' -fvisibility=hidden'
        cxx_visibility = ' -fvisibility=hidden -fvisibility-inlines-hidden'
    else:
        c_visibility = ''
        cxx_visibility = ''

    makefile = """\
CXX      = %(CXX)s
CC       = %(CC)s
CXXFLAGS = %(CXXFLAGS)s%(cxx_visibility)s
CFLAGS   = %(CFLAGS)s%(c_visibility)s

#CXXFLAGS += -DPC_DEBUG # make qpoases print out a bunch of debugging info

LDFLAGS = -lm

# qpDUNES library
LDFLAGS += -lqpdunes

UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
#    LDFLAGS += -L/opt/local/lib
#    INCLUDES += -I/opt/local/include
    CFLAGS += -I/usr/local/include
else
    LDFLAGS += -lrt
endif

CFLAGS += -I/usr/local/include/qpdunes/
CFLAGS += -I/usr/local/include/qpdunes/qp

CXXFLAGS += -I/usr/local/include/qpdunes/
CXXFLAGS += -I/usr/local/include/qpdunes/qp

CXX_SRC = \\
\teverything.cpp

C_SRC = \\
\tpython_interface.c \\
\tacado_auxiliary_functions.c

CXX_OBJ = $(CXX_SRC:%%.cpp=%%.o)
C_OBJ = $(C_SRC:%%.c=%%.o)

HEADERS = \\
\tacado_common.h \\
\tacado_auxiliary_functions.h

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
\t@echo ld $@ : ld -r $? $(LDFLAGS) -o $@
\t@ld -r $? -o $@

clean :
\t@echo rm -f ocp.a $(CXX_OBJ) $(C_OBJ) ocp.so
\t@rm -f ocp.a ocp.so ocp.o $(CXX_OBJ) $(C_OBJ)
""" % {'CXX': cgOptions['CXX'],
       'CC': cgOptions['CC'],
       'CXXFLAGS': cgOptions['CXXFLAGS'],
       'CFLAGS': cgOptions['CFLAGS'],
       'c_visibility': c_visibility,
       'cxx_visibility': cxx_visibility}
    return makefile

def exportPhase2(cgOptions, phase1src):
    
    genfiles = phase1src

    # add makefile
    genfiles['Makefile'] = mkMakefile( cgOptions )
    genfiles['workspace.c'] ='''\
#include "acado_common.h"
ACADOworkspace acadoWorkspace;
ACADOvariables acadoVariables;
'''

    # write things in user specified directory if 'export_without_build_path' is not None
    # then return without building
    if cgOptions['export_without_build_path'] is not None:
        codegen.writeDifferentFiles(cgOptions['export_without_build_path'], genfiles)
        return

    # write things in user specified directory if 'force_export_path' is not None
    # then compile as normal
    if cgOptions['force_export_path'] is not None:
        codegen.writeDifferentFiles(cgOptions['force_export_path'], genfiles)
        exportpath = cgOptions['force_export_path']
    # otherwise write things in memoized directory, then compile as normal
    else:
        exportpath = codegen.memoizeFiles(genfiles,prefix=cgOptions['hashPrefix']+'__')

    # compile!
    (ret, msgs) = subprocess_tee.call(['make',codegen.makeJobs()], cwd=exportpath)
    if ret != 0:
        raise Exception("ocp compilation failed:\n\n"+msgs)

    # return shared object
    return exportpath
