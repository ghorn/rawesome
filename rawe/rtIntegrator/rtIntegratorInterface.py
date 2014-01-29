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

from ..utils import pkgconfig

def phase1src(dae,options,measurements):
    ret = '''\
#include <string>
#include <iostream>

#include <acado_toolkit.hpp>

extern "C"{
  int export_integrator( const char * genPath);
}

using namespace std;
USING_NAMESPACE_ACADO

int export_integrator( const char * genPath)
{
  const double timestep = 1.0;
  const int numIntervals = 1;
  SIMexport sim(numIntervals, timestep);

  // set INTEGRATOR_TYPE
  sim.set( INTEGRATOR_TYPE, %(INTEGRATOR_TYPE)s);

  // set NUM_INTEGRATOR_STEPS
  sim.set( NUM_INTEGRATOR_STEPS, %(NUM_INTEGRATOR_STEPS)s );

  // 0 == rhs()
  sim.setModel( "model", "rhs", "rhsJacob" );
  sim.setDimensions( %(nx)d, %(nx)d, %(nz)d, %(nup)d );
''' % {'INTEGRATOR_TYPE': options['INTEGRATOR_TYPE'],
       'NUM_INTEGRATOR_STEPS': options['NUM_INTEGRATOR_STEPS'],
       'nx':len(dae.xNames()),
       'nz':len(dae.zNames()),
       'nup':len(dae.uNames())+len(dae.pNames())}

    if measurements is not None:
        ret += '''
  // set MEASUREMENT_GRID
  // sim.set( MEASUREMENT_GRID, EQUIDISTANT_GRID );

  // set output/measurements
  DVector Meas( 1 );
  Meas( 0 ) = 1;
  sim.addOutput("measurements", "measurementsJacob", %(outputDimension)d, Meas);
''' % {'outputDimension':C.densify(measurements).size()}
    ret +='''
  sim.set( GENERATE_MAKE_FILE, false );

  return sim.exportCode(genPath);
}
'''
    return ret

def phase1makefile():
    rpathAcado = pkgconfig.getRpath('acado')
    return '''\
CXX       = g++
CXXFLAGS  = -Wall
CXXFLAGS += -g
CXXFLAGS += -O2
CXXFLAGS += -fPIC
CXXFLAGS += -I..
CXXFLAGS += -Wno-unused-function
#CXXFLAGS += -std=c++11

LDFLAGS = -lstdc++

CXXFLAGS += `pkg-config --cflags acado`
LDFLAGS  += `pkg-config --libs   acado`

CPP_SRC = export_integrator.cpp

.PHONY: clean all export_integrator.so
all : $(OBJ) export_integrator.so

%%.o : %%.cpp
	@echo CPP $@ #: $(CXX) $(CXXFLAGS) -c $< -o $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@

OBJ = $(CPP_SRC:%%.cpp=%%.o)

export_integrator.so::LDFLAGS+=-Wl,-rpath,%(rpathAcado)s

export_integrator.so : $(OBJ)
	@echo LD $@ #: $(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)
	@$(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)

clean :
	rm -f *.o *.so
''' % {'rpathAcado':rpathAcado}
