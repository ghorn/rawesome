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

#include <acado_code_generation.hpp>

extern "C"{
  int export_integrator( const char * genPath);
}

using namespace std;
USING_NAMESPACE_ACADO

int export_integrator( const char * genPath)
{
    string path( genPath );
    string _stdout = path + "/_stdout.txt";
    
    std::ofstream out( _stdout.c_str() );
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); // redirect std::cout to the text file

    Logger::instance().setLogLevel( LVL_DEBUG );

  const double timestep = 1.0;
  const int numIntervals = 1;
  SIMexport sim(numIntervals, timestep);

  // set INTEGRATOR_TYPE
  sim.set( INTEGRATOR_TYPE, %(INTEGRATOR_TYPE)s);

  // set NUM_INTEGRATOR_STEPS
  sim.set( NUM_INTEGRATOR_STEPS, %(NUM_INTEGRATOR_STEPS)s );
  
  sim.set( DYNAMIC_SENSITIVITY, %(dynamic_sensitivity)s );

  // 0 == rhs()
  sim.setModel( "model", "rhs", "rhsJacob" );
  sim.setDimensions( %(nx)d, %(nx)d, %(nz)d, %(nu)d, %(nod)d, 0 );
''' % {'INTEGRATOR_TYPE': options['INTEGRATOR_TYPE'],
       'NUM_INTEGRATOR_STEPS': options['NUM_INTEGRATOR_STEPS'],
       'nx': len(dae.xNames()),
       'nz': len(dae.zNames()),
       'nu': len(dae.uNames()),
       'nod': len(dae.pNames()),
       'dynamic_sensitivity': options['DYNAMIC_SENSITIVITY']}

    if measurements is not None:
        ret += '''
  // set MEASUREMENT_GRID
  // sim.set( MEASUREMENT_GRID, EQUIDISTANT_GRID );

  // set output/measurements
  DVector Meas( 1 );
  Meas( 0 ) = 1;
  sim.addOutput("measurements", "measurementsJacob", %(outputDimension)d, Meas);
''' % {'outputDimension':C.dense(measurements).size()}
    ret +='''
  sim.set( GENERATE_MAKE_FILE, false );
    
    returnValue status = sim.exportCode(genPath);
    
    std::cout.rdbuf(coutbuf); //reset to standard output again

    return status;
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
