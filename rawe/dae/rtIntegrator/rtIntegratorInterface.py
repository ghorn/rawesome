import casadi as C

from ...utils import pkgconfig

def phase1src(dae,options,measurements):
    ret = '''\
#include <string>
#include <iostream>

#include <acado_code_generation.hpp>
#include <acado/utils/acado_types.hpp>

extern "C"{
  int makeRtIntegrator( const char * genPath);
}

using namespace std;
USING_NAMESPACE_ACADO

int makeRtIntegrator( const char * genPath)
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
  sim.set( MEASUREMENT_GRID, EQUIDISTANT_GRID );

  // set output/measurements
  sim.addOutput( "measurements", "measurementsJacob", %(outputDimension)d );
  Vector Meas(1);
  Meas(0) = 1;
  sim.setMeasurements( Meas );
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

CPP_SRC = rtIntegratorInterface.cpp

.PHONY: clean all rtIntegratorInterface.so
all : $(OBJ) rtIntegratorInterface.so

%%.o : %%.cpp
	@echo CPP $@ #: $(CXX) $(CXXFLAGS) -c $< -o $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@

OBJ = $(CPP_SRC:%%.cpp=%%.o)

rtIntegratorInterface.so::LDFLAGS+=-Wl,-rpath,%(rpathAcado)s

rtIntegratorInterface.so : $(OBJ)
	@echo LD $@ #: $(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)
	@$(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)

clean :
	rm -f *.o *.so
''' % {'rpathAcado':rpathAcado}
