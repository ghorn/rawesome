from ..utils import pkgconfig

def phase1src():
    return '''\
#include <string>
#include <iostream>

#include <acado_code_generation.hpp>
#include <acado/utils/acado_types.hpp>

extern "C"{
  int makeRienIntegrator( const char * genPath,
			  int numIntervals,
			  double timestep,
			  const char * integratorType,
			  const char * integratorGrid,
			  int numIntegratorSteps,
			  int nx, int nz, int nu);
}

using namespace std;
USING_NAMESPACE_ACADO

map<string, int> makeIntegratorTypeMap(void){
  map<string, int> it;
  it["INT_EX_EULER"] = INT_EX_EULER;
  it["INT_RK2"] = INT_RK2;
  it["INT_RK3"] = INT_RK3;
  it["INT_RK4"] = INT_RK4;
  it["INT_IRK_GL2"] = INT_IRK_GL2;
  it["INT_IRK_GL4"] = INT_IRK_GL4;
  it["INT_IRK_GL6"] = INT_IRK_GL6;
  it["INT_IRK_GL8"] = INT_IRK_GL8;
  it["INT_IRK_RIIA1"] = INT_IRK_RIIA1;
  it["INT_IRK_RIIA3"] = INT_IRK_RIIA3;
  it["INT_IRK_RIIA5"] = INT_IRK_RIIA5;
  it["INT_DIRK3"] = INT_DIRK3;
  it["INT_DIRK4"] = INT_DIRK4;
  it["INT_DIRK5"] = INT_DIRK5;
//  it["INT_RK12"] = INT_RK12;
//  it["INT_RK23"] = INT_RK23;
//  it["INT_RK45"] = INT_RK45;
//  it["INT_RK78"] = INT_RK78;
//  it["INT_BDF"] = INT_BDF;
//  it["INT_DISCRETE"] = INT_DISCRETE;
//  it["INT_LYAPUNOV45"] = INT_LYAPUNOV45;
  return it;
}

map<string, int> makeIntegratorGridMap(void){
  map<string, int> ig;
  ig["EQUIDISTANT_SUBGRID"] = EQUIDISTANT_SUBGRID;
  ig["EQUIDISTANT_GRID"]    = EQUIDISTANT_GRID;
  ig["ONLINE_GRID"]	    = ONLINE_GRID;
  return ig;
}

int makeRienIntegrator( const char * genPath,
			const int numIntervals,
			const double timestep,
			const char * integratorType,
			const char * integratorGrid,
			const int numIntegratorSteps,
			const int nx, const int nz, const int nu)
{
  SIMexport sim(numIntervals, timestep);

  // set INTEGRATOR_TYPE
  string integratorTypeStr(integratorType);
  map<string, int> itMap = makeIntegratorTypeMap();
  map<string,int>::const_iterator it = itMap.find(integratorTypeStr);
  if( it != itMap.end() ) {
    sim.set( INTEGRATOR_TYPE, it->second);
    //cout << "using integrator type \\"" << integratorTypeStr << "\\"\\n";
  } else {
    cerr << "unrecognized integrator type \\"" << integratorTypeStr << "\\"\\n";
    return -1;
  }

  // set INTEGRATOR_GRID
  if (integratorGrid != NULL) {
    string integratorGridStr(integratorGrid);
    map<string, int> igMap = makeIntegratorGridMap();
    it = igMap.find(integratorGridStr);
    if ( it != igMap.end() ) {
      sim.set( MEASUREMENT_GRID, it->second);
      //cout << "using measurement grid \\"" << integratorGridStr << "\\"\\n";
    } else {
      cerr << "unrecognized measurement grid \\"" << integratorGridStr << "\\"\\n";
      return -1;
    }
  }
  
  // set NUM_INTEGRATOR_STEPS
  sim.set( NUM_INTEGRATOR_STEPS, numIntegratorSteps );

  // 0 == rhs()
  sim.setModel( "model", "rhs", "rhs_jac" );
  sim.setDimensions( nx, nx, nz, nu );

  sim.set( GENERATE_MAKE_FILE, false );

//    sim.addOutput( "out", "out_jac", 2 );
//    sim.setMeasurements( Meas );
//    sim.setTimingSteps( 10000 );
//    sim.exportAndRun( "externModel_export", "init_externModel.txt", "controls_externModel.txt" );

  return sim.exportCode(genPath);
}
'''

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

CPP_SRC = rienIntegratorInterface.cpp

.PHONY: clean all rienIntegratorInterface.so
all : $(OBJ) rienIntegratorInterface.so

%%.o : %%.cpp
	@echo CPP $@ #: $(CXX) $(CXXFLAGS) -c $< -o $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@

OBJ = $(CPP_SRC:%%.cpp=%%.o)

rienIntegratorInterface.so::LDFLAGS+=-Wl,-rpath,%(rpathAcado)s

rienIntegratorInterface.so : $(OBJ)
	@echo LD $@ #: $(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)
	@$(CXX) -shared -o $@ $(OBJ) $(LDFLAGS)

clean :
	rm -f *.o *.so
''' % {'rpathAcado':rpathAcado}
