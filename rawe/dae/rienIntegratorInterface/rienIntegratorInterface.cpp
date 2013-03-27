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
			int numIntervals,
			double timestep,
			const char * integratorType,
			const char * integratorGrid,
			int numIntegratorSteps,
			int nx, int nz, int nu)
{
  SIMexport sim(numIntervals, timestep);

  // set INTEGRATOR_TYPE
  string integratorTypeStr(integratorType);
  map<string, int> itMap = makeIntegratorTypeMap();
  map<string,int>::const_iterator it = itMap.find(integratorTypeStr);
  if( it != itMap.end() ) {
    sim.set( INTEGRATOR_TYPE, it->second);
    //cout << "using integrator type \"" << integratorTypeStr << "\"\n";
  } else {
    cerr << "unrecognized integrator type \"" << integratorTypeStr << "\"\n";
    return -1;
  }

  // set INTEGRATOR_GRID
  if (integratorGrid != NULL) {
    string integratorGridStr(integratorGrid);
    map<string, int> igMap = makeIntegratorGridMap();
    it = igMap.find(integratorGridStr);
    if ( it != igMap.end() ) {
      sim.set( MEASUREMENT_GRID, it->second);
      //cout << "using measurement grid \"" << integratorGridStr << "\"\n";
    } else {
      cerr << "unrecognized measurement grid \"" << integratorGridStr << "\"\n";
      return -1;
    }
  }
  
  // set NUM_INTEGRATOR_STEPS
  sim.set( NUM_INTEGRATOR_STEPS, numIntegratorSteps );

  sim.setModel( "model", "rhs", "rhs_jac" );
  sim.setDimensions( nx, nx, nz, nu );

  sim.set( GENERATE_MAKE_FILE, false );

//    sim.addOutput( "out", "out_jac", 2 );
//    sim.setMeasurements( Meas );
//    sim.setTimingSteps( 10000 );
//    sim.exportAndRun( "externModel_export", "init_externModel.txt", "controls_externModel.txt" );

  sim.exportCode(genPath);
  return 0;
}
