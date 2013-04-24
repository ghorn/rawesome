
ocg_interface = '''\
#include "acado_common.h"
#include "acado_auxiliary_functions.h"

#include <string.h>

void py_initialize(void){
  memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
  memset(&acadoVariables, 0, sizeof( acadoVariables ));
  initializeSolver();
}

int memcpyMat(real_t * const dest, real_t const * const src,
               const int nr1, const int nc1,
               const int nr2, const int nc2){
  int exactMatch = (nr1 == nr2) && (nc1 == nc2);
  int transposeMatch = (nr1 == nc2) && (nr2 == nc1);
  if (exactMatch || transposeMatch){
    memcpy(dest, src, sizeof(real_t)*nr1*nc1);
    return 0;
  } else {
    return 1;
  }
}

int py_set_x(real_t * val, const int nr, const int nc){
  return memcpyMat(acadoVariables.x, val, nr, nc, ACADO_N + 1, ACADO_NX); }
int py_get_x(real_t * val, const int nr, const int nc){
  return memcpyMat(val, acadoVariables.x, nr, nc, ACADO_N + 1, ACADO_NX); }

int py_set_u(real_t * val, const int nr, const int nc){
  return memcpyMat(acadoVariables.u, val, nr, nc, ACADO_N, ACADO_NU); }
int py_get_u(real_t * val, const int nr, const int nc){
  return memcpyMat(val, acadoVariables.u, nr, nc, ACADO_N, ACADO_NU); }

#if ACADO_NP
int py_set_p(real_t * val, const int nr, const int nc){
  return memcpyMat(acadoVariables.p, val, nr, nc, ACADO_NP, 1); }
int py_get_p(real_t * val, const int nr, const int nc){
  return memcpyMat(val, acadoVariables.p, nr, nc, ACADO_NP, 1); }
#endif

int py_set_y(real_t * val, const int nr, const int nc){
  return memcpyMat(acadoVariables.y, val, nr, nc, ACADO_N, ACADO_NY); }
int py_get_y(real_t * val, const int nr, const int nc){
  return memcpyMat(val, acadoVariables.y, nr, nc, ACADO_N, ACADO_NY); }

#if ACADO_NYN
int py_set_yN(real_t * val, const int nr, const int nc){
  return memcpyMat(acadoVariables.yN, val, nr, nc, ACADO_NYN, 1); }
int py_get_yN(real_t * val, const int nr, const int nc){
  return memcpyMat(val, acadoVariables.yN, nr, nc, ACADO_NYN, 1); }
#endif /* ACADO_NYN */

#if ACADO_INITIAL_STATE_FIXED
int py_set_x0(real_t * val, const int nr, const int nc){
  return memcpyMat(acadoVariables.x0, val, nr, nc, ACADO_NX, 1); }
int py_get_x0(real_t * val, const int nr, const int nc){
  return memcpyMat(val, acadoVariables.x0, nr, nc, ACADO_NX, 1); }
#endif /* ACADO_INITIAL_STATE_FIXED */

#if ACADO_WEIGHTING_MATRICES_TYPE == 1
int py_set_S(real_t * val, const int nr, const int nc){
  return memcpyMat(acadoVariables.S, val, nr, nc, ACADO_NY, ACADO_NY); }
int py_get_S(real_t * val, const int nr, const int nc){
  return memcpyMat(val, acadoVariables.S, nr, nc, ACADO_NY, ACADO_NY); }
int py_set_SN(real_t * val, const int nr, const int nc){
  return memcpyMat(acadoVariables.SN, val, nr, nc, ACADO_NYN, ACADO_NYN); }
int py_get_SN(real_t * val, const int nr, const int nc){
  return memcpyMat(val, acadoVariables.SN, nr, nc, ACADO_NYN, ACADO_NYN); }
#elif ACADO_WEIGHTING_MATRICES_TYPE == 2
int py_set_S(real_t * val, const int nr, const int nc){
  return memcpyMat(acadoVariables.S, val, nr, nc, ACADO_N * ACADO_NY, ACADO_NY); }
int py_get_S(real_t * val, const int nr, const int nc){
  return memcpyMat(val, acadoVariables.S, nr, nc, ACADO_N * ACADO_NY, ACADO_NY); }
int py_set_SN(real_t * val, const int nr, const int nc){
  return memcpyMat(acadoVariables.SN, val, nr, nc, ACADO_NYN, ACADO_NYN); }
int py_get_SN(real_t * val, const int nr, const int nc){
  return memcpyMat(val, acadoVariables.SN, nr, nc, ACADO_NYN, ACADO_NYN); }
#endif /* ACADO_WEIGHTING_MATRICES_TYPE */


/** Number of control/estimation intervals. */
int py_get_ACADO_N(void){ return ACADO_N; }
/** Number of differential variables. */
int py_get_ACADO_NX(void){ return ACADO_NX; }
/** Number of differential derivative variables. */
int py_get_ACADO_NXD(void){ return ACADO_NXD; }
/** Number of algebraic variables. */
int py_get_ACADO_NXA(void){ return ACADO_NXA; }
/** Number of control variables. */
int py_get_ACADO_NU(void){ return ACADO_NU; }
/** Number of parameters (which are NOT optimization variables). */
int py_get_ACADO_NP(void){ return ACADO_NP; }
/** Number of references/measurements per node on the first N nodes. */
int py_get_ACADO_NY(void){ return ACADO_NY; }
/** Number of references/measurements on the last (N + 1)st node. */
int py_get_ACADO_NYN(void){ return ACADO_NYN; }

/** qpOASES QP solver indicator. */
int py_get_ACADO_QPOASES(void){ return ACADO_QPOASES; }
/** FORCES QP solver indicator.*/
int py_get_ACADO_FORCES(void){ return ACADO_FORCES; }
/** qpDUNES QP solver indicator.*/
int py_get_ACADO_QPDUNES(void){ return ACADO_QPDUNES; }
/** Indicator for determining the QP solver used by the ACADO solver code. */
int py_get_ACADO_QP_SOLVER(void){ return ACADO_QP_SOLVER; }
/** Indicator for fixed initial state. */
int py_get_ACADO_INITIAL_STATE_FIXED(void){ return ACADO_INITIAL_STATE_FIXED; }
/** Indicator for type of fixed weighting matrices. */
int py_get_ACADO_WEIGHTING_MATRICES_TYPE(void){ return ACADO_WEIGHTING_MATRICES_TYPE; }
/** Flag indicating whether constraint values are hard-coded or not. */
int py_get_ACADO_HARDCODED_CONSTRAINT_VALUES(void){ return ACADO_HARDCODED_CONSTRAINT_VALUES; }
'''
