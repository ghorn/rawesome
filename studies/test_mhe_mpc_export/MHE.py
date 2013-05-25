import rawe
import casadi as C

intOpts = rawe.RtIntegratorOptions()
intOpts['INTEGRATOR_TYPE'] = 'INT_IRK_GL2'
intOpts['NUM_INTEGRATOR_STEPS'] = 40
intOpts['IMPLICIT_INTEGRATOR_NUM_ITS'] = 3
intOpts['IMPLICIT_INTEGRATOR_NUM_ITS_INIT'] = 0
intOpts['LINEAR_ALGEBRA_SOLVER'] = 'HOUSEHOLDER_QR'
intOpts['UNROLL_LINEAR_SOLVER'] = False
intOpts['IMPLICIT_INTEGRATOR_MODE'] = 'IFTR'

def makeMhe(dae,N,dt):
    from rawe.ocp import Ocp
    mhe = Ocp(dae, N=N, ts=dt)
    
    ocpOpts = rawe.OcpExportOptions()
    ocpOpts['HESSIAN_APPROXIMATION'] = 'GAUSS_NEWTON'
    ocpOpts['DISCRETIZATION_TYPE'] = 'MULTIPLE_SHOOTING'
    ocpOpts['QP_SOLVER'] = 'QP_QPOASES'
    ocpOpts['HOTSTART_QP'] = False
    ocpOpts['SPARSE_QP_SOLUTION'] = 'CONDENSING'
#   ocpOpts['SPARSE_QP_SOLUTION'] = 'FULL_CONDENSING_U2'
#   ocpOpts['AX_NUM_QP_ITERATIONS'] = '30'
    ocpOpts['FIX_INITIAL_STATE'] = False
               
#    mhe.minimizeLsq(C.veccat([mhe['x'],mhe['u']]))
#    mhe.minimizeLsqEndTerm(C.veccat([mhe['x']]))
    mhe.minimizeLsq(mhe['measurements'])
    mhe.minimizeLsqEndTerm(mhe['measurementsN'])

    cgOpts = {'CXX':'g++', 'CC':'gcc'}
    print "makeMhe calling rawe.OcpRT"
    return rawe.OcpRT(mhe, ocpOptions=ocpOpts, integratorOptions=intOpts,
                       codegenOptions=cgOpts)
