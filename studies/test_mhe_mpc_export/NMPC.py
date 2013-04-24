import rawe
import casadi as C

def makeNmpc(dae,N,dt):
    from rawe.ocp import Ocp
    mpc = Ocp(dae, N=N, ts=dt)
    
    xref = C.veccat( [0,0] )
    uref = C.veccat( [0] )

    acadoOpts=[('HESSIAN_APPROXIMATION','GAUSS_NEWTON'),
               ('DISCRETIZATION_TYPE','MULTIPLE_SHOOTING'),
               ('QP_SOLVER','QP_QPOASES'),
               ('HOTSTART_QP','NO'),
               ('INTEGRATOR_TYPE','INT_IRK_GL2'),
               ('NUM_INTEGRATOR_STEPS',str(40*N)),
               ('IMPLICIT_INTEGRATOR_NUM_ITS','3'),
               ('IMPLICIT_INTEGRATOR_NUM_ITS_INIT','0'),
               ('LINEAR_ALGEBRA_SOLVER','HOUSEHOLDER_QR'),
               ('UNROLL_LINEAR_SOLVER','NO'),
               ('IMPLICIT_INTEGRATOR_MODE','IFTR'),
               ('SPARSE_QP_SOLUTION','CONDENSING'),
#               ('SPARSE_QP_SOLUTION','FULL_CONDENSING_U2'),
#               ('AX_NUM_QP_ITERATIONS','30'),
               ('FIX_INITIAL_STATE','YES'),
               ('GENERATE_TEST_FILE','NO'),
               ('GENERATE_SIMULINK_INTERFACE','NO'),
               ('GENERATE_MAKE_FILE','NO'),
               ('CG_USE_C99','YES')]

    mpc.minimizeLsq(C.veccat([mpc['x'],mpc['v'],mpc['u']]))
    mpc.minimizeLsqEndTerm(C.veccat([mpc['x'],mpc['v']]))

    cgOpts = {'CXX':'g++', 'CC':'gcc'}
    mpcRT = mpc.exportCode(codegenOptions=cgOpts,acadoOptions=acadoOpts)
    
    
    mpcRT.S[0,0] = 1.0#/(xRms*N)**2
    mpcRT.S[1,1] = 1.0#/(vRms*N)**2
    mpcRT.S[2,2] = 1.0#/(fRms*N)**2

    mpcRT.SN[0,0] = 1.0#/(xRms*N)**2
    mpcRT.SN[1,1] = 1.0#/(vRms*N)**2
    
    return mpcRT

if __name__=='__main__':
    from highwind_carousel_conf import conf
    dae = rawe.models.carousel(conf)

    OcpRt = makeNmpc(dae,10,0.1)
