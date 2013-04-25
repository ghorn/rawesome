import rawe
import casadi as C

def makeMhe(dae,N,dt):
    from rawe.ocp import Ocp
    mhe = Ocp(dae, N=N, ts=dt)
    
    intOpts = {'numsteps':40}
    intOpts = [('INTEGRATOR_TYPE','INT_IRK_GL2'),
               ('NUM_INTEGRATOR_STEPS',str(intOpts['numsteps']*N)),
               ('IMPLICIT_INTEGRATOR_NUM_ITS','3'),
               ('IMPLICIT_INTEGRATOR_NUM_ITS_INIT','0'),
               ('LINEAR_ALGEBRA_SOLVER','HOUSEHOLDER_QR'),
               ('UNROLL_LINEAR_SOLVER','NO'),
               ('IMPLICIT_INTEGRATOR_MODE','IFTR')]
    
    acadoOpts=[('HESSIAN_APPROXIMATION','GAUSS_NEWTON'),
               ('DISCRETIZATION_TYPE','MULTIPLE_SHOOTING'),
               ('QP_SOLVER','QP_QPOASES'),
               ('HOTSTART_QP','NO'),
               ('SPARSE_QP_SOLUTION','CONDENSING'),
#               ('SPARSE_QP_SOLUTION','FULL_CONDENSING_U2'),
#               ('AX_NUM_QP_ITERATIONS','30'),
               ('FIX_INITIAL_STATE','NO'),
               ('GENERATE_TEST_FILE','NO'),
               ('GENERATE_SIMULINK_INTERFACE','NO'),
               ('GENERATE_MAKE_FILE','NO'),
               ('CG_USE_C99','YES')]
               
    acadoOpts += intOpts

#    mhe.minimizeLsq(C.veccat([mhe['x'],mhe['u']]))
#    mhe.minimizeLsqEndTerm(C.veccat([mhe['x']]))
    mhe.minimizeLsq(mhe['measurements'])
    mhe.minimizeLsqEndTerm(mhe['measurementsN'])

    cgOpts = {'CXX':'g++', 'CC':'gcc'}
    mheRT = mhe.exportCode(codegenOptions=cgOpts,acadoOptions=acadoOpts)
    
    
    return mheRT, intOpts

#if __name__=='__main__':
#    from highwind_carousel_conf import conf
#    dae = rawe.models.carousel(conf)
#
#    OcpRt = makeMhe(dae,10,0.1)
