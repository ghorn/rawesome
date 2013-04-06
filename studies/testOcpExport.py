import rawe
import casadi as C

if __name__=='__main__':
    # specify the dae
    dae = rawe.dae.Dae()

    [pos,vel] = dae.addX( ["pos","vel"] )
    force = dae.addU( "force" )
    
    dae.setResidual([dae.ddt('pos') - vel,
                     dae.ddt('vel') - (force - 0*3.0*pos - 0.2*vel)])

    # specify the Ocp
    N = 20
    mpc = rawe.ocp.Ocp(dae, N=N, ts=0.2)
    mpc.constrain(mpc['pos'], '==', 1, when='AT_START')
    mpc.constrain(mpc['vel'], '==', 0, when='AT_START')

    mpc.constrain(mpc['pos'], '==', 0, when='AT_END')
    mpc.constrain(mpc['vel'], '==', 0, when='AT_END')

    mpc.constrain(-5, '<=', mpc['vel'], '<=', 6)
    mpc.constrain(mpc['force']**2, '<=', 4)

    mpc.minimizeLsq(C.veccat([mpc['pos'],mpc['vel'],mpc['force']]))
    mpc.minimizeLsqEndTerm(C.veccat([mpc['pos']]))

    cgOptions = {'CXX':'clang++', 'CC':'clang'}
#    cg_options = {'CXX':'g++', 'CC':'gcc'}
    acadoOptions = [("HESSIAN_APPROXIMATION",     "GAUSS_NEWTON"),
                    ("DISCRETIZATION_TYPE",       "MULTIPLE_SHOOTING"),
                    ("INTEGRATOR_TYPE",           "INT_IRK_RIIA3"),
                    ("NUM_INTEGRATOR_STEPS",      "N * 5"),
                    ("LINEAR_ALGEBRA_SOLVER",     "GAUSS_LU"),
                    ("SPARSE_QP_SOLUTION",        "FULL_CONDENSING"),
                    ("QP_SOLVER",                 "QP_QPOASES"),
                    ("HOTSTART_QP",               "YES"),
                    ("GENERATE_MATLAB_INTERFACE", "YES")]
    ocpRt = mpc.exportCode(cgOptions=cgOptions,acadoOptions=acadoOptions,qpSolver='QP_OASES')


    # run the generated code
    print '='*80
    for k in range(ocpRt.u.size):
        ocpRt.u[k] = 0.1*(k+1)
    ocpRt.S[0,0] = 5;
    ocpRt.S[1,1] = 0.1;
    ocpRt.S[2,2] = 3;
    for k in range(ocpRt.SN.shape[0]):
        ocpRt.SN[k,k] = 10.0

    print "u:"
    print ocpRt.u
    ocpRt.initializeNodesByForwardSimulation()
    print "x:"
    print ocpRt.x
    print "x0:"
    print ocpRt.x0

    for k in range(ocpRt.S.shape[0]):
        ocpRt.S[k,k] = 0.01

    for k in range(100):
        ocpRt.preparationStep()
        prepret = ocpRt.feedbackStep()
        if prepret != 0:
            raise Exception("feedbackStep returned error code "+str(prepret))
    print "x after feedback:"
    print ocpRt.x
    print "u after feedback:"
    print ocpRt.u
    print "x0 after feedback:"
    print ocpRt.x0

    print "kkt: "+str(ocpRt.getKKT())
    print "kkt type: "+str(type(ocpRt.getKKT()))
    print "objective: "+str(ocpRt._lib.getObjective())
    print "objective type: "+str(type(ocpRt._lib.getObjective()))
