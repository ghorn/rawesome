import matplotlib.pyplot as plt
import numpy

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
    N = 100
    mpc = rawe.ocp.Ocp(dae, N=N, ts=0.2)

#    mpc.constrain(mpc['pos'], '==', 0, when='AT_START')
#    mpc.constrain(mpc['vel'], '==', 0, when='AT_START')

    mpc.constrain(mpc['pos'], '==', 0.38, when='AT_END')
    mpc.constrain(mpc['vel'], '==', 0, when='AT_END')

    mpc.constrain(mpc['vel'], '<=', 0.2)
    mpc.constrain(-0.1, '<=', mpc['force'], '<=', 0.1)

    mpc.minimizeLsq(C.veccat([mpc['pos'],mpc['vel'],mpc['force']]))
    mpc.minimizeLsqEndTerm(C.veccat([mpc['pos'], mpc['vel']]))

#    cgOptions = {'CXX':'clang++', 'CC':'clang'}
    cgOptions = {'CXX':'g++', 'CC':'gcc'}
#    cgOptions = {'CXX':'icpc', 'CC':'icc'}
    phase1Opts = {'CXX':'g++'}
    acadoOptions = [("HESSIAN_APPROXIMATION",     "GAUSS_NEWTON"),
                    ("DISCRETIZATION_TYPE",       "MULTIPLE_SHOOTING"),
                    ("INTEGRATOR_TYPE",           "INT_IRK_RIIA3"),
                    ("NUM_INTEGRATOR_STEPS",      str(N*5)),
                    ("LINEAR_ALGEBRA_SOLVER",     "GAUSS_LU"),
#                    ("SPARSE_QP_SOLUTION",        "FULL_CONDENSING"),
                    ("SPARSE_QP_SOLUTION",        "CONDENSING"),
                    ("QP_SOLVER",                 "QP_QPOASES"),
#                    ("SPARSE_QP_SOLUTION",        "SPARSE_SOLVER"),
#                    ("QP_SOLVER",                 "QP_QPDUNES"),
                    ("FIX_INITIAL_STATE",         "YES"),
                    ("HOTSTART_QP",               "YES"),
                    ("GENERATE_MAKE_FILE",        "NO")]
#                    ("GENERATE_MATLAB_INTERFACE", "YES")]
    ocpRt = mpc.exportCode(cgOptions=cgOptions,acadoOptions=acadoOptions,
                           phase1Options=phase1Opts)
    print '='*80

    # set the cost hessians
#    xRms = 0.2
#    vRms = 0.2
#    fRms = 20
    ocpRt.S[0,0] = 1.0#/(xRms*N)**2
    ocpRt.S[1,1] = 1.0#/(vRms*N)**2
    ocpRt.S[2,2] = 1.0#/(fRms*N)**2

    ocpRt.SN[0,0] = 1.0#/(xRms*N)**2
    ocpRt.SN[1,1] = 1.0#/(vRms*N)**2

    # make an initial guess
    for k in range(N):
        ocpRt.u[k] = 0.0
    ocpRt.initializeNodesByForwardSimulation()

#    import pickle
#    (ocpRt.x, ocpRt.u) = pickle.load( open('rt_initial_guess.dat','rb'))

    # plot initial guess
    plt.figure()
    plt.subplot(331)
    plt.plot(ocpRt.x[:,0])
    plt.title('x initial guess')
    plt.subplot(334)
    plt.plot(ocpRt.x[:,1])
    plt.title('v initial guess')
    plt.subplot(337)
    plt.plot(ocpRt.u[:,0])
    plt.title('u initial guess')

    # set tracking trajectory
    ocpRt.yN[0,0] = 0
    ocpRt.yN[1,0] = 0

    # iterate
    kkts = []
    objs = []
    for k in range(5):
        ocpRt.preparationStep()
        fbret = ocpRt.feedbackStep()
        if fbret != 0:
            raise Exception("feedbackStep returned error code "+str(fbret))
        print "sqp iteration",k,"\tkkts:",ocpRt.getKKT(),"\tobjective:",ocpRt.getObjective()
        kkts.append(ocpRt.getKKT() + 1e-200)
        objs.append(ocpRt.getObjective() + 1e-200)

#    import pickle
#    print "saving to rt_initial_guess.dat'
#    pickle.dump( (ocpRt.x, ocpRt.u), open('rt_initial_guess.dat','wb'))

    # plot results
    plt.subplot(332)
    plt.plot(ocpRt.x[:,0])
    plt.title('x after feedback')
    plt.subplot(335)
    plt.plot(ocpRt.x[:,1])
    plt.title('v after feedback')
    plt.subplot(338)
    plt.plot(ocpRt.u[:,0])
    plt.title('u after feedback')

    plt.subplot(233)
    plt.semilogy(kkts)
    plt.title('kkt convergence')
    plt.subplot(236)
    plt.semilogy(objs)
    plt.title('objective convergence')
    plt.show()
