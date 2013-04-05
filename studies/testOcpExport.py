import rawe
import casadi as C

if __name__=='__main__':
    dae = rawe.dae.Dae()

    [pos,vel] = dae.addX( ["pos","vel"] )
    dae.addZ('dummyZ')
    force = dae.addU( "force" )
    endTime = dae.addP( 'someRandomParameter' )
    
    # specify the dae residual
    dae.setResidual([dae.ddt('pos') - vel,
                     dae.ddt('vel') - (force - 3.0*pos - 0.2*vel),
                     dae['dummyZ']])

    from rawe.ocp import Ocp
    mpc = Ocp(dae, N=10, ts=0.2)
    mpc.constrain(mpc['pos'], '==', 3, when='AT_START')
    mpc.constrain(mpc['vel'], '==', 0, when='AT_START')

    mpc.constrain(mpc['pos'], '==', 0, when='AT_END')
    mpc.constrain(mpc['vel'], '==', 0, when='AT_END')

    mpc.constrain(mpc['force']**2, '<=', 4)
    mpc.constrain(mpc['someRandomParameter']/4.0, '<=', 2)

    mpc.minimizeLsq(C.veccat([mpc['pos'],mpc['vel'],mpc['someRandomParameter']]))
    mpc.minimizeLsqEndTerm(C.veccat([mpc['pos']]))

#    options = {'CXX':'clang++', 'CC':'clang'}
    options = {'CXX':'g++', 'CC':'gcc'}
    mpcRt = mpc.exportCode(options,qpSolver='QP_OASES')

    mpcRt.initializeNodesByForwardSimulation()
    print "objective: "+str(mpcRt.getObjective())
    print "objective type: "+str(type(mpcRt.getObjective()))
