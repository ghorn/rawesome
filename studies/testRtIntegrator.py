import rawe
import casadi as C

if __name__=='__main__':
    dae = rawe.dae.Dae()

    [pos,vel] = dae.addX( ["pos","vel"] )
    force = dae.addU( "force" )
    endTime = dae.addP( 'endTime' )

    # specify the dae residual
    dae.setResidual([dae.ddt('pos') - vel,
                     dae.ddt('vel') - (force - 3.0*pos - 0.2*vel)])

    from rawe import RtIntegrator, RtIntegratorOptions
    intOpts = RtIntegratorOptions()
    intOpts['INTEGRATOR_TYPE'] = 'INT_IRK_RIIA3'
    intOpts['NUM_INTEGRATOR_STEPS'] = 5
    intOpts['LINEAR_ALGEBRA_SOLVER'] = 'GAUSS_LU'
    integrator = RtIntegrator(dae,ts=endTime, options=intOpts,
                              measurements=C.veccat([dae.ddt('pos'), vel, force]))
    
    x = {'pos':5.3, 'vel':0.6}
    u = {'force':-4.2}
    p = {'endTime':0.2}

    #xdot = {'pos':9.4,'vel':1.2}
    #print integrator.rhs(xdot,x,{},u,{}, compareWithSX=True)
    #print integrator.rhsJac(xdot,x,{},u,{}, compareWithSX=True)

    log = {'pos':[x['pos']], 'vel':[x['vel']]}
    print "running integrator..."
    for k in range(400):
        x = integrator.step(x,u,p)
        log['pos'].append(x['pos'])
        log['vel'].append(x['vel'])

    import matplotlib.pylab as plt
    plt.figure()
    plt.subplot(211)
    plt.plot(log['pos'])
    plt.title('position')
    plt.subplot(212)
    plt.plot(log['vel'])
    plt.title('velocity')
    plt.show()
