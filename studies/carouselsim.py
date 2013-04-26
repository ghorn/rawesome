import casadi as C

import time

import rawe
import rawekite
from rawekite.carouselSteadyState import getSteadyState

if __name__=='__main__':
    # create the model
    from highwind_carousel_conf import conf
    dae = rawe.models.carousel(conf)
    
    # compute the steady state
    steadyState, ssDot = getSteadyState(dae,conf,2*C.pi,1.2,-0.1)
    print steadyState
    print ssDot
    
    
    # create the sim
    dt = 0.01
    sim = rawe.sim.Sim(dae,dt)
    communicator = rawekite.communicator.Communicator()
#    js = joy.Joy()

    # set the initial state from steadyState
    x = {}
    for name in dae.xNames():
        x[name] = steadyState[name]
    u = {}
    for name in dae.uNames():
        u[name] = steadyState[name]
    p = {}
    for name in dae.pNames():
        p[name] = steadyState[name]

    # set up the sim timer
    timer = rawe.sim.Timer(dt)
    timer.start()
    # loop through and simulate, if there's an error close the communicator and throw exception
    try:
        while True:
#        for k in range(100):
            # sleep for dt
            timer.sleep()
#            time.sleep(0.1)
            # send message to visualizer/plotter
            outs = sim.getOutputs(x,u,p)
            outs['delta'] = C.arctan2(x['sin_delta'], x['cos_delta'])
            communicator.sendKite(x,u,p,outs,conf)
            # try to take a simulation step of dt
            try:
                x = sim.step(x,u,p)
            except RuntimeError:
                # problem simulating, close the communicator
                communicator.close()
                raise Exception('OH NOES, IDAS CHOKED')
    except KeyboardInterrupt:
        print "closing..."
        communicator.close()
        pass
