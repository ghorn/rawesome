import casadi as C

from carouselSteadyState import getSteadyState
import rawe

if __name__=='__main__':
    print "creating model"
    from highwind_carousel_conf import conf
    dae = rawe.models.carousel(conf)
    steadyState = getSteadyState(dae,conf,2*C.pi,1.2)

    dt = 0.02
    sim = rawe.sim.Sim(dae,dt)
    communicator = rawe.simutils.Communicator()
#    js = joy.Joy()
    x = {}
    for name in dae.xNames():
        x[name] = steadyState[name]
    u = {}
    for name in dae.uNames():
        u[name] = steadyState[name]
    p = {}
    for name in dae.pNames():
        p[name] = steadyState[name]

    print "simulating..."
    timer = rawe.simutils.Timer(dt)
    timer.start()
    try:
        while True:
            timer.sleep()
            outs = sim.getOutputs(x,u,p)
            communicator.sendKite(x,u,p,outs,conf)
            try:
                x = sim.step(x,u,p)
            except RuntimeError:
                communicator.close()
                raise Exception('OH NOES, IDAS CHOKED')
    except KeyboardInterrupt:
        print "closing..."
        communicator.close()
        pass
