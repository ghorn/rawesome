import casadi as C
from carouselSteadyState import getSteadyState
from config import readConfig
import models
from sim import Sim
import simutils
import joy

if __name__=='__main__':
    print "creating model"
    conf = readConfig('config.ini','configspec.ini')
    dae = models.carousel(conf)
    steadyState = getSteadyState(dae,conf,2*C.pi,1.2)

    dt = 0.02
    sim = Sim(dae,dt)
    communicator = simutils.Communicator()
    js = joy.Joy()
    x = steadyState
    u = {'daileron':0,'delevator':0,'motor torque':steadyState['motor torque'],'ddr':0}
    p = {'w0':0}

    print "simulating..."
    timer = simutils.Timer(dt)
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
