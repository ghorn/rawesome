import casadi as C
from collocation import Coll,LagrangePoly
import models

from newton.newton2 import Newton

if __name__ == '__main__':
    print "creating model"
    dae = models.pendulum2()
#    dae.addP('endTime')

    nk = 30
    nicp = 1
    deg = 4
    newton = Newton(LagrangePoly,dae,nk,nicp,deg,'RADAU')
    endTime = 0.01*nk
    newton.setupStuff(endTime)

    r = 0.3
    x0 = [r,0,0,0]

#    newton.isolver.input(0).set(x0)
    newton.isolver.setInput(x0,0)
    newton.isolver.setInput(0,1) # torque
    newton.isolver.setInput(0.3,2) # mass
    newton.isolver.evaluate()
    newton.isolver.evaluate()

    j = newton.isolver.jacobian(0,1)
    j.init()
    j.setInput(x0,0)
    j.setInput(0,1)
    j.setInput(0.3,2)
    j.evaluate()
    print j.output(0)
    print j.output(1)
