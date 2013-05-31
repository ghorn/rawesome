import LLT_solver
from geometry import Geometry
import numpy
import casadi as C
import pylab

operAlphaDegLst = list(numpy.linspace(5,40,100))
geomRoot = [1,0,0]
geomTip = [0.4,5,0]
aeroCLaRoot = numpy.array([0.0, 0.0, 2*numpy.pi, 0.0])
aeroCLaTip  = numpy.array([0.0, 0.0, 2*numpy.pi, 0.0])
#aeroCLaRoot = numpy.array([ -28.2136, 16.4140, 0.9568,-0.4000])
#aeroCLaTip = numpy.array([ -28.2136, 16.4140, 0.9568, -0.4000])

n = 30
geom = Geometry(geomRoot, geomTip, aeroCLaRoot, aeroCLaTip, n)
LLT_solver.LLT_solver(operAlphaDegLst, geom)


# set up the LLT symbolics
alpha = C.ssym('alpha')
An = C.ssym('A',geom.n)
(makeMeZero, alphaiLoc) = LLT_solver.setupImplicitFunction(alpha, An, geom)

# make the solver
f = C.SXFunction([An,alpha], [makeMeZero])
f.init()

solver = C.NLPImplicitSolver(f)
solver.setOption('nlp_solver', C.IpoptSolver)
#solver.setOption('nlp_solver_options',{'suppress_all_output':'yes','print_time':False})
solver.init()

fAlphai = C.SXFunction([An], [alphaiLoc])
fAlphai.init()

################ now lets make a run woo ##############
operAlpha = numpy.radians(5)
solver.setInput(operAlpha)
guessAn = numpy.zeros(n)
guessAn[0] = 0.014
solver.setOutput(guessAn)
solver.evaluate()
An = numpy.squeeze(numpy.array(solver.output()))

fAlphai.setInput(An)
fAlphai.evaluate()
alphaiLoc = numpy.squeeze(numpy.array(fAlphai.output()))
alphaiLocDeg = numpy.degrees(alphaiLoc)

# plot alpha
pylab.plot(geom.yLoc, numpy.degrees(operAlpha)-alphaiLocDeg)
pylab.legend(['alpha(deg)'])

#for name in ['chordLoc','aIncGeometricLoc','thetaLoc']:
#    pylab.figure()
#    pylab.plot(geom.yLoc, getattr(geom,name))
#    pylab.legend([name])


# plot Cl
pylab.figure()
Cl = geom.clLoc(geom.alphaGeometric(operAlpha) - alphaiLoc)
pylab.plot(geom.yLoc, Cl)
Clellipse = Cl[-1]*numpy.sqrt(1 - geom.yLoc**2/25)
pylab.plot(geom.yLoc, Clellipse/geom.chordLoc)
pylab.legend(['Cl','Cl ellipse'])

## plot Gamma
#Gamma = numpy.dot(numpy.sin(numpy.outer(geom.thetaLoc, geom.sumN)), An)
#pylab.figure()
#pylab.plot(geom.thetaLoc, Gamma)
#pylab.legend(['gamma'])
#pylab.xlabel('theta')
#
#pylab.figure()
#pylab.plot(geom.yLoc, Gamma)
#pylab.legend(['gamma'])
#pylab.xlabel('y')
#
#Gamma = numpy.sin(numpy.outer(geom.thetaLoc, geom.sumN*An))
#pylab.figure()
#pylab.plot(geom.yLoc, Gamma)
#pylab.legend(['A'+str(k) for k in geom.sumN])
#pylab.xlabel('y')


# plot An
pylab.figure()
pylab.semilogy(solver.output(),'.')
pylab.legend(['An'])

pylab.show()
