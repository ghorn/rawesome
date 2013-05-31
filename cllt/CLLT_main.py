import LLT_solver
from geometry import Geometry
import numpy

operAlphaDegLst = range(-5,15,1)
operRates = [0,0,0]
geomRoot = [1,0,0]
geomTip = [1,5,0]
aeroCLaRoot = numpy.array([0.0, 0.0, 2*numpy.pi, 0.0])
aeroCLaTip  = numpy.array([0.0, 0.0, 2*numpy.pi, 0.0])
#aeroCLaRoot = numpy.array([ -28.2136, 16.4140, 0.9568,-0.4000])
#aeroCLaTip = numpy.array([ -28.2136, 16.4140, 0.9568, -0.4000])

n = 30

geom = Geometry(geomRoot, geomTip, aeroCLaRoot, aeroCLaTip, n)

LLT_solver.LLT_sovler(operAlphaDegLst, operRates, geom)
