import LLT_solver
import numpy

operAlphaDegLst = range(-5,36,1)
operRates=[0,0,0]
geomRoot=[1,0,0]
geomTip=[1,5,-0]
aeroCLaRoot = numpy.array([ -0.0, 0.0, 2*numpy.pi,-0.0000])
aeroCLaTip  = numpy.array([ -0.0, 0.0, 2*numpy.pi,-0.0000])
#aeroCLaRoot[:]=[ -28.2136, 16.4140, 0.9568,-0.4000]
#aeroCLaTip[:]=[ -28.2136, 16.4140, 0.9568, -0.4000]


LLT_solver.LLT_sovler(operAlphaDegLst, operRates, geomRoot, geomTip, aeroCLaRoot, aeroCLaTip)
