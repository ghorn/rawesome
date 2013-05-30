import LLT_solver
import numpy
operAseq=[-5,35,1]
operRates=[0,0,0]
geomRoot=[1,0,0]
geomTip=[1,5,-0]
aeroCLaRoot=numpy.zeros(4)
aeroCLaTip=numpy.zeros(4)
aeroCLaRoot[:]=[ -0.0, 0.0, 2*numpy.pi,-0.0000]
aeroCLaTip[:]=[ -0.0, 0.0, 2*numpy.pi,-0.0000]
#aeroCLaRoot[:]=[ -28.2136, 16.4140, 0.9568,-0.4000]
#aeroCLaTip[:]=[ -28.2136, 16.4140, 0.9568, -0.4000]
LLT_solver.LLT_sovler(operAseq, operRates, geomRoot, geomTip, aeroCLaRoot, aeroCLaTip)