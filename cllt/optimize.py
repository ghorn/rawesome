import casadi as C
import casadi.tools as CT
import numpy
from LLT_solver import setupImplicitFunction
import geometry
import matplotlib.pyplot as plt

# geomRootChord = geomRoot[0]
# geomRootY     = geomRoot[1]
# geomRootAinc  = numpy.radians(geomRoot[2])
# geomTipChord  = geomTip[0]
# geomTipY      = geomTip[1]
# geomTipAinc   = numpy.radians(geomTip[2])
n = 50

dvs = CT.struct_ssym(['operAlpha',
                      CT.entry('chordLoc',shape=n),
                      CT.entry('aIncGeometricLoc',shape=n),
                      CT.entry('An',shape=n)])

thetaLoc = numpy.linspace(1,n,n)*numpy.pi/(2.0*n)
bref = 10.0
yLoc = 0.5*bref*numpy.cos(thetaLoc)
chordLoc = dvs['chordLoc']
aIncGeometricLoc = dvs['aIncGeometricLoc']

aeroCLaRoot = numpy.array([0.0, 0.0, 2*numpy.pi, 0.0])
aeroCLaTip  = numpy.array([0.0, 0.0, 2*numpy.pi, 0.0])
#aeroCLaRoot = numpy.array([ -28.2136, 16.4140, 0.9568,-0.4000])
#aeroCLaTip = numpy.array([ -28.2136, 16.4140, 0.9568, -0.4000])
clPolyLoc = numpy.zeros((4,n))    #setup lift slope curves for each station
for i in range(n):
    clPolyLoc[:,i] = aeroCLaRoot[:] + 2.0*(aeroCLaTip[:] - aeroCLaRoot[:])*yLoc[i]/bref

geom = geometry.Geometry(thetaLoc, chordLoc, yLoc, aIncGeometricLoc, clPolyLoc, bref, n)

(makeMeZero, alphaiLoc, CL, CDi) = setupImplicitFunction(dvs['operAlpha'], dvs['An'], geom)
outputsFcn = C.SXFunction([dvs.cat], [alphaiLoc, CL, CDi, geom.sref])
outputsFcn.init()

g = CT.struct_SX([ CT.entry("makeMeZero",expr=makeMeZero),
                   CT.entry('alphaiLoc',expr=alphaiLoc),
                   CT.entry("CL",expr=CL),
                   CT.entry("CDi",expr=CDi),
                   CT.entry("sref",expr=geom.sref)])

obj = -CL / (CDi + 0.01)

nlp = C.SXFunction(C.nlpIn(x=dvs),C.nlpOut(f=obj,g=g))
solver = C.IpoptSolver(nlp)
solver.setOption('tol',1e-11)
solver.setOption('linear_solver','ma27')
#solver.setOption('linear_solver','ma57')
solver.init()

lbg = g(solver.input("lbg"))
ubg = g(solver.input("ubg"))
lbx = dvs(solver.input("lbx"))
ubx = dvs(solver.input("ubx"))
x0 = dvs(solver.input("x0"))

lbg["makeMeZero"] = 0.0
ubg["makeMeZero"] = 0.0

#lbg['alphaiLoc']  = numpy.radians(-10)
#ubg['alphaiLoc']  = numpy.radians( 15)

lbg['sref']  = 6.0
#ubg['sref']  = 6.0

lbg['CL']  = 1e-6
#lbg['CDi']  = 1e-6
#lbg['sref']  = 6.0
#ubg['sref']  = 6.0

lbx['aIncGeometricLoc'] = numpy.radians(-4)
ubx['aIncGeometricLoc'] = numpy.radians(0)

lbx['operAlpha'] = numpy.radians(6)
#ubx['operAlpha'] = numpy.radians(6)

x0['An'] = 0.014
x0['chordLoc'] = 1.0

lbx['chordLoc'] = 0.45
ubx['chordLoc'] = 1.0
#lbx['chordLoc'] = 1.0
#ubx['chordLoc'] = 1.0
#print lbx['chordLoc']
#print ubx['chordLoc']

solver.solve()
ret = solver.getStat('return_status')
assert ret in ['Solve_Succeeded','Solved_To_Acceptable_Level'], 'Solver failed: '+ret

result = dvs(solver.output("x"))
plt.figure()
plt.plot(yLoc, result['chordLoc'])
plt.xlabel('span [m]')
plt.ylabel('chord [m]')
plt.legend(['chord'])
plt.axes().set_aspect('equal')

plt.figure()
plt.plot(yLoc, result['aIncGeometricLoc']*180.0/numpy.pi)
plt.legend(['twist (deg)'])
plt.xlabel('span [m]')
plt.ylabel('twist [deg]')

outputsFcn.setInput(result.cat)
outputsFcn.evaluate()
alphaiLoc = outputsFcn.output(0)
CL = outputsFcn.output(1)
CDi = outputsFcn.output(2)
sref = outputsFcn.output(3)

plt.figure()
plt.plot(yLoc, alphaiLoc*180.0/numpy.pi)
plt.legend(['alphaiLoc(deg)'])
plt.xlabel('span [m]')
plt.ylabel('induced AOA [deg]')

print "oper alpha: "+str(result['operAlpha']*180/C.pi)+" degrees"
print "sref:",sref
print "CL:",CL
print "CDi:",CDi
print "CL/CDi:",CL/CDi
plt.show()
