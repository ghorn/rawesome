import casadi as C
import casadi.tools as CT
import numpy
from LLT_solver import setupImplicitFunction
import geometry
import matplotlib.pyplot as plt
import os

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

# callback
class MyCallback:
    def __init__(self):
        self.iters = []
    def __call__(self,f,*args):
        self.iters.append(numpy.array(f.getInput("x")))
mycallback = MyCallback()
pyfun = C.PyFunction( mycallback, C.nlpSolverOut(x=C.sp_dense(dvs.size,1),
                                                 f=C.sp_dense(1,1),
                                                 lam_x=C.sp_dense(dvs.size,1),
                                                 lam_g = C.sp_dense(g.size,1),
                                                 lam_p = C.sp_dense(0,1),
                                                 g = C.sp_dense(g.size,1) ),
                      [C.sp_dense(1,1)] )
pyfun.init()

solver = C.IpoptSolver(nlp)
solver.setOption('tol',1e-11)
solver.setOption('linear_solver','ma27')
#solver.setOption('linear_solver','ma57')
solver.setOption("iteration_callback",pyfun)
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

lbg['sref']  = 9.0
ubg['sref']  = 9.0

lbg['CL']  = 1e-6
#lbg['CDi']  = 1e-6
#lbg['sref']  = 6.0
#ubg['sref']  = 6.0

lbx['aIncGeometricLoc'] = numpy.radians(-4)
ubx['aIncGeometricLoc'] = numpy.radians(4)
x0['aIncGeometricLoc'] = numpy.radians(-2)

lbx['operAlpha'] = numpy.radians(6)
#ubx['operAlpha'] = numpy.radians(6)

x0['An'] = 0.014
x0['chordLoc'] = 1.0

lbx['chordLoc'] = 0.3
ubx['chordLoc'] = 1.5
#lbx['chordLoc'] = 1.0
#ubx['chordLoc'] = 1.0
#print lbx['chordLoc']
#print ubx['chordLoc']

solver.solve()
ret = solver.getStat('return_status')
assert ret in ['Solve_Succeeded','Solved_To_Acceptable_Level'], 'Solver failed: '+ret

result = dvs(solver.output("x"))

def make_movie(draw=False):
    if draw:
        plt.ion()
    files = []
    fig = plt.figure(figsize=(10,8))
    fig.subplots_adjust(hspace=0.05,wspace=0.05)
    ax = fig.add_subplot(211)
    a2 = fig.add_subplot(212)
    ax.set_aspect('equal')
    for i,x_iter in enumerate(mycallback.iters):
        ax.cla()
        a2.cla()
        dvs_iter = dvs(x_iter)
        chordLoc = numpy.squeeze(numpy.array(dvs_iter['chordLoc']))
        aIncGeometricLoc = numpy.squeeze(numpy.array(dvs_iter['aIncGeometricLoc']*180.0/3.1415))
        ax.plot(numpy.append(numpy.flipud(yLoc),yLoc),
                numpy.append(-numpy.flipud(chordLoc)*0.5,chordLoc*0.5),
        )
        ax.set_ylim([-0.9,0.9])
        ax.set_xlim([0,5.5])
        a2.plot(numpy.flipud(yLoc),numpy.flipud(aIncGeometricLoc))
        a2.set_ylim([-5,4])
        a2.set_xlim([0,5.5])
        ax.legend(['wing shape'])
        a2.legend(['wing twist'])
        ax.set_xlabel('span [m]')
        a2.set_xlabel('span [m]')
        a2.set_ylabel('twist [deg]')
        fname = 'data/_tmp%03d.png'%i
        print 'Saving frame', fname
        if not draw:
            fig.savefig(fname)
            files.append(fname)
        if draw:
            plt.draw()
        #import time; time.sleep(1000)
    if not draw:
        print 'Making movie animation.mpg - this make take a while'
        os.system("mencoder 'mf://data/_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o data/animation.mpg")
        os.system('rm data/_tmp*.png')

#make_movie(draw=True); import sys; sys.exit()
#make_movie(draw=False); import sys; sys.exit()

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
