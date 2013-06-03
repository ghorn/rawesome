#This solver is meant to provide a lifting line theory based aerodynamic prediction for a supplied alpha range and gerometry.
import pylab
import numpy
import casadi as C

def setupImplicitFunction(operAlpha, An, geom):

    #setup function to obtain B and RHS for the Fourier series calc
    def getBandRHS(alphaiLoc):
        alphaLoc = geom.alphaGeometric(operAlpha) - alphaiLoc
        clLoc = geom.clLoc(alphaLoc)
        RHS = clLoc*geom.chordLoc
        B = numpy.sin(numpy.outer(geom.thetaLoc, geom.sumN))*(4.0*geom.bref)
        return (B,RHS)

    alphaiLoc = []
    for i in range(geom.n):
        alphaiLoc.append(
            C.inner_prod(geom.sumN*numpy.sin(geom.sumN*geom.thetaLoc[i])/numpy.sin(geom.thetaLoc[i]), An))
    alphaiLoc = C.veccat(alphaiLoc)

    (B,RHS) = getBandRHS(alphaiLoc)
    makeMeZero = RHS - C.mul(B, An)

    CL  = An[0]*numpy.pi*geom.AR
    CDi = 0.0
    for i in range(geom.n):
        k = geom.sumN[i]
        CDi += k * An[i]**2 / An[0]**2
    CDi *= numpy.pi*geom.AR*An[0]**2

    return (makeMeZero, alphaiLoc, CL, CDi)


def LLT_solver(operAlphaDegLst, geom):
    operAlphaLst = [numpy.radians(x) for x in operAlphaDegLst]
    operCLLst   = []
    operCDiLst  = []
    #
    #Now let's iterate over the Alphas and solve LLT equations to obtain lift for each one
    #
    alpha = C.ssym('alpha')
    An = C.ssym('A',geom.n)
    (makeMeZero, alphaiLoc, _, _) = setupImplicitFunction(alpha, An, geom)
    
    # make the solver
    f = C.SXFunction([An,alpha], [makeMeZero])
    f.init()
    
    solver = C.NLPImplicitSolver(f)
    solver.setOption('nlp_solver', C.IpoptSolver)
    #solver.setOption('nlp_solver_options',{'suppress_all_output':'yes','print_time':False})
    solver.init()

    for operAlpha in operAlphaLst:
        print 'Current Alpha = ', numpy.degrees(operAlpha)
        guess = numpy.zeros(geom.n)
        guess[0] = 0.01
        solver.setOutput(guess,0)
        solver.setInput(operAlpha)
        solver.evaluate()
        iterAn = numpy.squeeze(numpy.array(solver.output()))

        #Compute CL/CDi
        CL  = iterAn[0]*numpy.pi*geom.AR
        CD0 = 0
        if iterAn[0] != 0:
            for i in range(geom.n):
                k = geom.sumN[i]
                CD0 += k*iterAn[i]**2/(iterAn[0]**2)
        CDi = numpy.pi*geom.AR*iterAn[0]**2*CD0

        operCLLst.append(CL)
        operCDiLst.append(CDi)
    operCLLst = numpy.array(operCLLst)
    operCDiLst = numpy.array(operCDiLst)

    #
    # plot everything
    #
    pylab.plot(numpy.degrees(operAlphaLst),operCLLst,'r',numpy.degrees(operAlphaLst), numpy.array(operAlphaLst)*2*numpy.pi*10.0/12.0)
    pylab.xlabel('Alpha')
    pylab.ylabel('CL')
    pylab.legend(['llt CL','flat plate CL'])

    pylab.figure()
    pylab.plot(numpy.degrees(operAlphaLst),operCDiLst,'r',numpy.degrees(operAlphaLst),operCLLst[:]**2/(numpy.pi*geom.AR),'b')
    pylab.legend(['llt CDi','flat plate CDi'])
    pylab.xlabel('Alpha')
    pylab.ylabel('CDi')

    pylab.figure()
    pylab.plot(operCDiLst,operCLLst)
    pylab.xlabel('CDi')
    pylab.ylabel('CL')
    pylab.show()
