#This solver is meant to provide a lifting line theory based aerodynamic prediction for a supplied alpha range and gerometry.
import pylab
import numpy
import casadi as C

def setupImplicitFunction(geom):
    operAlpha      = C.ssym('alpha')
    alphaGeometric = operAlpha + geom.aIncGeometricLoc

    n = geom.n

    #setup function to obtain B and RHS for the Fourier series calc
    def getBandRHS(iterAlphaiLoc):
        alphaLoc = alphaGeometric - iterAlphaiLoc
        clLoc = geom.clPolyLoc[0,:]*alphaLoc**3 + \
                geom.clPolyLoc[1,:]*alphaLoc**2 + \
                geom.clPolyLoc[2,:]*alphaLoc + \
                geom.clPolyLoc[3,:]
        RHS = clLoc*geom.chordLoc/(4.0*geom.bref)
        B = numpy.sin(numpy.outer(geom.thetaLoc, geom.sumN))
        return (B,RHS)

    iterAn = C.ssym('iterAn',geom.n)
    iterAlphaiLoc = []
    for i in range(n):
        iterAlphaiLoc.append(
            C.inner_prod(geom.sumN*numpy.sin(geom.sumN*geom.thetaLoc[i])/numpy.sin(geom.thetaLoc[i]), iterAn))
    iterAlphaiLoc = C.veccat(iterAlphaiLoc)

    (B,RHS) = getBandRHS(iterAlphaiLoc)
    makeMeZero = RHS - C.mul(B, iterAn)
#    makeMeZero = C.solve(B,RHS) - iterAn

#    fIterAlphaiLoc = C.SXFunction([iterAn],[iterAlphaiLoc])
#    fIterAlphaiLoc.init()
#    fIterAlphaiLoc.setInput(iterAn)
#    fIterAlphaiLoc.evaluate()
#    iterAlphaiLoc = numpy.array(fIterAlphaiLoc.output())

    print "setting up function"
    f = C.SXFunction([iterAn,operAlpha], [makeMeZero])
    print "initializing function"
    f.init()

    print "setting up solver"
    i = C.NLPImplicitSolver(f)
    i.setOption('nlp_solver', C.IpoptSolver)
    print "initializing solver"
    i.setOption('nlp_solver_options',{'suppress_all_output':'yes',
                                      'print_time':False})
    i.init()

    print "done"
    return i


def LLT_sovler(operAlphaDegLst, operRates, geom):
    #
    #operAlpha is a 1x3 vector containing initial alpha, final alpha and alpha increment for a sweep.
    #
#    operAseqMin     = numpy.radians(operAseq[0])
#    operAseqMax     = numpy.radians(operAseq[1])
#    operAseqIncr    = numpy.radians(operAseq[2])
    #
    #now let's setup the range of alphas we are going to run this calc over
#    operAlphaLst = []
#    tempAlpha = operAseqMin
#    while tempAlpha < (operAseqMax + operAseqIncr / 10):
#        operAlphaLst.append(tempAlpha)
#        tempAlpha = tempAlpha + operAseqIncr
    operAlphaLst = [numpy.radians(x) for x in operAlphaDegLst]
    #################Alpha hard set option
    #operAlphaLst=[0]
    ##########################################
    #
    #Before we start iterating, let's setup lists for storing lift and induced drag
    #
    operCLLst   = []
    operCDiLst  = []
    #
    #Next, lets choose the number of spanwise stations (and hence number of Fourier series terms) that we will divide each half wing into
    #We will also pick a suitable convergence criterion for the non-linear iteration and a relaxation factor to aid convergence
    #
    #Finally, let's create n span stations at which LLT equations will be solved and define some fixed wing properties at each station
    #
    #
    #Now let's iterate over the Alphas and solve LLT equations to obtain lift for each one
    #
    f = setupImplicitFunction(geom)
    for operAlpha in operAlphaLst:
        print 'Current Alpha = ', numpy.degrees(operAlpha)
        f.setOutput(numpy.zeros(geom.n))
        f.setInput(operAlpha)
        f.evaluate()
        iterAn = numpy.squeeze(numpy.array(f.output()))

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
    pylab.show()
    pylab.plot(numpy.degrees(operAlphaLst),operCDiLst,'r',numpy.degrees(operAlphaLst),operCLLst[:]**2/(numpy.pi*geom.AR),'b')
    pylab.legend(['llt CDi','flat plate CDi'])
    pylab.xlabel('Alpha')
    pylab.ylabel('CDi')
    pylab.show()
    pylab.plot(operCDiLst,operCLLst)
    pylab.xlabel('CDi')
    pylab.ylabel('CL')
    pylab.show()
