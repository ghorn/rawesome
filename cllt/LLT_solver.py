#This solver is meant to provide a lifting line theory based aerodynamic prediction for a supplied alpha range and gerometry.
import pylab
import numpy
import scipy.linalg

def woohoo(operAlpha, aIncFromZeroLift, aIncGeometric, lltCLa, lltTheta, lltChord, n, geomAR, geomBref):
    alphaFromZeroLift = operAlpha + aIncFromZeroLift
    alphaGeometric    = operAlpha + aIncGeometric

    #
    #setup function to obtain B and RHS for the Fourier series calc
    #
    def getBandRHS(iterAlphaiLoc):
        iterAlphaLoc    = alphaFromZeroLift - iterAlphaiLoc
        lltAlphaLoc     = alphaGeometric    - iterAlphaiLoc
        iterCLALoc   = (lltCLa[0,:]*lltAlphaLoc[:]**3 + lltCLa[1,:]*lltAlphaLoc[:]**2 + lltCLa[2,:]*lltAlphaLoc[:] + lltCLa[3,:])/(iterAlphaLoc[:])
        RHS             = iterAlphaLoc*numpy.sin(lltTheta)*lltChord*iterCLALoc/(4.*geomBref)
        B               = numpy.zeros((n,n))
        mu              = lltChord * iterCLALoc/(4.*geomBref)
        for i in range(n):
            l       = float(2*i+1)
            B[:,i]  = numpy.sin(l*lltTheta)*(mu*l+numpy.sin(lltTheta))
        return (B,RHS)

    iterAlphaiLoc   = (lltCLa[0,:]*alphaGeometric[:]**3 + lltCLa[1,:]*alphaGeometric[:]**2 + lltCLa[2,:]*alphaGeometric[:] + lltCLa[3,:])/(numpy.pi*geomAR)

    #
    #iterate to solve for Fourier Coefficients
    #
    iterConvCrit    = 1e-6
    iterResAlphai   = 10
    while iterResAlphai > iterConvCrit:
        iterAlphaiLocOld = numpy.copy(iterAlphaiLoc)
        (iterB, iterRHS)    = getBandRHS(iterAlphaiLoc)
        iterAn              = scipy.linalg.solve(iterB, iterRHS)
        #
        #calculate induced alphas
        #
        sumN = numpy.linspace(1,2*n+1,num=n,endpoint=True)
        for i in range(n):
            iterAlphaiLoc[i] = sum(sumN*iterAn*numpy.sin(sumN*lltTheta[i])/numpy.sin(lltTheta[i]))
        #
        #update local alphas and calculate residuals
        #
        iterResAlphai   = numpy.max(numpy.absolute(iterAlphaiLocOld - iterAlphaiLoc))
        print 'Residual = ', iterResAlphai
    #
    #Compute CL
    #
    CL  = iterAn[0]*numpy.pi*geomAR
    #print CL, iterAn[0]
    #
    #Compute CDi
    #
    CD0 = 1
    if iterAn[0] != 0:
        for i in range(1,n):
            CD0 = CD0 + (2.*i+1)*iterAn[i]**2/(iterAn[0]**2)
    CDi     = numpy.pi*geomAR*iterAn[0]**2*CD0

    return (CL, CDi, iterAn)


def LLT_sovler(operAlphaDegLst, operRates, geomRoot, geomTip, aeroCLaRoot, aeroCLaTip):
    #setup defining variables from inputs
    #
    #geomRoot is a 1x3 vector containing chord, span location and angle of incidence
    #
    geomRootChord   = geomRoot[0]
    geomRootY       = geomRoot[1]
    geomRootAinc    = numpy.radians(geomRoot[2])
    geomTipChord    = geomTip[0]
    geomTipY        = geomTip[1]
    geomTipAinc     = numpy.radians(geomTip[2])
    geomBref        = 2*(geomTipY-geomRootY) #reference span for A/C
    geomCref        = 0.5*(geomRootChord+geomTipChord)  #reference chord for A/c
    geomSref        = (geomBref*geomTipChord)+(0.5*geomBref*abs(geomRootChord-geomTipChord)) #reference wing surface area (projected)
    geomAR          = (geomBref**2)/geomSref  #Aspect ratio Bref^2/Sref
    #
    #aeroCLaRoot and aeroCLaTip are 1x4 vectors, containing coefficients for a cubic polynomial describing a curve representing
    #the function CL(Alpha) for a 2D airfoil. Once this polynomial is setup, we perform an analytical differentiation to obtain
    #curves for CL_alpha(Alpha), ie, the lift slope as a function of Alpha
    #
    aeroCLAaRoot = numpy.array([3*aeroCLaRoot[0], 2*aeroCLaRoot[1], 1*aeroCLaRoot[2]])
    aeroCLAaTip  = numpy.array([3*aeroCLaTip[0], 2*aeroCLaTip[1], 1*aeroCLaTip[2]])
    #
    #Now the minimum absolute root of the polynomial for CLa provides us the zero-lift andgles
    #
    aeroAlphaCLoRoot    = min(numpy.absolute(numpy.roots(aeroCLaRoot)))
    aeroAlphaCLoTip     = min(numpy.absolute(numpy.roots(aeroCLaTip)))
    print "Zero lift alpha at the root =", numpy.degrees(aeroAlphaCLoRoot)
    print "Zero lift alpha at the tip =", numpy.degrees(aeroAlphaCLoTip)
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
    operAlphaLst = [numpy.radians(x)+1e-18 for x in operAlphaDegLst]
    #################Alpha hard set option
    #operAlphaLst=[0]
    ##########################################
    print "Alphas= ", operAlphaLst, numpy.degrees(operAlphaLst)
    #
    #Before we start iterating, let's setup lists for storing lift and induced drag
    #
    operCLLst   = []
    operCDiLst  = []
    #
    #Next, lets choose the number of spanwise stations (and hence number of Fourier series terms) that we will divide each half wing into
    #We will also pick a suitable convergence criterion for the non-linear iteration and a relaxation factor to aid convergence
    #
    n = 20
    #
    #Finally, let's create n span stations at which LLT equations will be solved and define some fixed wing properties at each station
    #
    lltTheta    = numpy.linspace(1,n,num=n,endpoint=True)*numpy.pi/(2.*n)
    lltY        = (0.5*geomBref)*numpy.cos(lltTheta)
    lltChord    = geomRootChord + 2.*(geomTipChord - geomRootChord)*lltY/geomBref
    iterAinc    = geomRootAinc - (aeroAlphaCLoRoot + 2.*(aeroAlphaCLoTip - aeroAlphaCLoRoot + geomRootAinc - geomTipAinc)*lltY/geomBref)
    lltAinc     = geomRootAinc - (2.*(geomRootAinc - geomTipAinc)*lltY/geomBref)
    lltCLa      = numpy.zeros((4,n))    #setup lift slope curves for each station
    for i in range(n):
        lltCLa[:,i]    = aeroCLaRoot[:] + 2.*(aeroCLaTip[:] - aeroCLaRoot[:])*lltY[i]/geomBref
    #
    #Now let's iterate over the Alphas and solve LLT equations to obtain lift for each one
    #
    for operAlpha in operAlphaLst:
        print 'Current Alpha = ', numpy.degrees(operAlpha)
        (CL, CDi, iterAn) = woohoo(operAlpha, iterAinc, lltAinc, lltCLa, lltTheta, lltChord, n, geomAR, geomBref)
        operCLLst.append(CL)
        operCDiLst.append(CDi)
    operCLLst = numpy.array(operCLLst)
    operCDiLst = numpy.array(operCDiLst)
    #
    #plot polar
    #
    pylab.plot(numpy.degrees(operAlphaLst),operCLLst)
    pylab.xlabel('Alpha')
    pylab.ylabel('CL')
    pylab.show()
    pylab.plot(numpy.degrees(operAlphaLst),operCDiLst,'r',numpy.degrees(operAlphaLst),operCLLst[:]**2/(numpy.pi*geomAR),'b')
    pylab.legend(['llt CDi','flat plate CDi'])
    pylab.xlabel('Alpha')
    pylab.ylabel('CDi')
    pylab.show()
    pylab.plot(operCDiLst,operCLLst)
    pylab.xlabel('CDi')
    pylab.ylabel('CL')
    pylab.show()
