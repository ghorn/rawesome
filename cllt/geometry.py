import numpy
import casadi as C

class Geometry(object):
    @property
    def aIncGeometricLoc(self):
        return self._aIncGeometricLoc

    @property
    def thetaLoc(self):
        return self._thetaLoc

    @property
    def chordLoc(self):
        return self._chordLoc

    @property
    def AR(self):
        return self._AR

    @property
    def bref(self):
        return self._bref

    @property
    def n(self):
        return self._n

    @property
    def sumN(self):
        return self._sumN

    @property
    def clPolyLoc(self):
        return self._clPolyLoc

    @property
    def yLoc(self):
        return self._yLoc

    def __init__(self, geomRoot, geomTip, aeroCLaRoot, aeroCLaTip, n):
        '''
        n is the number of spanwise stations (and hence number of Fourier series terms)
        that we will divide each half wing into
        '''
        #setup defining variables from inputs
        #
        #geomRoot is a 1x3 vector containing chord, span location and angle of incidence
        #
        self._n = n
        self._sumN = numpy.linspace(1, 2*self.n-1, n)

        geomRootChord = geomRoot[0]
        geomRootY     = geomRoot[1]
        geomRootAinc  = numpy.radians(geomRoot[2])
        geomTipChord  = geomTip[0]
        geomTipY      = geomTip[1]
        geomTipAinc   = numpy.radians(geomTip[2])

        #reference span for A/C
        self._bref = 2*(geomTipY-geomRootY)
        #reference chord for A/c
        self._cref = 0.5*(geomRootChord+geomTipChord)
        #reference wing surface area (projected)
        geomSref = (self.bref*geomTipChord)+(0.5*self.bref*abs(geomRootChord-geomTipChord))

        #Aspect ratio Bref^2/Sref
        self._AR = (self.bref**2)/geomSref

        #aeroCLaRoot and aeroCLaTip are 1x4 vectors, containing coefficients for a cubic polynomial describing a curve representing
        #the function CL(Alpha) for a 2D airfoil. Once this polynomial is setup, we perform an analytical differentiation to obtain
        #curves for CL_alpha(Alpha), ie, the lift slope as a function of Alpha
        #

        # let's create n span stations at which LLT equations will be solved and define some fixed wing properties at each station
        self._thetaLoc    = numpy.linspace(1,n,num=n,endpoint=True)*numpy.pi/(2.*n)
        self._yLoc = (0.5*self._bref)*numpy.cos(self.thetaLoc)
        self._chordLoc    = geomRootChord + 2.*(geomTipChord - geomRootChord)*self.yLoc/self.bref
        self._aIncGeometricLoc = geomRootAinc - (2.*(geomRootAinc - geomTipAinc)*self.yLoc/self.bref)
        self._clPolyLoc      = numpy.zeros((4,n))    #setup lift slope curves for each station
        for i in range(n):
            self._clPolyLoc[:,i] = aeroCLaRoot[:] + 2.*(aeroCLaTip[:] - aeroCLaRoot[:])*self.yLoc[i]/self.bref

    def alphaGeometric(self, alpha, p=None, q=None):
        assert p is None and q is None, 'angular velocity not yet implemented'
        return alpha + self.aIncGeometricLoc
        
    def clLoc(self,alphaLoc):
        return self.clPolyLoc[0,:]*alphaLoc**3 + \
               self.clPolyLoc[1,:]*alphaLoc**2 + \
               self.clPolyLoc[2,:]*alphaLoc + \
               self.clPolyLoc[3,:]
