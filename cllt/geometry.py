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

    def __init__(self, thetaLoc, chordLoc, yLoc, aIncGeometricLoc, clPolyLoc, bref, n):
        '''
        n is the number of spanwise stations (and hence number of Fourier series terms)
        that we will divide each half wing into
        '''
        self._n = n
        self._sumN = numpy.linspace(1, 2*self.n-1, self.n)
        self._yLoc = yLoc
        self._thetaLoc = thetaLoc
        self._chordLoc = chordLoc
        self._aIncGeometricLoc = aIncGeometricLoc
        self._clPolyLoc = clPolyLoc

        #reference span for A/C
        self._bref = bref

        #reference wing surface area (projected)
        self.sref = 0
        for k in range(self.n-1):
            self.sref += (self.yLoc[k]-self.yLoc[k+1])*0.5*(self.chordLoc[k+1]+self.chordLoc[k])
        self.sref *= 2 # half-span to full-span

        #reference chord for A/c
        self._cref = self.sref / self.bref

        # Aspect ratio Bref^2/Sref
        self._AR = self.bref**2 / self.sref

    def alphaGeometric(self, alpha, p=None, q=None):
        assert p is None and q is None, 'angular velocity not yet implemented'
        return alpha + self.aIncGeometricLoc

    def clLoc(self,alphaLoc):
        return self.clPolyLoc[0,:]*alphaLoc**3 + \
               self.clPolyLoc[1,:]*alphaLoc**2 + \
               self.clPolyLoc[2,:]*alphaLoc + \
               self.clPolyLoc[3,:]


def simpleGeometry(geomRoot, geomTip, aeroCLaRoot, aeroCLaTip, n):
    geomRootChord = geomRoot[0]
    geomRootY     = geomRoot[1]
    geomRootAinc  = numpy.radians(geomRoot[2])
    geomTipChord  = geomTip[0]
    geomTipY      = geomTip[1]
    geomTipAinc   = numpy.radians(geomTip[2])

    bref = 2*(geomTipY-geomRootY)
    thetaLoc = numpy.linspace(1,n,num=n,endpoint=True)*numpy.pi/(2.0*n)
    yLoc = 0.5*bref*numpy.cos(thetaLoc)
    chordLoc = geomRootChord + 2.*(geomTipChord - geomRootChord)*yLoc/bref
    aIncGeometricLoc = geomRootAinc - 2.0*(geomRootAinc - geomTipAinc)*yLoc/bref

    clPolyLoc = numpy.zeros((4,n))    #setup lift slope curves for each station
    for i in range(n):
        clPolyLoc[:,i] = aeroCLaRoot[:] + 2.0*(aeroCLaTip[:] - aeroCLaRoot[:])*yLoc[i]/bref

    return Geometry(thetaLoc, chordLoc, yLoc, aIncGeometricLoc, clPolyLoc, bref, n)
