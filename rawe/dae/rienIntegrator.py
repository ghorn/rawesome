import ctypes
import os
import subprocess

rawesomeDataPath = os.path.expanduser("~/.rawesome")

def runExporter(dae):
    # get the filename of the shared object
    filename = __file__.rstrip('.pyc')
    filename = filename.rstrip('.py')
    filename = filename.rstrip('rienIntegrator')
    filename += 'rienIntegratorInterface/rienIntegratorInterface.so'

    # load the shared object
    print 'loading "'+filename+'"'
    lib = ctypes.cdll.LoadLibrary(filename)

#    import numpy
#    x = numpy.array([[1,2,3,4],[5,6,7,8]], dtype=numpy.double)
#    y = numpy.zeros(x.shape)
#    print x
#    print y
#    xIn  = ctypes.c_void_p(x.ctypes.data)
#    xOut = ctypes.c_void_p(y.ctypes.data)

    # set some options
    numIntervals = 1
    timestep = 0.1
    numIntegratorSteps = 100
    integratorType = ctypes.c_char_p("INT_IRK_RIIA3")
    integratorGrid = None

    nx = len(dae.xNames())
    nz = len(dae.zNames())
    nup = len(dae.uNames()) + len(dae.pNames())

    timestep = ctypes.c_double(timestep)

    # set up the path for the rawesome directory
    integratorExportPath = rawesomeDataPath + "/exportIntegrator"
    if not os.path.exists(integratorExportPath):
        os.makedirs(integratorExportPath)

    # call makeRienIntegrator
    ret = lib.makeRienIntegrator(ctypes.c_char_p(integratorExportPath),
                                 numIntervals, 
                                 timestep, 
                                 integratorType,
                                 integratorGrid,
                                 numIntegratorSteps,
                                 nx, nz, nup)
    if ret != 0:
        raise Exception("rien integrator creater, what goon set bad options?")

    # write the model file
    modelFile = dae.acadoSimGen()[0]
    f = open(integratorExportPath+"/model.c",'w')
    f.write(modelFile)
    f.close()

    # compile the code
    p = subprocess.Popen(["make"], stdout=subprocess.PIPE, cwd=integratorExportPath)
    p.wait()
    ret = p.stdout.read()
    print "ret: " +str(ret)
