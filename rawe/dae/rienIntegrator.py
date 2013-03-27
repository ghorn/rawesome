import ctypes
import os
import subprocess
import tempfile
import shutil
import hashlib

rawesomeDataPath = os.path.expanduser("~/.rawesome")

def loadIntegratorInterface():
    # get the filename of the shared object
    filename = __file__.rstrip('.pyc')
    filename = filename.rstrip('.py')
    filename = filename.rstrip('rienIntegrator')
    filename += 'rienIntegratorInterface/rienIntegratorInterface.so'

    # load the shared object
    return ctypes.cdll.LoadLibrary(filename)

def makeMakefile(cfiles):
    return """\
CC      = gcc
CFLAGS  = -O3 -fPIC -finline-functions -I.
LDFLAGS = -lm

C_SRC = %(cfiles)s

.PHONY: clean lib obj

%%.o : %%.c acado.h
\t@echo CC $@: $(CC) $(CFLAGS) -c $< -o $@
\t@$(CC) $(CFLAGS) -c $< -o $@

OBJ = $(C_SRC:%%.c=%%.o)
integrator.so : $(OBJ)
\t@echo LD $@: $(CC) -shared -Wl,-soname,$@ -o $@ $(OBJ) $(LDFLAGS)
\t@$(CC)   -shared -Wl,-soname,$@ -o $@ $(OBJ) $(LDFLAGS)

# aliases
obj : $(OBJ)
lib : integrator.so
clean :
\trm -f *.o *.so
""" % {'cfiles':' '.join(cfiles)}


def writeRienIntegrator(dae, path):
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

    # call makeRienIntegrator
    lib = loadIntegratorInterface()
    ret = lib.makeRienIntegrator(ctypes.c_char_p(path),
                                 numIntervals, 
                                 timestep, 
                                 integratorType,
                                 integratorGrid,
                                 numIntegratorSteps,
                                 nx, nz, nup)
    if ret != 0:
        raise Exception("rien integrator creater, what goon set bad options?")


def generateRienIntegrator(dae):
    # make temporary directory and generate integrator.c and acado.h there
    # then read those files
    tmppath = tempfile.mkdtemp()
    try:
        writeRienIntegrator(dae,tmppath)

        integratorPath = tmppath+'/integrator.c'
        acadoHeaderPath = tmppath+'/acado.h'

        f = open(integratorPath,'r')
        integrator = f.read()
        f.close()

        f = open(acadoHeaderPath,'r')
        acadoHeader = f.read()
        f.close()

    finally:
        shutil.rmtree(tmppath)

    return (integrator, acadoHeader)
    

def exportIntegrator(dae):
    # make ~/.rawesome if it doesn't exist
    if not os.path.exists(rawesomeDataPath):
        os.makedirs(rawesomeDataPath)

    # get the exported integrator files
    (integrator, acadoHeader) = generateRienIntegrator(dae)

    # model file
    modelFile = dae.acadoSimGen()[0]

    # write the makefile
    makefile = makeMakefile(['workspace.c', 'model.c', 'integrator.c'])

    # write the static workspace file (temporary)
    workspace = """\
#include <acado.h>
ACADOworkspace acadoWorkspace;
ACADOvariables acadoVariables;
"""

    genfiles = [('integrator.c', integrator),
                ('acado.h', acadoHeader),
                ('model.c', modelFile),
                ('workspace.c', workspace),
                ('Makefile', makefile)]

    # hash the files
    exportpath = rawesomeDataPath + '/' + \
        hashlib.md5(''.join([''.join(x) for x in genfiles])).hexdigest()

    def writeFiles():
        for (filename, filestring) in genfiles:
            f = open(exportpath+'/'+filename,'w')
            f.write(filestring)
            f.close()
    def compileFiles():
        # compile the code
        p = subprocess.Popen(['make'], cwd=exportpath)
        ret = p.wait()
        if ret != 0:
            raise Excepetion("integrator compilation failed, return code "+str(ret))

        #p = subprocess.Popen(['make'], stdout=subprocess.PIPE, cwd=exportpath)
        #p.wait()
        #ret = p.stdout.read()
        #print "ret: " +str(ret)

    # if no directory named by this hash exists, create it and compile the project there
    print "integrator export path:  "+exportpath
    if not os.path.exists(exportpath):
        os.makedirs(exportpath)
        writeFiles()
        compileFiles()

    # if the directory already exists, check if the contents match
    else:
        unmatched = []
        for (filename, filestring) in genfiles:
            f = open(exportpath+'/'+filename, 'r')
            contents = f.read()
            f.close()
            if contents != filestring:
                unmatched += filename
        # if the files match, run "make" to ensure everything is built
        if len(unmatched) == 0:
            compileFiles()
        # if the files don't match, print message and write new files
        else:
            print "shit! got hash collision"
            # remove any existing files
            for f in os.listdir(exportpath):
                os.remove(exportpath+'/'+f)
            # write and compile our new files
            writeFiles()
            compileFiles()

    return exportpath+'/integrator.so'


def runExporter(dae):
    libpath = exportIntegrator(dae)

    print 'loading '+libpath
    lib = ctypes.cdll.LoadLibrary(libpath)
#    import numpy
#    x = numpy.array([[1,2,3,4],[5,6,7,8]], dtype=numpy.double)
#    y = numpy.zeros(x.shape)
#    print x
#    print y
#    xIn  = ctypes.c_void_p(x.ctypes.data)
#    xOut = ctypes.c_void_p(y.ctypes.data)

