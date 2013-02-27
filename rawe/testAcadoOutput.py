import models
import casadi as C
from casadi import *
import numpy
from config import readConfig

if __name__=='__main__':
    print "reading config..."
    conf = readConfig('config.ini','configspec.ini')
    
    print "creating model..."
    dae = models.crosswind(conf)
    
    (modelFile, simExportFile) = dae.acadoSimGen()
    print modelFile
    print simExportFile

#    f = open('auto_model.c','w')
#    f.write(modelFile)
#    f.close()
#
#    f = open('auto_sim_export.c','w')
#    f.write(simExportFile)
#    f.close()

    modelFile = dae.acadoModelGen()
    print modelFile
    print simExportFile

#    f = open('auto_model.cpp','w')
#    f.write(modelFile)
#    f.close()
