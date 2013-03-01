import models
import casadi as C
from casadi import *
import numpy
from config import readConfig

if __name__=='__main__':
    print "reading config..."
    conf = readConfig('config.ini','configspec.ini')
    
    print "creating model..."
    dae = models.carousel(conf)
    
    modelFile = dae.acadoModelGen()
    
    print "saving model..."
    f = open('TestCarousel/model.cpp','w')
    f.write(modelFile)
    f.close()
