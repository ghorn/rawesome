import numpy

import casadi as C

import rawe

if __name__=='__main__':
    from conf import conf
    
    print "creating model..."
    dae = rawe.models.carousel(conf)
    
    blah = dae.octaveSimGen('carouselOde')
    f = open('carouselOde_modelAndJacob.m','w')
    f.write(blah)
    f.close()
#    (modelFile, simExportFile) = 
#    print modelFile
#    print simExportFile
    import sys; sys.exit()

    (modelFile, simExportFile) = dae.acadoSimGen()
    print modelFile
    print simExportFile

    f = open('auto_model.c','w')
    f.write(modelFile)
    f.close()

    f = open('auto_sim_export.c','w')
    f.write(simExportFile)
    f.close()

    modelFile = dae.acadoModelGen()
    print modelFile
    print simExportFile

    f = open('auto_model.cpp','w')
    f.write(modelFile)
    f.close()
