# Copyright 2012-2013 Greg Horn
#
# This file is part of rawesome.
#
# rawesome is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rawesome is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with rawesome.  If not, see <http://www.gnu.org/licenses/>.

import rawe

if __name__=='__main__':
    from highwind_carousel_conf import conf
    
    print "creating model..."
    dae = rawe.models.carousel(conf)
   
    blah = dae.octaveSimGen('carouselOde')
    f = open('carouselOde_modelAndJacob.m','w')
    f.write(blah)
    f.close()

    import sys; sys.exit()
    (modelFile, simExportFile) = dae.acadoSimGen()

    f = open('data/auto_model.cpp','w')
    f.write(modelFile)
    f.close()

    f = open('data/auto_sim_export.c','w')
    f.write(simExportFile)
    f.close()

    modelFile = dae.acadoModelGen()
    print modelFile
    print simExportFile

    f = open('auto_model.cpp','w')
    f.write(modelFile)
    f.close()
