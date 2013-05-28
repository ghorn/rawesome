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

import pickle
import matplotlib.pyplot as plt

if __name__=='__main__':
    filename = "data/crosswind_opt.dat"
    
    f=open(filename,'r')
    traj = pickle.load(f)
    f.close()

    print "differential states: "+str(traj.xNames)
    print "algebraic states:    "+str(traj.zNames)
    print "actions:             "+str(traj.uNames)
    print "parameters:          "+str(traj.pNames)
    print "outputs:             "+str(traj.outputNames)

    traj.subplot([['aileron','elevator'],['daileron','delevator']],title='control surfaces')
    traj.subplot(['r','dr','ddr'])
    traj.subplot(['wind at altitude','dr','dx'])
    traj.subplot(['c','cdot','cddot'],title="invariants")
    traj.plot('airspeed',title='airspeed')
    traj.subplot([['alpha(deg)','alphaTail(deg)'],['beta(deg)','betaTail(deg)']])
    traj.subplot(['cL','cD','L/D'])
    traj.subplot(['winch power', 'tether tension'])
    traj.plot(["loyd's limit","loyd's limit (exact)","-(winch power)"])
    traj.subplot([['daileronCost','delevatorCost','ddrCost'],['winch power']])
    
    plt.show()
