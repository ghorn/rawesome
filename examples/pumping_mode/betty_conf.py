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

def makeConf():
    return \
        {'g': 9.81,  #  gravitational constant #  [ m /s^2]
         'rho': 1.23,  #  density of the air #  [ kg/m^3]

         #[aero]
         'alpha0deg': 0,

         #ROLL DAMPING
         'rD': 0.3e-2 ,
         'pD': 0,#1e-3,
         'yD': 0,#1e-3,

         #WIND-TUNNEL PARAMETERS
         #Lift (report p. 67)
         'cLA': 5.786876,
         'cLe': -1.924*0,
         'cL0': 0.203530,

         #Drag (report p. 70)
         'cDA': 0.018751,
         'cDA2': 1.529989,
         'cDB2': 0,

         'cD0': 0.006952,

         #Roll (report p. 72)
         'cRB': -0.062,
         'cRAB': -0.271 ,
         'cRr': -5.637e-1,

         #Pitch (report p. 74)
         'cPA': 0.293,
         'cPe': -4.9766e-1,

         'cP0': 0.03,

         #Yaw (report p. 76)
         'cYB': 0.05,
         'cYAB': 0.229,
         'cYrudder': 0.5,

         #[kite]
         'mass':  4.5,

         #TAIL LENGTH
         'lT': 0.4,

         # how alpha/beta computed
         #     'alpha_beta_computation':'first_order',
         'alpha_beta_computation':'closed_form',

         'sref': 0.684,
         'bref': 2.904, #sqrt(sref*AR),
         'cref': 0.2512, #sqrt(sref/AR),

         'zt': -0.01,

         #INERTIA MATRIX (Kurt's direct measurements)
         'j1': 0.0163,
         'j31': 0.0006,
         'j2': 0.0078,
         'j3': 0.0229,

         #[carousel]
         #'rArm': 1.085 #(dixit Kurt),
         'rArm': 2,

         #Carousel Friction & inertia
         'jCarousel': 1e2,
         'cfric': 0,#100,

         # wind model
         'wind_model':{'name':'wind_shear',
                       'z0':100,
                       'zt_roughness': 0.1}
#         'wind_model':{'name':'constant'}
        }
