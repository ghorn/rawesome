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

from numpy import pi

def makeConf():
    return \
        {'g': 9.81,  #  gravitational constant #  [ m /s^2]
         'rho': 1.23,  #  density of the air #  [ kg/m^3]

         'alpha0deg': 0,

         # forces
         'cL0': 0.203530,
         'cL_A':5.786876,

         'cD_A': 0.018751,
         'cD_A2': 1.529989,
         'cD_B2': -0.16247,
         'cD0': 0.008767,

         'cY_B':-0.239789,

         # control surface forces
         'cL_elev':-0.0105*180/pi,
         'cL_flaps':0.0184*180/pi,
         'cY_rudder':0.0035*180/pi,

         # stability derivatives
         'cl_p':-0.576, 'cl_q':  0.0, 'cl_r': 0.0707,
         'cm_p':   0.0, 'cm_q':-15.5, 'cm_r':    0.0,
         'cn_p':-0.036, 'cn_q':  0.0, 'cn_r':-0.0667,

         'cl_B':-0.051808,
         'cl_AB':-0.208344,
         'cm_A':-0.450643,
         # cm0 valid for CG/bridle location 0.1 meters behind main wing leading edge
         'cm0':0.028980,
         'cn_B':0.037183,
         'cn_AB':-0.028933,

         # control surface moments
         'cl_ail':0.0073*180/pi,
         'cm_elev':0.0352*180/pi,
         'cm_flaps':0.0026*180/pi,
         'cn_rudder':0.001176*180/pi,

         #[kite]
         'mass':  4.5,

         # how alpha/beta computed
         #'alpha_beta_computation':'first_order',
         'alpha_beta_computation':'closed_form',

         'sref': 0.684,
         'bref': 2.904, #sqrt(sref*AR),
         'cref': 0.2512, #sqrt(sref/AR),

         'zt': 0.10,

         # INERTIA MATRIX: THIS IS A WILD GUESS
         'j1': 0.565,
         'j31': 0.0,
         'j2': 0.161,
         'j3': 0.723,

         #[carousel]
         #'rArm': 1.085 #(dixit Kurt),
         'rArm': 2,

         #Carousel Friction & inertia
         'jCarousel': 1e2,
         'cfric': 0,#100,

         # carousel
         'delta_parameterization':'cos_sin',

         # wind model
#         'wind_model':{'name':'wind_shear',
#                       'z0':100,
#                       'zt_roughness': 0.1}
         'wind_model':{'name':'constant'}
        }
