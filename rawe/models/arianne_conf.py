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
         # these are calculated by a matlab fitting script based on vlm alpha/beta sweeps
         'cL0': 0.145115,
         'cL_A':  5.474703,

         'cD_A': 0.027310,
         'cD_A2': 1.948549,
         'cD_B2': -0.013039,
         'cD0':  0.006600,

         'cY_B': -0.131871,

         # control surface forces
         #note sign change from avl to local coords
         'cL_elev': -0.0074*180/pi,
         
         # these are calculated by a matlab fitting script based on vlm alpha/beta sweeps
         'cD_elev2': 3.52135e-05, 'cD_A_elev': -0.000101006, 'cD_elev': -6.67268e-06,
         'cD_ail2':0.000120247, 'cD_B_ail':-1.89122e-05, 'cD_ail':0,



         # stability derivatives
         # these are average values across the vlm alpha/beta sweep
         'cl_p': -0.4928, 'cl_q':    0.0, 'cl_r': 0.0566,
         'cm_p':     0.0, 'cm_q':-16.737, 'cm_r':    0.0,
         'cn_p': -0.0223, 'cn_q':    0.0, 'cn_r':-0.0636,

         # these are calculated by a matlab fitting script based on vlm alpha/beta sweeps
         'cl_B': -0.057367,
         'cl_AB': -0.437572,

         'cm_A':-2.285041,
         # cm0 valid for CG/bridle location 0.1 meters behind main wing leading edge
         'cm0':-0.056816,

         'cn_B': 0.038296,
         'cn_AB':-0.052142,

         # control surface moments
         # these are average values across the vlm alpha/beta sweep
         #note sign change from avl to local coords
         'cl_ail':0.0056*180/pi,  
         'cm_elev':0.0262*180/pi,


         #[kite]
         'mass':  0.626,  #  mass of the kite   #  [ kg]
         'tether_mass': 0.0,

         # use aerodynamic approximations instead of triginometry
         'alpha_beta_computation':'first_order',
#         'alpha_beta_computation':'closed_form',

         # use cos(delta), sin(delta) as states
         'delta_parameterization':'cos_sin',
 	      	
         'sref': 0.1051,
         'bref': 0.96, #sqrt(sref*AR)
         'cref': 0.128, #sqrt(sref/AR)

         'zt': 0.01,

         #INERTIA MATRIX (Kurt's direct measurements)
         'j1': 0.0163,
         'j31': -0.0006,
         'j2': 0.0078,
         'j3': 0.0229,

         #[carousel]
         #'rArm': 1.085, #(dixit Kurt)
         'rArm': 2,

         #Carousel Friction & inertia
         'jCarousel': 1e2,
         'cfric': 0,#100

#         'wind_model':{'name':'wind_shear',
#                       'z0':100,
#                       'zt_roughness': 0.1}
#         'wind_model':{'name':'constant'}
        }
