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

import casadi as C

def getWindAnglesFrom_v_bw_b(airspeed, v_bw_b):
    alpha =  C.arctan2(v_bw_b[2], v_bw_b[0] )
    beta  =  C.arcsin (v_bw_b[1] / airspeed )
    return (alpha, beta)

def aeroForcesTorques(dae, conf, v_bw_n, v_bw_b, (w1,w2,w3), (eTe1, eTe2, eTe3), (aileron,elevator)):
    rho = conf['rho']
    alpha0 = conf['alpha0deg']*C.pi/180

    sref = conf['sref']
    bref = conf['bref']
    cref = conf['cref']

    #Airfoil speed in carousel frame
    v_bw_n_x = v_bw_n[0]
    v_bw_n_y = v_bw_n[1]
    v_bw_n_z = v_bw_n[2]

    vKite2 = C.mul(v_bw_n.trans(), v_bw_n) #Airfoil speed^2
    vKite = C.sqrt(vKite2) #Airfoil speed
    dae['airspeed'] = vKite

    # Lift axis ** Normed to we @>@> **
    eLe = C.cross(C.veccat([eTe1,eTe2,eTe3]), v_bw_n)
    eLe1 = eLe[0]
    eLe2 = eLe[1]
    eLe3 = eLe[2]

    # Relative wind speed in Airfoil's referential 'E'
    v_bw_b_x = v_bw_b[0]
    v_bw_b_y = v_bw_b[1]
    v_bw_b_z = v_bw_b[2]

    if conf['alpha_beta_computation'] == 'first_order':
        alpha = alpha0 + v_bw_b_z / v_bw_b_x
        beta = v_bw_b_y / v_bw_b_x
#        beta = v_bw_b_y / C.sqrt(v_bw_b_x*v_bw_b_x + v_bw_b_z*v_bw_b_z)
    elif conf['alpha_beta_computation'] == 'closed_form':
        (alpha, beta) = getWindAnglesFrom_v_bw_b(vKite, v_bw_b)
        alpha += alpha0
    else:
        raise ValueError('config "alpha_beta_compuation" value '+str(conf['alpha_beta_computation'])+' not recognized, use "first_order" or "closed_form"')

    dae['alpha_deg'] = alpha*180/C.pi
    dae['beta_deg'] = beta*180/C.pi

    ########### force coefficients ###########
    # with alpha/beta
    cL = conf['cL_A']*alpha + conf['cL0']
    cD = conf['cD_A']*alpha + conf['cD_A2']*alpha*alpha + conf['cD_B2']*beta*beta + conf['cD0']
    cY = conf['cY_B']*beta

    # with control surfaces
    cL += conf['cL_elev']*elevator
    print "PUT CD_{FLAPS,RUDDER,AILERONS,ELEVATOR} IN (please)"
    if 'flaps' in dae:
        cL += conf['cL_flaps']*dae['flaps']
    if 'rudder' in dae:
        cY += conf['cY_rudder']*dae['rudder']

    dae['cL'] = cL
    dae['cD'] = cD
    dae['cD_tether'] = 0.25*dae['r']*0.001/sref
    dae['L_over_D'] = cL/cD
    cD = dae['cD'] + dae['cD_tether']
    dae['L_over_D_with_tether'] = cL/cD


    ######## moment coefficients #######
    # offset
    momentCoeffs = C.DMatrix([0, conf['cm0'], 0])

    # with roll rates
    momentCoeffs += C.mul(C.vertcat([C.horzcat([conf['cl_p'],conf['cl_q'],conf['cl_r']]),
                                     C.horzcat([conf['cm_p'],conf['cm_q'],conf['cm_r']]),
                                     C.horzcat([conf['cn_p'],conf['cn_q'],conf['cn_r']])]),
                          dae['w_bn_b'])

    # with alpha beta
    momentCoeffs += C.mul(C.vertcat([C.horzcat([           0, conf['cl_B'], conf['cl_AB']]),
                                     C.horzcat([conf['cm_A'],            0,             0]),
                                     C.horzcat([           0, conf['cn_B'], conf['cn_AB']])]),
                          C.vertcat([alpha, beta, alpha*beta]))

    # with control surfaces
    momentCoeffs[0] += conf['cl_ail']*dae['aileron']
    momentCoeffs[1] += conf['cm_elev']*dae['elevator']
    if 'flaps' in dae:
        momentCoeffs[1] += conf['cm_flaps']*dae['flaps']
    if 'rudder' in dae:
        momentCoeffs[2] += conf['cn_rudder']*dae['rudder']

    dae['cl'] = momentCoeffs[0]
    dae['cm'] = momentCoeffs[1]
    dae['cn'] = momentCoeffs[2]



    # LIFT :
    # ###############################
    dae['fL'] = sref*rho*cL*vKite2/2.0
    fL1 =  sref*rho*cL*eLe1*vKite/2.0
    fL2 =  sref*rho*cL*eLe2*vKite/2.0
    fL3 =  sref*rho*cL*eLe3*vKite/2.0

    # DRAG :
    # #############################
    dae['fD'] = sref*rho*vKite*cD*2.0
    fD1 = -sref*rho*vKite*cD*v_bw_n_x/2.0
    fD2 = -sref*rho*vKite*cD*v_bw_n_y/2.0
    fD3 = -sref*rho*vKite*cD*v_bw_n_z/2.0

    # aero forces
    f1 = fL1 + fD1
    f2 = fL2 + fD2
    f3 = fL3 + fD3

    # aero torques
    t1 =  0.5*rho*vKite2*sref*bref*dae['cl']
    t2 =  0.5*rho*vKite2*sref*cref*dae['cm']
    t3 =  0.5*rho*vKite2*sref*bref*dae['cn']

    dae['aero_fx'] = f1
    dae['aero_fy'] = f2
    dae['aero_fz'] = f3
    dae['aero_mx'] = t1
    dae['aero_my'] = t2
    dae['aero_mz'] = t3
    return (f1, f2, f3, t1, t2, t3)

#def tetherDragLumped( conf, r_n2b_n, v_bn_n, get_v_wn_n, tetherLength, N=10 ):
#    Cd = 1.1 # meh
#
#    forceTotal = C.DMatrix([0.0, 0.0, 0.0])
#    sumMe = 0
#    for k in range(N):
#        scaling = (k + 0.5)/N
#        sumM += scaling
#        r_n2t_n = r_n2b_n*scaling # local tether position
#        v_tn_n  =  v_bn_n*scaling
#        v_tw_b = v_tn_n - get_v_wn_n(r_n2t_n) # local tether velocity in wind frame
#        localV2 = C.mul(v_tw_b.T, v_tw_b)
#        localV = C.sqrt(localV2)
#
#        deltaH = scaling*r_n2b_n
#        hProjected = C.sqrt( C.mul(deltaH.T, deltaH) - C.mul(deltaH.T, v_tw_b)/localV
#
#        localForce = 0.5*rho*localV2*Cd*hProjected
#    raise Exception('summe: ',sumMe)
