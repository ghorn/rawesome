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

def aeroForcesTorques(dae, conf, v_bw_n, v_bw_b, (w1,w2,w3), (eTe1, eTe2, eTe3)):
    rho = conf['rho']
    alpha0 = conf['alpha0deg']*C.pi/180

    sref = conf['sref']
    bref = conf['bref']
    cref = conf['cref']

    # airfoil speed in wind frame
    v_bw_n_x = v_bw_n[0]
    v_bw_n_y = v_bw_n[1]
    v_bw_n_z = v_bw_n[2]

    vKite2 = C.mul(v_bw_n.trans(), v_bw_n) #Airfoil speed^2
    vKite = C.sqrt(vKite2) #Airfoil speed
    dae['airspeed'] = vKite

    # Lift axis, normed to airspeed
    eLe_v = C.cross(C.veccat([eTe1,eTe2,eTe3]), v_bw_n)

    # sideforce axis, normalized to airspeed^2
    eYe_v2 = C.cross(eLe_v, -v_bw_n)

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
    cL += conf['cL_elev']*dae['elevator']
    print "PUT CD_{FLAPS,RUDDER,AILERONS,ELEVATOR} IN (please)"
    if 'flaps' in dae:
        cL += conf['cL_flaps']*dae['flaps']
    if 'rudder' in dae:
        cY += conf['cY_rudder']*dae['rudder']

    dae['cL'] = cL
    dae['cD'] = cD
    dae['cY'] = cY
    dae['cD_tether'] = 0.25*dae['r']*0.001/sref
    dae['L_over_D'] = cL/cD
    cD = dae['cD'] + dae['cD_tether']
    dae['L_over_D_with_tether'] = cL/cD


    ######## moment coefficients #######
    # offset
    dae['momentCoeffs0'] = C.DMatrix([0, conf['cm0'], 0])

    # with roll rates
    # non-dimensionalized angular velocity
    w_bn_b_hat = C.veccat([0.5*conf['bref']/dae['airspeed']*dae['w_bn_b_x'],
                           0.5*conf['cref']/dae['airspeed']*dae['w_bn_b_y'],
                           0.5*conf['bref']/dae['airspeed']*dae['w_bn_b_z']])
    momentCoeffs_pqr = C.mul(C.vertcat([C.horzcat([conf['cl_p'],conf['cl_q'],conf['cl_r']]),
                                        C.horzcat([conf['cm_p'],conf['cm_q'],conf['cm_r']]),
                                        C.horzcat([conf['cn_p'],conf['cn_q'],conf['cn_r']])]),
                             w_bn_b_hat)
    dae['momentCoeffs_pqr'] = momentCoeffs_pqr

    # with alpha beta
    momentCoeffs_AB = C.mul(C.vertcat([C.horzcat([           0, conf['cl_B'], conf['cl_AB']]),
                                       C.horzcat([conf['cm_A'],            0,             0]),
                                       C.horzcat([           0, conf['cn_B'], conf['cn_AB']])]),
                            C.vertcat([alpha, beta, alpha*beta]))
    dae['momentCoeffs_AB'] = momentCoeffs_AB

    # with control surfaces
    momentCoeffs_surf = C.SXMatrix(3,1,0)
    momentCoeffs_surf[0] += conf['cl_ail']*dae['aileron']
    momentCoeffs_surf[1] += conf['cm_elev']*dae['elevator']
    if 'flaps' in dae:
        momentCoeffs_surf[1] += conf['cm_flaps']*dae['flaps']
    if 'rudder' in dae:
        momentCoeffs_surf[2] += conf['cn_rudder']*dae['rudder']
    dae['momentCoeffs_surf'] = momentCoeffs_surf

    momentCoeffs = dae['momentCoeffs0'] + dae['momentCoeffs_pqr'] + \
                   dae['momentCoeffs_AB'] + dae['momentCoeffs_surf']
    dae['cl'] = momentCoeffs[0]
    dae['cm'] = momentCoeffs[1]
    dae['cn'] = momentCoeffs[2]



    # LIFT :
    dae['fL'] = 0.5*rho*vKite2*sref*cL
    fLx =  0.5*rho*vKite*sref*cL*eLe_v[0]
    fLy =  0.5*rho*vKite*sref*cL*eLe_v[1]
    fLz =  0.5*rho*vKite*sref*cL*eLe_v[2]

    # DRAG :
    dae['fD'] = 0.5*rho*vKite2*sref*cD
    fDx = -0.5*rho*sref*vKite*cD*v_bw_n_x
    fDy = -0.5*rho*sref*vKite*cD*v_bw_n_y
    fDz = -0.5*rho*sref*vKite*cD*v_bw_n_z

    # sideforce
    dae['fY'] = 0.5*rho*vKite2*sref*dae['cY']
    fYx = 0.5*rho*sref*cY*eYe_v2[0]
    fYy = 0.5*rho*sref*cY*eYe_v2[1]
    fYz = 0.5*rho*sref*cY*eYe_v2[2]

    # aero forces
    f1 = fLx + fDx + fYx
    f2 = fLy + fDy + fYy
    f3 = fLz + fDz + fYz

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
