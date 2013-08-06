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

def aeroForcesTorques(dae, conf, v_bw_f, v_bw_b, w_bn_b, (eTe_f_1, eTe_f_2, eTe_f_3)):
    """
    
    Parameters
    ----------
    
    v_bw (vector): Linear velocity of aircraft w.r.t wind carrying frame
        a.k.a. airfoil velocity
        a.k.a. relative wind velocity
    
    v_bw_f (projected vector): Linear velocity of aircraft w.r.t wind carrying frame, expressed in a generic frame (f)
    v_bw_b (projected vector): Linear velocity of aircraft w.r.t wind carrying frame, expressed in body frame
    w_bn_b (projected vector): Rotational velocity of aircraft w.r.t 
    
    (eTe_f_1, eTe_f_2, eTe_f_3) : Transverse axis of the airplane (y-axis) expressed in a generic frame (f)
        Used to define what component of aeroforces amounts to sideforce
        
        
    Returns
    --------
      (f_f_1, f_f_2, f_f_3, t_b_1, t_b_2, t_b_3)
      
      First three entries denote forces, expressed in frame (f)
      Last three entries denote moments, expressed in body frame (b)
       
    """
    eTe_f = C.veccat([eTe_f_1, eTe_f_2, eTe_f_3])
    rho = conf['rho']
    alpha0 = conf['alpha0deg']*C.pi/180

    sref = conf['sref']
    bref = conf['bref']
    cref = conf['cref']

    vKite2 = C.mul(v_bw_f.T, v_bw_f) #Airfoil speed^2
    vKite = C.sqrt(vKite2) #Airfoil speed
    dae['airspeed'] = vKite

    # Lift axis, normed to airspeed
    eLe_v_f = C.cross(eTe_f, v_bw_f)

    # sideforce axis, normalized to airspeed^2
    eYe_v2_f = C.cross(eLe_v_f, -v_bw_f)

    if conf['alpha_beta_computation'] == 'first_order':
        alpha = alpha0 + v_bw_b[2] / v_bw_b[0]
        beta = v_bw_b[1] / v_bw_b[0]
#        beta = v_bw_b_y / C.sqrt(v_bw_b[0]*v_bw_b[0] + v_bw_b[2]*v_bw_b[2])
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
    cD += conf['cD_elev2']*dae['elevator']*dae['elevator'] + conf['cD_A_elev']*alpha*dae['elevator'] + conf['cD_elev']*dae['elevator']
    cD += conf['cD_ail2']*dae['aileron']*dae['aileron'] + conf['cD_B_ail']*beta*dae['aileron'] + conf['cD_ail']*dae['aileron']
    if 'flaps' in dae:
        cL += conf['cL_flaps']*dae['flaps']
        cD += conf['cD_flaps2']*dae['flaps']*dae['flaps'] + conf['cD_A_flaps']*alpha*dae['flaps'] + conf['cD_flaps']*dae['flaps']
    if 'rudder' in dae:
        cY += conf['cY_rudder']*dae['rudder']
        cD += conf['cD_rudder2']*dae['rudder']*dae['rudder'] + conf['cD_B_rudder']*beta*dae['rudder'] + conf['cD_rudder']*dae['rudder']

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
    w_bn_b_hat = C.veccat([conf['bref'],conf['cref'],conf['bref']])*0.5/dae['airspeed']*w_bn_b

    momentCoeffs_pqr = C.mul(C.blockcat([[conf['cl_p'], conf['cl_q'], conf['cl_r']],
                                         [conf['cm_p'], conf['cm_q'], conf['cm_r']],
                                         [conf['cn_p'], conf['cn_q'], conf['cn_r']]]),
                             w_bn_b_hat)
    dae['momentCoeffs_pqr'] = momentCoeffs_pqr

    # with alpha beta
    momentCoeffs_AB = C.mul(C.blockcat([[           0, conf['cl_B'], conf['cl_AB']],
                                        [conf['cm_A'],            0,             0],
                                        [           0, conf['cn_B'], conf['cn_AB']]]),
                            C.vertcat([alpha, beta, alpha*beta]))
    dae['momentCoeffs_AB'] = momentCoeffs_AB

    # with control surfaces
    momentCoeffs_surf = C.SXMatrix(3, 1, 0)
    momentCoeffs_surf[0] += conf['cl_ail']*dae['aileron']
    momentCoeffs_surf[1] += conf['cm_elev']*dae['elevator']
    if 'flaps' in dae:
        momentCoeffs_surf[1] += conf['cm_flaps']*dae['flaps']
    if 'rudder' in dae:
        momentCoeffs_surf[2] += conf['cn_rudder']*dae['rudder']
    dae['momentCoeffs_surf'] = momentCoeffs_surf

    momentCoeffs = dae['momentCoeffs0'] + dae['momentCoeffs_pqr'] + \
                   dae['momentCoeffs_AB'] + dae['momentCoeffs_surf']
    dae['cl_small'] = momentCoeffs[0]
    dae['cm_small'] = momentCoeffs[1]
    dae['cn_small'] = momentCoeffs[2]

    # Common subexpressionf for aerodynamics
    rho_sref    = 0.5*rho*sref
    rho_sref_v2 = rho_sref*vKite2
    rho_sref_v  = rho_sref*vKite
    
    # LIFT :
    dae['fL'] = rho_sref_v2*cL
    fL_f = rho_sref_v*cL*eLe_v_f

    # DRAG :
    dae['fD'] = rho_sref_v2*cD
    fD_f = -rho_sref_v*cD*v_bw_f

    # sideforce
    dae['fY'] = rho_sref_v2*dae['cY']
    fY_f = rho_sref*cY*eYe_v2_f

    # aero forces
    f_f = fL_f + fD_f + fY_f

    # aero torques expressed in body frame
    t_b_1 =  rho_sref_v2*bref*dae['cl_small']
    t_b_2 =  rho_sref_v2*cref*dae['cm_small']
    t_b_3 =  rho_sref_v2*bref*dae['cn_small']

    dae['aero_fx'] = f_f[0]
    dae['aero_fy'] = f_f[1]
    dae['aero_fz'] = f_f[2]
    dae['aero_mx'] = t_b_1
    dae['aero_my'] = t_b_2
    dae['aero_mz'] = t_b_3
    return (f_f[0], f_f[1], f_f[2], t_b_1, t_b_2, t_b_3)

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
