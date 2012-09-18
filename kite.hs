{-# OPTIONS_GHC -Wall #-}

import Data.Packed
import Data.Packed.ST
import Numeric.Container hiding ( conj )
import Graphics.Gnuplot.Simple

import Numeric.GSL.ODE

modelInteg :: Double -> Vector Double -> Vector Double -> (Vector Double, Vector Double)
modelInteg r state u = (sys, fromList [c, cdot, cddot])
  where
    --  PARAMETERS OF THE KITE :
    --  -----------------------------
    m =  0.626      --  mass of the kite               --  [ kg    ]
                 
    --   PHYSICAL CONSTANTS :
    --  -----------------------------
    g =    9.81      --  gravitational constant         --  [ m /s^2]
    rho =    1.23      --  density of the air             --  [ kg/m^3]
    
    --  PARAMETERS OF THE CABLE :
    --  -----------------------------
     
    --CAROUSEL ARM LENGTH
    rA = 1.085 --(dixit Kurt)
    
    _ZT = -0.01
    
    --INERTIA MATRIX (Kurt's direct measurements)
    _J1 = 0.0163
    _J31 = 0.0006
    _J2 = 0.0078
    _J3 = 0.0229
    
    --Carousel Friction & inertia
    _I = 1e2
    _Cfric = 100
    
    alpha0 = -0*pi/180 
    
    --TAIL LENGTH
    _LT = 0.4
    
    --ROLL DAMPING
    _RD = 1e-2 
    _PD = 0*1e-3
    _YD = 0*1e-3
    
    --WIND-TUNNEL PARAMETERS
    --Lift (report p. 67)
    _CLA = 5.064
    
    _CLe = -1.924
    
    _CL0 = 0.239
    
    --Drag (report p. 70)
    _CDA = -0.195
    _CDA2 = 4.268
    _CDB2 = 5
    {-
    _CDe = 0.044
    _CDr = 0.111
    -}
    _CD0 = 0.026
    
    --Roll (report p. 72)
    _CRB = -0.062
    _CRAB = -0.271 
    _CRr = -5.637e-1
    
    --Pitch (report p. 74)
    _CPA = 0.293
    _CPe = -4.9766e-1
    
    _CP0 = 0.03
    
    --Yaw (report p. 76)
    _CYB = 0.05
    _CYAB = 0.229
    
    _SPAN = 0.96
    _CHORD = 0.1
    
    -----------------------     model integ ---------------------------------------
    x =   state @> 0
    y =   state @> 1
    z =   state @> 2
    
    e11 = state @> 3
    e12 = state @> 4
    e13 = state @> 5
    
    e21 = state @> 6
    e22 = state @> 7
    e23 = state @> 8
    
    e31 = state @> 9
    e32 = state @> 10
    e33 = state @> 11
    
    dx =  state @> 12
    dy =  state @> 13
    dz =  state @> 14
    
    w1 =  state @> 15
    w2 =  state @> 16
    w3 =  state @> 17

    {-
    delta = state @> 18
    -}
    ddelta = state @> 19
    
    {- u = linint(P.tu,t) -}
    _Tc = u @> 0 --Carousel motor torque
    u1 =  u @> 1
    u2 =  u @> 2
    

    {-
     r = P.r
     dr = 0
    -}
    
    --------------- kinfile -------------
    {-
    p = fromList [ cos(delta)*(rA + x) - y*sin(delta)
                 , sin(delta)*(rA + x) + y*cos(delta)
                 , z
                 ]
    dp = fromList [ dx*cos(delta) - ddelta*(sin(delta)*(rA + x) + y*cos(delta)) - dy*sin(delta)
                  , ddelta*(cos(delta)*(rA + x) - y*sin(delta)) + dy*cos(delta) + dx*sin(delta)
                  , dz
                  ]
    -}
    dpE = fromList [ dx*e11 + dy*e12 + dz*e13 + ddelta*e12*rA + ddelta*e12*x - ddelta*e11*y
                   , dx*e21 + dy*e22 + dz*e23 + ddelta*e22*rA + ddelta*e22*x - ddelta*e21*y
                   , dx*e31 + dy*e32 + dz*e33 + ddelta*e32*rA + ddelta*e32*x - ddelta*e31*y
                   ]
    dp_carousel_frame = fromList [ dx - ddelta*y
                                 , dy + ddelta*rA + ddelta*x
                                 , dz
                                 ]
    -----------------------------------------
    
    
    
    ---------- more model_integ ----------------------
    
    -- EFFECTIVE WIND IN THE KITE`S SYSTEM :
    -- ---------------------------------------------------------------
    
    --Airfoil speed in carousel frame
    we1 = dp_carousel_frame @> 0
    we2 = dp_carousel_frame @> 1
    we3 = dp_carousel_frame @> 2
    
    _VKite2 = dp_carousel_frame <.> dp_carousel_frame --Airfoil speed^2 
    _VKite = sqrt(_VKite2) --Airfoil speed
    
    -- CALCULATION OF THE FORCES :
    -- ---------------------------------------------------------------
    --
    --   FORCE ARE COMPUTED IN THE CAROUSEL FRAME !!!
    
    --Aero coeff.
    
            
    -- LIFT DIRECTION VECTOR
    -- -------------------------
    
    -- Relative wind speed in Airfoil's referential 'E'
    
    wE1 = dpE @> 0
    wE2 = dpE @> 1
    wE3 = dpE @> 2
    
    -- Airfoil's transversal axis in carousel referential 'e'
    eTe1 = e21
    eTe2 = e22
    eTe3 = e23
    
    
    -- Lift axis ** Normed to we !! **
    eLe1 = - eTe2*we3 + eTe3*we2
    eLe2 = - eTe3*we1 + eTe1*we3
    eLe3 = - eTe1*we2 + eTe2*we1
       
    
    -- AERODYNAMIC COEEFICIENTS
    -- ----------------------------------
    _VT1 =           wE1
    _VT2 = -_LT*w3 + wE2
    _VT3 =  _LT*w2 + wE3
    
    
    alpha = alpha0-wE3/wE1
    
    --NOTE: beta & alphaTail are compensated for the tail motion induced by
    --omega !!
    betaTail = _VT2/sqrt(_VT1*_VT1 + _VT3*_VT3)
    beta = wE2/sqrt(wE1*wE1 + wE3*wE3)
    alphaTail = alpha0-_VT3/_VT1
    
    -- _CL = _CLA*alpha + _CLe*u2   + _CL0
    -- _CD = _CDA*alpha + _CDA2*alpha*alpha + _CDB2*beta*beta + _CDe*u2 + _CDr*u1 + _CD0
    -- _CR = -_RD*w1 + _CRB*beta + _CRAB*alphaTail*beta + _CRr*u1
    -- _CP = _CPA*alphaTail + _CPe*u2  + _CP0
    -- _CY = _CYB*beta + _CYAB*alphaTail*beta
    
    _CL' = _CLA*alpha + _CLe*u2 + _CL0
    _CD' = _CDA*alpha + _CDA2*alpha*alpha + _CDB2*beta*beta + _CD0
    _CR = -_RD*w1 + _CRB*betaTail + _CRr*u1 + _CRAB*alphaTail*betaTail
    _CP = -_PD*w2 + _CPA*alphaTail + _CPe*u2 + _CP0
    _CY = -_YD*w3 + _CYB*betaTail + _CYAB*alphaTail*betaTail
    
    
    -- LIFT :
    -- ---------------------------------------------------------------
    _CL = 0.2*_CL'
    _CD = 0.5*_CD'
    _FL1 =  rho*_CL*eLe1*_VKite/2.0
    _FL2 =  rho*_CL*eLe2*_VKite/2.0
    _FL3 =  rho*_CL*eLe3*_VKite/2.0
    
    -- DRAG :
    -- -----------------------------------------------------------
    _FD1 = -rho*_VKite*_CD*we1/2.0
    _FD2 = -rho*_VKite*_CD*we2/2.0 
    _FD3 = -rho*_VKite*_CD*we3/2.0 
    
    
    
    -- FORCES (AERO)
    -- ---------------------------------------------------------------
    
    
    _F1 = _FL1 + _FD1
    _F2 = _FL2 + _FD2
    _F3 = _FL3 + _FD3
    
    --_F = _F-_FT
       
    -- TORQUES (AERO)
    -- ---------------------------------------------------------------
     
    _T1 =  0.5*rho*_VKite2*_SPAN*_CR
    _T2 =  0.5*rho*_VKite2*_CHORD*_CP
    _T3 =  0.5*rho*_VKite2*_SPAN*_CY
     
     
    -- ATTITUDE DYNAMICS
    -- -----------------------------------------------------------
    
    -- DynFile          -- Call DAE
    _MM :: Matrix Double
    _MM = runSTMatrix $ do
      _MM' <- newMatrix 0 8 8
      let writeMatrix' (row,col) = writeMatrix _MM' (row-1) (col-1)
      writeMatrix' (1,1) $ _I + m*rA*rA + m*x*x + m*y*y + 2*m*rA*x 
      writeMatrix' (1,2) $ -m*y 
      writeMatrix' (1,3) $ m*(rA + x) 
      writeMatrix' (1,4) $ 0 
      writeMatrix' (1,5) $ 0 
      writeMatrix' (1,6) $ 0 
      writeMatrix' (1,7) $ 0 
      writeMatrix' (1,8) $ 0
      
      writeMatrix' (2,1) $ -m*y 
      writeMatrix' (2,2) $ m 
      writeMatrix' (2,3) $ 0 
      writeMatrix' (2,4) $ 0 
      writeMatrix' (2,5) $ 0 
      writeMatrix' (2,6) $ 0 
      writeMatrix' (2,7) $ 0 
      writeMatrix' (2,8) $ x + _ZT*e31
      
      writeMatrix' (3,1) $ m*(rA + x) 
      writeMatrix' (3,2) $ 0 
      writeMatrix' (3,3) $ m 
      writeMatrix' (3,4) $ 0 
      writeMatrix' (3,5) $ 0 
      writeMatrix' (3,6) $ 0 
      writeMatrix' (3,7) $ 0 
      writeMatrix' (3,8) $ y + _ZT*e32
      
      writeMatrix' (4,1) $ 0 
      writeMatrix' (4,2) $ 0 
      writeMatrix' (4,3) $ 0 
      writeMatrix' (4,4) $ m 
      writeMatrix' (4,5) $ 0 
      writeMatrix' (4,6) $ 0 
      writeMatrix' (4,7) $ 0 
      writeMatrix' (4,8) $ z + _ZT*e33
      
      writeMatrix' (5,1) $ 0 
      writeMatrix' (5,2) $ 0 
      writeMatrix' (5,3) $ 0 
      writeMatrix' (5,4) $ 0 
      writeMatrix' (5,5) $ _J1 
      writeMatrix' (5,6) $ 0 
      writeMatrix' (5,7) $ _J31 
      writeMatrix' (5,8) $ -_ZT*(e21*x + e22*y + e23*z + _ZT*e21*e31 + _ZT*e22*e32 + _ZT*e23*e33)
      
      writeMatrix' (6,1) $ 0 
      writeMatrix' (6,2) $ 0 
      writeMatrix' (6,3) $ 0 
      writeMatrix' (6,4) $ 0 
      writeMatrix' (6,5) $ 0 
      writeMatrix' (6,6) $ _J2 
      writeMatrix' (6,7) $ 0 
      writeMatrix' (6,8) $ _ZT*(e11*x + e12*y + e13*z + _ZT*e11*e31 + _ZT*e12*e32 + _ZT*e13*e33)
      
      writeMatrix' (7,1) $ 0 
      writeMatrix' (7,2) $ 0 
      writeMatrix' (7,3) $ 0 
      writeMatrix' (7,4) $ 0 
      writeMatrix' (7,5) $ _J31 
      writeMatrix' (7,6) $ 0 
      writeMatrix' (7,7) $ _J3 
      writeMatrix' (7,8) $ 0
      
      writeMatrix' (8,1) $ -_ZT*(e11*e23*x - e13*e21*x + e12*e23*y - e13*e22*y + _ZT*e11*e23*e31 - _ZT*e13*e21*e31 + _ZT*e12*e23*e32 - _ZT*e13*e22*e32) 
      writeMatrix' (8,2) $ x + _ZT*e31 
      writeMatrix' (8,3) $ y + _ZT*e32 
      writeMatrix' (8,4) $ z + _ZT*e33 
      writeMatrix' (8,5) $ -_ZT*(e21*x + e22*y + e23*z + _ZT*e21*e31 + _ZT*e22*e32 + _ZT*e23*e33) 
      writeMatrix' (8,6) $ _ZT*(e11*x + e12*y + e13*z + _ZT*e11*e31 + _ZT*e12*e32 + _ZT*e13*e33) 
      writeMatrix' (8,7) $ 0 
      writeMatrix' (8,8) $ 0
      return _MM'
    
    conj = id
    _RHS :: Vector Double
    _RHS = fromList
           [ _Tc - _Cfric*ddelta - _F1*y + _F2*(rA + x) + dy*m*(dx - 2*ddelta*y) - dx*m*(dy + 2*ddelta*rA + 2*ddelta*x) 
           , _F1 + ddelta*m*(dy + ddelta*rA + ddelta*x) + ddelta*dy*m 
           , _F2 - ddelta*m*(dx - ddelta*y) - ddelta*dx*m 
           , _F3 - g*m 
           , _T1 - w2*(_J3*w3 + _J31*w1) + _J2*w2*w3 
           , _T2 + w1*(_J3*w3 + _J31*w1) - w3*(_J1*w1 + _J31*w3) 
           , _T3 + w2*(_J1*w1 + _J31*w3) - _J2*w1*w2 
           , (w1 - ddelta*e13)*(e21*(_ZT*dx - _ZT2*e21*(conj(w1) - ddelta*e13) + _ZT2*e11*(conj(w2) - ddelta*e23)) + e22*(_ZT*dy - _ZT2*e22*(conj(w1) - ddelta*e13) + _ZT2*e12*(conj(w2) - ddelta*e23)) + _ZT*e33*(z*conj(w1) + ddelta*e11*x + ddelta*e12*y + _ZT*e33*conj(w1) + _ZT*ddelta*e11*e31 + _ZT*ddelta*e12*e32) + _ZT*e23*(dz + _ZT*e13*conj(w2) - _ZT*e23*conj(w1)) + _ZT*e31*(conj(w1) - ddelta*e13)*(x + _ZT*e31) + _ZT*e32*(conj(w1) - ddelta*e13)*(y + _ZT*e32)) - dz*(dz + _ZT*e13*w2 - _ZT*e23*w1) - dx*(dx - _ZT*e21*(w1 - ddelta*e13) + _ZT*e11*(w2 - ddelta*e23)) - dy*(dy - _ZT*e22*(w1 - ddelta*e13) + _ZT*e12*(w2 - ddelta*e23)) - (_ZT*conj(w1)*(e11*x + e12*y + e13*z + _ZT*e11*e31 + _ZT*e12*e32 + _ZT*e13*e33) + _ZT*conj(w2)*(e21*x + e22*y + e23*z + _ZT*e21*e31 + _ZT*e22*e32 + _ZT*e23*e33))*(w3 - ddelta*e33) - (w2 - ddelta*e23)*(e11*(_ZT*dx - _ZT2*e21*(conj(w1) - ddelta*e13) + _ZT2*e11*(conj(w2) - ddelta*e23)) + e12*(_ZT*dy - _ZT2*e22*(conj(w1) - ddelta*e13) + _ZT2*e12*(conj(w2) - ddelta*e23)) - _ZT*e33*(z*conj(w2) + ddelta*e21*x + ddelta*e22*y + _ZT*e33*conj(w2) + _ZT*ddelta*e21*e31 + _ZT*ddelta*e22*e32) + _ZT*e13*(dz + _ZT*e13*conj(w2) - _ZT*e23*conj(w1)) - _ZT*e31*(conj(w2) - ddelta*e23)*(x + _ZT*e31) - _ZT*e32*(conj(w2) - ddelta*e23)*(y + _ZT*e32)) 
           ]
      where
        _ZT2 = _ZT*_ZT
    
    
    dRexp :: Matrix Double
    dRexp =
      fromLists [ [ e21*(w3 - ddelta*e33) - e31*(w2 - ddelta*e23) 
                  , e31*(w1 - ddelta*e13) - e11*(w3 - ddelta*e33) 
                  , e11*(w2 - ddelta*e23) - e21*(w1 - ddelta*e13) 
                  ]
                , [ e22*(w3 - ddelta*e33) - e32*(w2 - ddelta*e23) 
                  , e32*(w1 - ddelta*e13) - e12*(w3 - ddelta*e33) 
                  , e12*(w2 - ddelta*e23) - e22*(w1 - ddelta*e13) 
                  ]
                , [ e23*w3 - e33*w2 
                  , e33*w1 - e13*w3 
                  , e13*w2 - e23*w1
                  ]
                ]
    
    ddq_nu :: Vector Double
    ddq_nu = _MM <\> _RHS -- Backsolve DAE
    
    ----Extract dynamics
    --dddelta = ddq_nu(1)       -- Carousel acceleration
    --ddX = ddq_nu(2:4)         -- x,y,z acceleration
    --dw = ddq_nu(5:end-1)      -- Angular accelerations
    --nu = ddq_nu(end)          -- Algebraic state
    
{-    [dddelta, ddX, dw, nu] = takesV [1,3,3,1] ddq_nu -}
    [dddelta', ddX, dw, _] = takesV [1,3,3,1] ddq_nu
    
    --sys = [dx;dy;dz;reshape(dRexp,9,1);ddX;dw;ddelta;dddelta];
    sys = join [ fromList [dx,dy,dz]
               , flatten (trans dRexp)
               , ddX
               , dw
               , fromList [ddelta]
               , dddelta'
               ]

    (c,cdot,cddot) = (c',cdot', cddot')
      where
        dw1 = dw @> 0
        dw2 = dw @> 1
        {-
        dw3 = dw @> 2
        -}
        ddx = ddX @> 0
        ddy = ddX @> 1
        ddz = ddX @> 2
        dddelta = dddelta' @> 0
        
        c' =(x + _ZT*e31)**2/2 + (y + _ZT*e32)**2/2 + (z + _ZT*e33)**2/2 - r**2/2
        
        cdot' =dx*(x + _ZT*e31) + dy*(y + _ZT*e32) + dz*(z + _ZT*e33) + _ZT*(w2 - ddelta*e23)*(e11*x + e12*y + e13*z + _ZT*e11*e31 + _ZT*e12*e32 + _ZT*e13*e33) - _ZT*(w1 - ddelta*e13)*(e21*x + e22*y + e23*z + _ZT*e21*e31 + _ZT*e22*e32 + _ZT*e23*e33)
        cddot' =(_ZT*conj(w1)*(e11*x + e12*y + e13*z + _ZT*e11*e31 + _ZT*e12*e32 + _ZT*e13*e33) + _ZT*conj(w2)*(e21*x + e22*y + e23*z + _ZT*e21*e31 + _ZT*e22*e32 + _ZT*e23*e33))*(w3 - ddelta*e33) + dx*(dx + _ZT*e11*w2 - _ZT*e21*w1 - _ZT*ddelta*e11*e23 + _ZT*ddelta*e13*e21) + dy*(dy + _ZT*e12*w2 - _ZT*e22*w1 - _ZT*ddelta*e12*e23 + _ZT*ddelta*e13*e22) + dz*(dz + _ZT*e13*w2 - _ZT*e23*w1) + ddx*(x + _ZT*e31) + ddy*(y + _ZT*e32) + ddz*(z + _ZT*e33) - (w1 - ddelta*e13)*(e21*(_ZT*dx - _ZT**2*e21*(conj(w1) - ddelta*e13) + _ZT**2*e11*(conj(w2) - ddelta*e23)) + e22*(_ZT*dy - _ZT**2*e22*(conj(w1) - ddelta*e13) + _ZT**2*e12*(conj(w2) - ddelta*e23)) + _ZT*e33*(z*conj(w1) + ddelta*e11*x + ddelta*e12*y + _ZT*e33*conj(w1) + _ZT*ddelta*e11*e31 + _ZT*ddelta*e12*e32) + _ZT*e23*(dz + _ZT*e13*conj(w2) - _ZT*e23*conj(w1)) + _ZT*e31*(conj(w1) - ddelta*e13)*(x + _ZT*e31) + _ZT*e32*(conj(w1) - ddelta*e13)*(y + _ZT*e32)) + (w2 - ddelta*e23)*(e11*(_ZT*dx - _ZT**2*e21*(conj(w1) - ddelta*e13) + _ZT**2*e11*(conj(w2) - ddelta*e23)) + e12*(_ZT*dy - _ZT**2*e22*(conj(w1) - ddelta*e13) + _ZT**2*e12*(conj(w2) - ddelta*e23)) - _ZT*e33*(z*conj(w2) + ddelta*e21*x + ddelta*e22*y + _ZT*e33*conj(w2) + _ZT*ddelta*e21*e31 + _ZT*ddelta*e22*e32) + _ZT*e13*(dz + _ZT*e13*conj(w2) - _ZT*e23*conj(w1)) - _ZT*e31*(conj(w2) - ddelta*e23)*(x + _ZT*e31) - _ZT*e32*(conj(w2) - ddelta*e23)*(y + _ZT*e32)) + _ZT*(dw2 - dddelta*e23)*(e11*x + e12*y + e13*z + _ZT*e11*e31 + _ZT*e12*e32 + _ZT*e13*e33) - _ZT*(dw1 - dddelta*e13)*(e21*x + e22*y + e23*z + _ZT*e21*e31 + _ZT*e22*e32 + _ZT*e23*e33) - _ZT*dddelta*(e11*e23*x - e13*e21*x + e12*e23*y - e13*e22*y + _ZT*e11*e23*e31 - _ZT*e13*e21*e31 + _ZT*e12*e23*e32 - _ZT*e13*e22*e32)


main :: IO ()
main = do
  let x0 :: Vector Double
      x0 = fromList [ 1.154244772411
                    , -0.103540608242
                    , -0.347959211327
                    , 0.124930983341
                    , 0.991534857363
                    , 0.035367725910
                    , 0.316039689643
                    , -0.073559821379
                    , 0.945889986864
                    , 0.940484536806
                    , -0.106993361072
                    , -0.322554269411
                    , 0.000000000000
                    , 0.000000000000
                    , 0.000000000000
                    , 0.137035790811
                    , 3.664945343102
                    , -1.249768772258
                    , 0.000000000000
                    , 3.874600000000
                    ]
      u = fromList [0,0,0]
      r = 1.2

  let xdot _t x = fst $ modelInteg r x u
      h0 = 1e-12
      abstol = 1e-6
      reltol = 1e-4
      tf = 20
      ts = linspace 300 (0,tf)
      sol = odeSolveV RKf45 h0 abstol reltol xdot x0 ts
      constraints = map (snd . flip (modelInteg r) u) (toRows sol)

  plotPaths [] $ map (zip (toList ts)) (map toList (toColumns sol))
  plotPaths [] $ map (zip (toList ts)) (map toList (toColumns (fromRows constraints)))
