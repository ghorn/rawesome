{-# OPTIONS_GHC -Wall #-}

module Main where

import Prelude hiding (span)


--import Casadi.Dvda
--import Casadi.SXM
--import Casadi.DAE
import Casadi.Integrator
import Casadi.IdasOptions
import SpatialMath
import Vis

import Foreign.C.Types ( CDouble )
import Data.Vector.Storable ( Vector, (!) )
import qualified Data.Vector.Storable as V
import Numeric.Container hiding ( Vector )
import Text.Printf ( printf )

import Draw
import Joy
import Model


data State = State { sTrails :: [[Xyz Double]]
                   , sX :: Vector CDouble
                   , sU :: (CDouble, CDouble)
                   }

toNice :: Vector CDouble -> (Xyz Double, Quat Double, Xyz Double, Xyz Double)
toNice state' = (xyz, q'n'b, r'n0'a0, r'n0't0)
  where
    state = V.map realToFrac state'
    e11 = state ! 3
    e12 = state ! 4
    e13 = state ! 5

    e21 = state ! 6
    e22 = state ! 7
    e23 = state ! 8

    e31 = state ! 9
    e32 = state ! 10
    e33 = state ! 11

    delta = state ! 18

    q'nwu'ned = Quat 0 1 0 0

    q'n'a = Quat (cos(0.5*delta)) 0 0 (sin(-0.5*delta))

    q'aNWU'bNWU = quatOfDcmB2A $ fromLists [ [e11, e21, e31]
                                           , [e12, e22, e32]
                                           , [e13, e23, e33]
                                           ]
    q'a'b = q'nwu'ned * q'aNWU'bNWU * q'nwu'ned
    q'n'b = q'n'a * q'a'b
    q'n'aNWU = q'n'a * q'nwu'ned

    rArm = Xyz 1.085 0 0
    xyzArm = rArm + Xyz (state @> 0) (state @> 1) (state @> 2)
    xyz = rotVecByQuatB2A q'n'aNWU xyzArm

    zt = -0.01
    r'n0'a0 = rotVecByQuatB2A q'n'a rArm
    r'n0't0 = xyz + (rotVecByQuatB2A q'n'b $ Xyz 0 0 (-zt))

drawFun :: State -> VisObject Double
drawFun state = VisObjects $ [axes, txt, ac, plane, trailLines, arm, line]
  where
    (pos@(Xyz px py pz), quat, r'n0'a0, r'n0't0) = toNice (sX state)
    
    axes = VisAxes (0.5, 15) (Xyz 0 0 0) (Quat 1 0 0 0)
    arm  = VisLine [Xyz 0 0 0, r'n0'a0] $ makeColor 1 1 0 1
    line = VisLine [r'n0'a0, r'n0't0]   $ makeColor 0 1 1 1
    plane = VisPlane (Xyz 0 0 1) 2 (makeColor 1 1 1 1) (makeColor 0.2 0.3 0.32 1)
    ddelta = realToFrac $ (sX state) @> 19 :: Double
--    [c,cdot,cddot] = toList $ snd $ modelInteg 1.2 (sX state) (fromList [tc0,0,0])

    (u0',u1') = sU state
    u0 = realToFrac u0' :: Double
    u1 = realToFrac u1' :: Double
    txt = VisObjects
          [ Vis2dText (printf "x: %.3f" px) (30,90) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "y: %.3f" py) (30,60) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "z: %.3f" pz) (30,30) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "RPM: %.3f" (ddelta*60/(2*pi))) (30,120) TimesRoman24 (makeColor 1 1 1 1)
--          , Vis2dText (printf "c:   %.3g" c    ) (30,150) TimesRoman24 (makeColor 1 1 1 1)
--          , Vis2dText (printf "c':  %.3g" cdot ) (30,180) TimesRoman24 (makeColor 1 1 1 1)
--          , Vis2dText (printf "c'': %.3g" cddot) (30,210) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "u0: %.3g \t(*180/pi = %.3f)" u0 (u0*180/pi)) (30,260) TimesRoman24 (makeColor 1 1 1 1)
          , Vis2dText (printf "u1: %.3g \t(*180/pi = %.3f)" u1 (u1*180/pi)) (30,240) TimesRoman24 (makeColor 1 1 1 1)
          ]
    (ac,_) = drawAc pos quat

    trailLines = drawTrails (sTrails state)

--tc0 :: Double
tc0 = 2*389.970797939731

main :: IO ()
main = do
  let --x0 :: Vector Double
      x0 = V.fromList [ 1.154244772411
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

      r = 1.2 :: Double

  let ts = 0.02

  putStrLn "creating model"
  (daeIn, daeOut, c, cdot) <- sxModel

  putStrLn "creating integrator"
  let integratorOptions = [ Reltol 1e-6
                          , Abstol 1e-8
                          , T0 0
                          , Tf (ts)
                          ]
  int <- createIdasIntegrator daeIn daeOut [] integratorOptions

  let updateTrail :: [Xyz a] -> Xyz a -> [Xyz a]
      updateTrail trail0 trail
        | length trail0 < 65 = trail:trail0
        | otherwise = take 65 (trail:trail0)

  js <- setupJoy
  let --simFun :: Float -> State -> IO State
      simFun _ (State {sTrails = trails0, sX = x'}) = do
        j <- getJoy js
        let (u0':u1':_:_:_:_:throttle:_) = jsAxes j
            u0 = -u0'*0.02
            u1 = u1'*0.05
            tc = 600*(1 - throttle)
        x <- int x' (Just (V.fromList [tc,u0,u1]))
        let (_,trails) = drawAc pos q
            (pos,q,_,_) = toNice x
        return $
          State { sTrails = zipWith updateTrail trails0 (fmap (fmap realToFrac) trails)
                , sX = x
                , sU = (u0,u1)
                }

  let state0 = State { sX = x0
                     , sTrails = [[],[],[]]
                     , sU = (0,0)
                     }
  simulateIO ts state0 (return . drawFun) simFun

--  print sol
--  plotPaths [] $ map (zip (toList ts)) (map toList (toColumns sol))
--  plotPaths [] $ map (zip (toList ts)) (map toList (toColumns (fromRows constraints)))
