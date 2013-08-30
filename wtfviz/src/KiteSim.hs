-- Copyright 2012-2013 Greg Horn
--
-- This file is part of rawesome.
--
-- rawesome is free software: you can redistribute it and/or modify
-- it under the terms of the GNU Lesser General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
--
-- rawesome is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU Lesser General Public License for more details.
--
-- You should have received a copy of the GNU Lesser General Public License
-- along with rawesome.  If not, see <http://www.gnu.org/licenses/>.

{-# OPTIONS_GHC -Wall #-}
{-# Language DoAndIfThenElse #-}
{-# Language OverloadedStrings #-}
{-# Language CPP #-}

module Main ( main ) where

import Data.Foldable ( toList )
import Data.Maybe ( fromMaybe )
import System.Random ( randomRs, mkStdGen)
#if OSX
import qualified System.ZMQ3 as ZMQ
#else
import qualified System.ZMQ as ZMQ
#endif
import Control.Concurrent ( MVar, forkIO, modifyMVar_, newMVar, readMVar)
import Control.Monad ( forever )
import qualified Data.ByteString.Lazy as BL
import Data.Packed ( fromLists )
import Text.ProtocolBuffers ( messageGet )
import Text.ProtocolBuffers.Basic ( uToString )
--import System.Remote.Monitoring ( forkServer )

import qualified Kite.CarouselState as CS
import qualified Kite.Dcm as Dcm
import qualified Kite.Xyz as KiteXyz

import SpatialMath
import Vis
import DrawAC

data State = State { sTrails :: [[Xyz Double]]
                   , sCS :: Maybe CS.CarouselState
                   , sParticles :: [Xyz Double]
                   }

toNice :: CS.CarouselState -> (Xyz Double, Quat Double, Xyz Double, Xyz Double, Double)
toNice cs = (xyz, q'n'b, r'n0'a0, r'n0't0, fromMaybe 1 $ CS.visSpan cs)
  where
    x = KiteXyz.x $ CS.kiteXyz cs
    y = KiteXyz.y $ CS.kiteXyz cs
    z = KiteXyz.z $ CS.kiteXyz cs

    r11 = Dcm.r11 $ CS.kiteDcm cs
    r12 = Dcm.r12 $ CS.kiteDcm cs
    r13 = Dcm.r13 $ CS.kiteDcm cs

    r21 = Dcm.r21 $ CS.kiteDcm cs
    r22 = Dcm.r22 $ CS.kiteDcm cs
    r23 = Dcm.r23 $ CS.kiteDcm cs

    r31 = Dcm.r31 $ CS.kiteDcm cs
    r32 = Dcm.r32 $ CS.kiteDcm cs
    r33 = Dcm.r33 $ CS.kiteDcm cs

    delta = CS.delta cs

    q'nwu'ned = Quat 0 1 0 0

    q'n'a = Quat (cos(0.5*delta)) 0 0 (sin(-0.5*delta))

    q'aNWU'bNWU = quatOfDcm $ fromLists [ [r11, r12, r13]
                                        , [r21, r22, r23]
                                        , [r31, r32, r33]
                                        ]
    q'a'b = q'nwu'ned * q'aNWU'bNWU * q'nwu'ned
    q'n'b = q'n'a * q'a'b
    q'n'aNWU = q'n'a * q'nwu'ned

    rArm = Xyz (CS.rArm cs) 0 0
    xyzArm = rArm + Xyz x y z
    xyz = rotVecByQuatB2A q'n'aNWU xyzArm

    zt = CS.zt cs
    r'n0'a0 = rotVecByQuatB2A q'n'a rArm
    r'n0't0 = xyz + (rotVecByQuatB2A q'n'b $ Xyz 0 0 (-zt))

drawFun :: State -> VisObject Double
----drawFun state = VisObjects $ [axes] ++ (map text [-5..5]) ++ [boxText, ac, plane,trailLines]
drawFun (State {sCS=Nothing}) = VisObjects []
drawFun state@(State {sCS=Just cs}) =
  VisObjects [axes, txt, ac, plane, trailLines, arm, line, zLine, xyLine, points]
  where
    (pos@(Xyz x y z), quat, r'n0'a0, r'n0't0, visSpan) = toNice cs

    points = Points (sParticles state) (Just 2) $ makeColor 1 1 1 0.5
    zLine = Line [Xyz x y (planeZ-0.01), pos]            $ makeColor 0.1 0.2 1 0.5
    xyLine = Line [Xyz x y (planeZ-0.01), Xyz 0 0 (planeZ-0.01)] $ makeColor 0.2 0.7 1 0.5

    axes = Axes (0.5, 15)
    arm  = Line [Xyz 0 0 0, r'n0'a0] $ makeColor 1 1 0 lineAlpha
    line = Line [r'n0'a0, r'n0't0]   $ makeColor 0 1 1 lineAlpha
    plane = Trans (Xyz 0 0 planeZ) $ Plane (Xyz 0 0 1) (makeColor 1 1 1 1) (makeColor 0.2 0.3 0.32 (realToFrac planeAlpha))
    planeZ' = planeZ-0.5
    planeAlpha
      | z < planeZ' = 1
      | z < planeZ'+2 = (planeZ'+2-z)/2
      | otherwise = 0

    txt = VisObjects $
          zipWith (\s k -> Text2d (uToString s) (30,fromIntegral $ 30*k) TimesRoman24 (makeColor 1 1 1 1)) messages (reverse [1..length messages])
    messages = toList $ CS.messages cs

    ac = Trans pos $ Scale (visSpan,visSpan,visSpan) ac'
    (ac',_) = drawAc kiteAlpha (Xyz 0 0 0) quat
    lineAlpha = realToFrac $ fromMaybe 1 (CS.lineTransparency cs)
    kiteAlpha = realToFrac $ fromMaybe 1 (CS.kiteTransparency cs)

    trailLines = drawTrails (sTrails state)

planeZ :: Double
planeZ = 1

particleBox :: Double
particleBox = 8

state0 :: State
state0 = State { sCS = Nothing
               , sTrails = [[],[],[]]
               , sParticles = take 300 $ randomRs (Xyz (-particleBox) (-particleBox) (planeZ-2*particleBox),
                                                   Xyz particleBox particleBox planeZ)
                              (mkStdGen 0)
               }

updateTrail :: [Xyz a] -> Xyz a -> [Xyz a]
updateTrail trail0 xyz
  | length trail0 < 65 = xyz:trail0
  | otherwise = take 65 (xyz:trail0)

boundParticle :: Xyz Double -> Xyz Double
boundParticle xyz@(Xyz x y z)
  | x >  particleBox = boundParticle (Xyz (x-2*particleBox) y z)
  | x < -particleBox = boundParticle (Xyz (x+2*particleBox) y z)
  | y >  particleBox = boundParticle (Xyz x (y-2*particleBox) z)
  | y < -particleBox = boundParticle (Xyz x (y+2*particleBox) z)
  | z > planeZ               = boundParticle (Xyz x y (z-2*particleBox))
  | z < planeZ-2*particleBox = boundParticle (Xyz x y (z+2*particleBox))
  | otherwise = xyz

windShear :: Double -> Double -> Double
windShear w0 z
  | z' < zt = 0
  | otherwise = w0*log(z'/zt)/log(z0/zt)
  where
    z' = z + planeZ + zt + 2
    z0 = 100
    zt = 0.1

updateState :: CS.CarouselState -> State -> IO State
updateState cs x0 =
  return $ State { sCS = Just cs
                 , sTrails = zipWith updateTrail trails0 trails
                 , sParticles = map (\xyz@(Xyz _ _ z) -> boundParticle $ (Xyz (ts*(windShear w0 (-z))) 0 0) + xyz) (sParticles x0)
                 }
  where
    w0 = fromMaybe 0 (CS.w0 cs)
    trails0 = sTrails x0
    (pos,q,_,_,_) = toNice cs
    (_,trails) = drawAc 1 pos q

withContext :: (ZMQ.Context -> IO a) -> IO a
#if OSX
withContext = ZMQ.withContext
#else
withContext = ZMQ.withContext 1
#endif

sub :: MVar State -> IO ()
sub m = withContext $ \context -> do
#if OSX
  let receive = ZMQ.receive
#else
  let receive = flip ZMQ.receive []
#endif
  ZMQ.withSocket context ZMQ.Sub $ \subscriber -> do
    ZMQ.connect subscriber "tcp://localhost:5563"
    ZMQ.subscribe subscriber "carousel"
    forever $ do
      _ <- receive subscriber
      mre <- ZMQ.moreToReceive subscriber
      if mre
      then do
        msg <- receive subscriber
        let cs = case messageGet (BL.fromChunks [msg]) of
              Left err -> error err
              Right (cs',_) -> cs'
        modifyMVar_ m (updateState cs)
      else return ()

ts :: Double
ts = 0.02

main :: IO ()
main = do
--  _ <- forkServer "localhost" 8000
  m <- newMVar state0
  _ <- forkIO (sub m)

--  threadDelay 5000000
  let simFun _ _ = return ()
      df _ = fmap drawFun (readMVar m)
  simulateIO (Just ((1260,940),(1930,40))) "kite sim" ts () df simFun

