{-# OPTIONS_GHC -Wall #-}
{-# Language DoAndIfThenElse #-}
{-# Language CPP #-}

module Main where

import Data.Foldable ( toList )
import Data.Maybe ( isNothing )
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

import qualified MheMpc.MheMpcHorizons as MMH
import qualified MheMpc.Dae as Dae
import qualified MheMpc.DaePlus as DaePlus
import qualified MheMpc.DifferentialStates as DS

import SpatialMath
import qualified Xyz
import Vis
import DrawAC
import ParseArgs ( getip )

type State = Maybe MMH.MheMpcHorizons

data NiceKite = NiceKite { nk_xyz :: Xyz Double
                         , nk_q'n'b :: Quat Double
                         , nk_r'n0'a0 :: Xyz Double
                         , nk_r'n0't0 :: Xyz Double
                         , nk_lineAlpha :: Float
                         , nk_kiteAlpha :: Float
                         , nk_visSpan :: Double
                         }

toNice :: Double -> Double -> Double -> DaePlus.DaePlus -> NiceKite
toNice visSpan rArm' zt daeplus =
  NiceKite { nk_xyz = xyz
           , nk_q'n'b = q'n'b
           , nk_r'n0'a0 = r'n0'a0
           , nk_r'n0't0 = r'n0't0
           , nk_lineAlpha = realToFrac $ DaePlus.lineTransparency daeplus
           , nk_kiteAlpha = realToFrac $ DaePlus.kiteTransparency daeplus
           , nk_visSpan = visSpan
           }
  where
    ds = Dae.diffStates (DaePlus.dae daeplus)
    x = DS.x ds
    y = DS.y ds
    z = DS.z ds

    -- transpose going on here
    e11 = DS.e11 ds
    e12 = DS.e21 ds
    e13 = DS.e31 ds

    e21 = DS.e12 ds
    e22 = DS.e22 ds
    e23 = DS.e32 ds

    e31 = DS.e13 ds
    e32 = DS.e23 ds
    e33 = DS.e33 ds

    delta = DS.delta ds

    q'nwu'ned = Quat 0 1 0 0

    q'n'a = Quat (cos(0.5*delta)) 0 0 (sin(-0.5*delta))

    q'aNWU'bNWU = quatOfDcm $ fromLists [ [e11, e12, e13]
                                        , [e21, e22, e23]
                                        , [e31, e32, e33]
                                        ]
    q'a'b = q'nwu'ned * q'aNWU'bNWU * q'nwu'ned
    q'n'b = q'n'a * q'a'b
    q'n'aNWU = q'n'a * q'nwu'ned

    rArm = Xyz rArm' 0 0
    xyzArm = rArm + Xyz x y z
    xyz = rotVecByQuatB2A q'n'aNWU xyzArm

    r'n0'a0 = rotVecByQuatB2A q'n'a rArm
    r'n0't0 = xyz + (rotVecByQuatB2A q'n'b $ Xyz 0 0 (-zt))

drawOneKite :: Double -> NiceKite -> (VisObject Double, Double)
drawOneKite minLineLength niceKite
  | nk_lineAlpha niceKite > 0 = (VisObjects [ac, arm, line], z)
  | otherwise = (ac, z)
  where
    pos@(Xyz _ _ z) = nk_xyz niceKite
    quat = nk_q'n'b niceKite
    r'n0'a0 = nk_r'n0'a0 niceKite
    r'n0't0 = nk_r'n0't0 niceKite
    
    arm  = Line [Xyz 0 0 0, r'n0'a0] $ makeColor 1 1 0 (nk_lineAlpha niceKite)
    line = VisObjects [ line1 $ makeColor 1 0.2 0 (nk_lineAlpha niceKite)
                      , line2 $ makeColor 0 1 1 (nk_lineAlpha niceKite)
                      ]
      where
        line1 = Line [r'n0'a0, rMid] 
        line2 = Line [rMid, r'n0't0] 

        rMid = r'n0'a0 + fmap (* (1 - minLineLength/normDr)) dr
        dr = r'n0't0 - r'n0'a0
        normDr = Xyz.norm dr
        

    s = nk_visSpan niceKite
    ac = Trans pos $ Scale (s,s,s) ac'
    (ac',_) = drawAc (nk_kiteAlpha niceKite) (Xyz 0 0 0) quat

drawFun :: Bool -> State -> VisObject Double
drawFun _ Nothing = VisObjects []
drawFun followKite (Just mmh) = drawFun' kiteDelta mmh
  where
    kiteDelta =
      if followKite
      then Just $ DS.delta $ Dae.diffStates $ MMH.currentState mmh
      else Nothing

drawFun' :: Maybe Double -> MMH.MheMpcHorizons -> VisObject Double
drawFun' deltaRot mmh = cameraRot $ VisObjects $ [axes,txt] ++ maybeplane ++ [mheKites,Trans (Xyz 0 0 (-4)) mpcKites]
  where
    cameraRot = case deltaRot of
      Nothing -> id
      Just delta -> RotQuat $ Quat (cos (0.5*delta)) 0 0 (sin (0.5*delta))

    visSpan = MMH.visSpan mmh
    rArm = MMH.rArm mmh
    zt = MMH.zt mmh
    toNice' = toNice visSpan rArm zt
    mheKites = drawSomeKites $ map toNice' (toList (MMH.mheHorizon mmh))
    mpcKites = drawSomeKites $ map toNice' (toList (MMH.mpcHorizon mmh))
    
    axes = Axes (0.5, 15)
    maybeplane = if isNothing deltaRot then [plane] else []
    plane = Trans (Xyz 0 0 planeZ') $ Plane (Xyz 0 0 1) (makeColor 1 1 1 1) (makeColor 0.2 0.3 0.32 1)

    planeZ' = MMH.carouselArmHeight mmh
--    boxText = 3dText "I'm a plane" (Xyz 0 0 (x-0.2)) TimesRoman24 (makeColor 1 0 0 1)
    txt = VisObjects $
          zipWith (\s k -> Text2d (uToString s) (30,fromIntegral $ 30*k) TimesRoman24 (makeColor 1 1 1 1)) messages (reverse [1..length messages])
    messages = toList $ MMH.messages mmh

drawSomeKites :: [NiceKite] -> VisObject Double
drawSomeKites niceKites = VisObjects kitelist
  where
    (kitelist,_) = unzip $ map (drawOneKite minLineLength) niceKites
    minLineLength = minimum $ map lineLength niceKites
      where
        lineLength nk = Xyz.norm (nk_r'n0't0 nk - nk_r'n0'a0 nk)
    

particleBox :: Double
particleBox = 4

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
  | z >  particleBox = boundParticle (Xyz x y (z-2*particleBox))
  | z < -particleBox = boundParticle (Xyz x y (z+2*particleBox))
  | otherwise = xyz

updateState :: MMH.MheMpcHorizons -> State -> IO State
updateState ko _ = return $ Just ko

withContext :: (ZMQ.Context -> IO a) -> IO a
#if OSX
withContext = ZMQ.withContext
#else
withContext = ZMQ.withContext 1
#endif

sub :: String -> MVar State -> IO ()
sub ip m = withContext $ \context -> do
#if OSX
  let receive = ZMQ.receive
#else
  let receive = flip ZMQ.receive []
#endif
  ZMQ.withSocket context ZMQ.Sub $ \subscriber -> do
    ZMQ.connect subscriber ip
    ZMQ.subscribe subscriber "mhe-mpc-horizons"
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
  (ip,followkite) <- getip "multicarousel" "tcp://localhost:5563"
  putStrLn $ "using ip \""++ip++"\""

  m <- newMVar Nothing
  _ <- forkIO (sub ip m)
  
--  threadDelay 5000000
  let simFun _ _ = return ()
      df _ = fmap (drawFun followkite) (readMVar m)
  simulateIO (Just ((1260,940),(1930,40))) "multi-carousel" ts () df simFun
