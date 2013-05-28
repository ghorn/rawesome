{-# OPTIONS_GHC -Wall #-}
{-# Language DoAndIfThenElse #-}
{-# Language OverloadedStrings #-}
{-# Language CPP #-}

module Main ( main ) where

import Data.Foldable ( toList )
#if OSX
import qualified System.ZMQ3 as ZMQ
#else
import qualified System.ZMQ as ZMQ
#endif
import Control.Concurrent ( MVar, forkIO, modifyMVar_, newMVar, readMVar)
import Control.Monad ( forever )
import qualified Data.ByteString.Lazy as BL
import Text.ProtocolBuffers ( messageGet )
import Text.ProtocolBuffers.Basic ( uToString )
--import System.Remote.Monitoring ( forkServer )

import qualified Kite.PendulumOpt as PO

import SpatialMath
import Vis
import DrawAC

data State = State { sTrails :: [Xyz Double]
                   , sMessages :: [String]
                   }


drawFun :: Maybe State -> VisObject Double
----drawFun state = VisObjects $ [axes] ++ (map text [-5..5]) ++ [boxText, ac, plane,trailLines]
drawFun Nothing = VisObjects []
drawFun (Just state) = VisObjects $ [txt, plane, trailLines, yLines, zLines]
  where
--    axes = Trans (Xyz 0 (-0.1) 0) $ Axes (0.5, 15)
    plane = Trans (Xyz 0 0 1) $ Plane (Xyz 0 0 1) (makeColor 1 1 1 1) (makeColor 0.2 0.3 0.32 1)
    txt = VisObjects $ zipWith (\s k -> Text2d s (30,30*k) TimesRoman24 (makeColor 1 1 1 1)) (sMessages state) [1..]
    trailLines = drawTrail (sTrails state) (\a -> makeColor a (1-a) 0 1)

    yLines = VisObjects $ map mkLine (sTrails state)
      where
        mkLine (Xyz x y z) = Line [Xyz x y z, Xyz x 0 z] (makeColor 0 0.2 1 0.8)

    zLines = VisObjects $ map mkLine (sTrails state)
      where
        mkLine (Xyz x y z) = Line [Xyz x y z, Xyz x y 0] (makeColor 1 0.2 0 0.8)


updateState :: PO.PendulumOpt -> (Maybe State) -> IO (Maybe State)
updateState po _ =
  return $ Just $ State { sMessages = reverse $ map uToString $ toList $ PO.messages po
                        , sTrails = zipWith3 Xyz xs ys zs
                        }
  where
    xs = toList $ PO.x po
    zs = toList $ PO.z po
    ys = take len [0,dy..]
      where
        len = length xs
        dy = 0.3 / realToFrac len
    
withContext :: (ZMQ.Context -> IO a) -> IO a
#if OSX
withContext = ZMQ.withContext
#else
withContext = ZMQ.withContext 1
#endif

sub :: MVar (Maybe State) -> IO ()
sub m = withContext $ \context -> do
#if OSX
  let receive = ZMQ.receive
#else
  let receive = flip ZMQ.receive []
#endif
  ZMQ.withSocket context ZMQ.Sub $ \subscriber -> do
    ZMQ.connect subscriber "tcp://localhost:5563"
    ZMQ.subscribe subscriber "pendulum-opt"
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
  m <- newMVar Nothing
  _ <- forkIO (sub m)
  
  let simFun _ _ = return ()
      df _ = fmap drawFun (readMVar m)
  simulateIO (Just ((1260,940),(1930,40))) "pendulum optimization" ts () df simFun
