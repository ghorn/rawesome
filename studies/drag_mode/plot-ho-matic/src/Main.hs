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
{-# Language CPP #-}
{-# Language DoAndIfThenElse #-}
{-# Language TemplateHaskell #-}
--{-# OPTIONS_GHC -ddump-splices #-}
--{-# Language OverloadedStrings #-}

module Main ( main ) where

#if OSX
import qualified System.ZMQ3 as ZMQ
#else
import qualified System.ZMQ as ZMQ
#endif
import qualified Control.Concurrent as CC
import Control.Monad ( forever )
import qualified Data.ByteString.Lazy as BL
import qualified Text.ProtocolBuffers as PB

--import qualified System.Remote.Monitoring as EKG

import qualified Kite.KiteTrajectory as KT

import ParseArgs ( getip )
import Plotter ( runPlotter, newChannel, makeAccessors )

main :: IO ()
main = do
--  ekgTid <- fmap EKG.serverThreadId $ EKG.forkServer "localhost" 8000
  ip <- getip "plot-ho-matic" "tcp://localhost:5563"
  putStrLn $ "using ip \""++ip++"\""

  let zeromqChannel = "kite-optimization"
  (c0, chan0) <- newChannel zeromqChannel $(makeAccessors ''KT.KiteTrajectory)
  listenerTid0 <- CC.forkIO (sub ip chan0 zeromqChannel)
  
  runPlotter [c0] [listenerTid0]

withContext :: (ZMQ.Context -> IO a) -> IO a
#if OSX
withContext = ZMQ.withContext
#else
withContext = ZMQ.withContext 1
#endif

sub :: (PB.Wire a, PB.ReflectDescriptor a) => String -> CC.Chan a -> String -> IO ()
sub ip chan name = withContext $ \context -> do
#if OSX
  let receive = ZMQ.receive
#else
  let receive = flip ZMQ.receive []
#endif
  ZMQ.withSocket context ZMQ.Sub $ \subscriber -> do
    ZMQ.connect subscriber ip
    ZMQ.subscribe subscriber name
    forever $ do
      _ <- receive subscriber
      mre <- ZMQ.moreToReceive subscriber
      if mre
      then do
        msg <- receive subscriber
        let cs = case PB.messageGet (BL.fromChunks [msg]) of
              Left err -> error err
              Right (cs',_) -> cs'
        CC.writeChan chan cs
      else return ()
