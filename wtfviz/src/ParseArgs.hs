{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE DeriveDataTypeable #-}

module ParseArgs ( getip
                 ) where

import System.Console.CmdArgs ( (&=), Data, Typeable )
import qualified System.Console.CmdArgs as CA

data VisArgs = VisArgs { ipfile :: String
                       , ip :: String
                       , followkite :: Bool
                       } deriving (Show, Data, Typeable)

myargs :: VisArgs
myargs = VisArgs { ipfile = "" &= CA.help "file to read IP address out of" &= CA.typ "FILENAME"
                 , ip = ""     &= CA.help "an IP address" &= CA.typ "ADDRESS"
                 , followkite = False &= CA.help "rotate the camera to follow the current estimate"
                 } &= CA.summary "the kite visualizer program"

getip :: String -> String -> IO (String,Bool)
getip programname defaultip = do
  a <- CA.cmdArgs (myargs &= CA.program programname)
  ip' <- case (ipfile a,ip a) of
    ("","") -> return defaultip
    ("",x) -> return x
    (f,"") -> fmap (head . lines) (readFile f)
    (_,_) -> error "please only specify your ip address one way"
  return (ip', followkite a)
--  
--main :: IO ()
--main = do
--  ip' <- getip "defaultip"
--  print ip'
--  putStrLn "finished successfully"
