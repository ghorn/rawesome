{-# OPTIONS_GHC -Wall #-}
--{-# OPTIONS_GHC -ddump-splices #-}
{-# Language TemplateHaskell #-}
{-# Language DeriveDataTypeable #-}

module TestQuotes where

import Control.Concurrent

import Quotes --( f, MyType(..) )
import Data.Typeable
import Data.Data
import qualified Text.ProtocolBuffers.Header as P'

--import Kite.TestMessages ( TestMessages )
--import qualified Kite.CarouselState as CS
--import qualified Kite.Dcm as Dcm
--import qualified Kite.Xyz as KiteXyz

data Xyz = MkXyz { x_ :: Double
                 , y_ :: Double
                 , z_ :: Double
                 }
data Axyz = MkAxyz { a_ :: Double
                   , xyz_ :: Xyz
                   }

data TestMessages = TestMessages
                    { a_double :: !P'.Double
--                    , a_float :: !P'.Float
--                    , a_int32 :: !P'.Int32
--                    , a_int64 :: !P'.Int64
--                    , a_uint32 :: !P'.Word32
--                    , a_uint64 :: !P'.Word64
--                    , a_sint32 :: !P'.Int32
--                    , a_sint64 :: !P'.Int64
--                    , a_fixed32 :: !P'.Word32
--                    , a_fixed64 :: !P'.Word64
--                    , a_sfixed32 :: !P'.Int32
--                    , a_sfixed64 :: !P'.Int64
--                    , a_bool :: !P'.Bool
--                    , a_string :: !P'.Utf8
--                    , a_bytes :: !P'.ByteString
                    } deriving ( Show, Eq, Ord, Typeable, Data )

anAxyz :: Axyz
anAxyz = MkAxyz 7 (MkXyz 1 2 3)

increment :: Axyz -> Axyz
increment (MkAxyz a (MkXyz x y z)) = MkAxyz (a+1) (MkXyz (x + 0.1) (y + 0.2) (z + 0.3))

printVarInfo :: VarInfo -> IO ()
printVarInfo (VarInfo name type' mv) = do
  vals <- readMVar mv
  putStrLn $ "(" ++ show type' ++ ")  "++ name ++ ": " ++ show vals

go :: IO ()
go = do
  (receiveNewMessage, infos) <- $(setupTelem "position." ''Axyz)

  let printLog = mapM_ printVarInfo infos

      updateLoop 0 _ = return ()
      updateLoop n anAxyz' = do
        receiveNewMessage anAxyz'
        putStrLn ""
        printLog
        updateLoop (n-1::Int) (increment anAxyz')

  printLog
  updateLoop 4 anAxyz
  --woo <- $(f ''KiteXyz.Xyz)
  --woo <- $(f ''Dcm.Dcm)
  --woo <- $(f ''CS.CarouselState)

woo :: IO ()
woo = do
  (_receiveNewMessage, infos) <- $(setupTelem "position." ''TestMessages)
  mapM_ printVarInfo infos
