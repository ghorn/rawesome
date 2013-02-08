{-# OPTIONS_GHC -Wall #-}
{-# Language MultiWayIf #-}
{-# Language TemplateHaskell #-}

module Quotes where

import Control.Concurrent ( MVar, newMVar, modifyMVar_ )
import Language.Haskell.TH
import Language.Haskell.TH.Syntax ( Lift(..) )
import Data.Sequence ( Seq, (|>) )
import qualified Data.Sequence as S

-- keep this abstract so that we can use a Seq or Vector later
type ContainerType a = Seq a
emptyContainer :: ContainerType a
emptyContainer = S.empty

appendContainer :: a -> ContainerType a -> ContainerType a
appendContainer x xs = xs |> x

data MyType = MyDouble
            | MyFloat
--            | MyInt32
--            | MyString
              deriving Show

instance Lift MyType where
  lift MyDouble = [| MyDouble |]
  lift MyFloat = [| MyFloat |]
--  lift MyInt32  = [| MyInt32 |]
--  lift MyString = [| MyString |]

data VarInfo = VarInfo String MyType (MVar (ContainerType Double))
--             | VarInfoF String MyType (MVar (ContainerType Float))

data Output = Output { outputName :: Name
                     , outputMVarName :: Name
                     , outputString :: String
                     , outputMyType :: MyType
                     }
instance Show Output where
  show (Output a b c d) = "Output " ++ show a ++ " " ++ show b ++ " " ++ show c ++ " " ++ show d

-- | take a record syntax field and return usable stuff
handleField :: String -> (Name, Type) -> Q (PatQ, [Output])
handleField prefix (name, ConT type') = do
  (TyConI (DataD _ _dataName _ [constructor] _)) <- reify type'
  let msg = init $ unlines
            [ "---------------- handlefield: -----------------"
            , "    name: " ++ show name
            , "    dataName: " ++ show _dataName
            , "    constructor: " ++ show constructor
            ]
--  reportWarning msg
  case constructor of
    -- recursive protobuf
    (RecC {}) -> handleConstructor (prefix ++ nameBase name ++ ".") constructor
    -- normal constructor
    (NormalC {}) -> do
      -- make the mvar name and pattern name
      mn <- newName ("m_" ++ nameBase name)
      patternName <- newName (nameBase name)
      let pattern = varP patternName

      -- lookup some type names
      (Just doubleName) <- lookupTypeName "Double"
      (Just floatName) <- lookupTypeName "Float"
--      (Just int32Name) <- lookupTypeName "Int32"
--      reportWarning $ "typeName: " ++ show int32Name

      if | type' == doubleName -> return (pattern, [Output patternName mn (prefix ++ nameBase name) MyDouble])
         | type' == floatName -> return (pattern, [Output patternName mn (prefix ++ nameBase name) MyFloat])
         | otherwise -> error $ "handleField: unhandled type (" ++ show constructor ++ ")"++"\n    "++msg
    _ -> error $ "handleField: unrecognized constructor " ++ show constructor
handleField _ x = error $ "handleField: the \"impossible\" happened" ++ show x


-- | Take a constructor with multiple fields, call handleFields on each of them,
--   assemble the result
handleConstructor :: String -> Con -> Q (PatQ, [Output])
handleConstructor prefix (RecC conName varStrictTypes) = do
  let varTypes = map (\(x,_,z) -> (x,z)) varStrictTypes
  (fieldPatterns,outputs)  <- fmap unzip $ mapM (handleField prefix) varTypes
  let cOutputs = concat outputs
  let conPattern = conP conName fieldPatterns
--  let msg = init $ unlines
--            [ "---------------- handleConstructor: ----------------------"
--            , "    conName: " ++ show conName ++ " ("++nameBase conName++")"
--            , "    varTypes: " ++ show varTypes
--            , "    conPattern: " ++ show conPattern
--            , "    outputs: " ++ show outputs
--            , "    cOutputs: " ++ show cOutputs
--            , "    ----------------------------------------------"
--            ]
--  reportWarning msg
  return (conPattern, cOutputs)
handleConstructor _ x = fail $ "\"" ++ show x ++ "\" is not a record syntax constructor"


setupTelem :: String -> Name -> Q Exp
setupTelem prefix typ = do
  -- get the type info
  TyConI (DataD _ _typeName _ [constructor] _ ) <- reify typ

  -- get the pattern and the names in a nice list of outputs
  (pattern, outputs) <- handleConstructor prefix constructor

  -- split the outputs
  let outMVs = map outputMVarName outputs
  let outStrings = map outputString outputs
  let outMyTypes = map outputMyType outputs
  let outNames = map outputName outputs

  -- define all the newMVars
  let makeMVars = map (\mv -> bindS (varP mv) [| newMVar emptyContainer |]) outMVs

  -- define the function to take new data and update the MVars
  updateFunName <- newName ("update_" ++ nameBase _typeName)
  let defUpdate = letS [funD updateFunName [clause [pattern] (normalB updates) []]]
        where
          updates :: ExpQ
          updates = doE $ zipWith gg outMVs outNames
            where
              gg :: Name -> Name -> StmtQ
              gg mv x = noBindS [| modifyMVar_ $(varE mv) (return . appendContainer $(varE x)) |]

  -- return (...)
  let retStuff = [| return ( $(varE updateFunName)
                           , zipWith3 VarInfo outStrings outMyTypes $(listE (map varE outMVs))
                           ) |]
  doE $ makeMVars ++ [defUpdate, noBindS retStuff]
