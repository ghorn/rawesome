{-# OPTIONS_GHC -Wall #-}
{-# Language TemplateHaskell #-}

module Quotes where

import Control.Concurrent ( MVar, newMVar, modifyMVar_ )
import Language.Haskell.TH
import Language.Haskell.TH.Syntax ( Lift(..) )
import Data.Sequence ( Seq, (|>) )
import qualified Data.Sequence as S
--import qualified GHC.Types

-- keep this abstract so that we can use a Seq or Vector later
type ContainerType a = Seq a
emptyContainer :: ContainerType a
emptyContainer = S.empty

appendContainer :: a -> ContainerType a -> ContainerType a
appendContainer x xs = xs |> x

data MyType = MyDouble
            | MyString
              deriving Show
instance Lift MyType where
  lift MyDouble = [| MyDouble |]
  lift MyString = [| MyString |]

data VarInfo = VarInfo String MyType (MVar (ContainerType Double))

data Output = Output { outputName :: Name
                     , outputMVarName :: Name
                     , outputString :: String
                     , outputMyType :: MyType
                     }
instance Show Output where
  show (Output a b c d) = "Output " ++ show a ++ " " ++ show b ++ " " ++ show c ++ " " ++ show d

-- | take a record syntax field and return usable stuff
handleField :: String -> (Name, Type) -> Q (Pat, [Output])
handleField prefix (name, ConT type') = do
  (TyConI (DataD _ _dataName _ [constructor] _)) <- reify type'
--  let msg = init $ unlines
--            [ "---------------- handlefield: -----------------"
--            , "    name: " ++ show name
--            , "    dataName: " ++ show _dataName
--            , "    constructor: " ++ show constructor
--            ]
--  reportWarning msg
  case constructor of
    -- recursive protobuf
    (RecC {}) -> handleConstructor (prefix ++ nameBase name ++ ".") constructor
    -- recursive protobuf
    (NormalC {}) -> case show type' of
      "GHC.Types.Double" -> do
        mn <- newName ("m_" ++ nameBase name)
        patternName <- newName (nameBase name)
        let pattern = VarP patternName
        return (pattern, [Output patternName mn (prefix ++ nameBase name) MyDouble])
      _ -> error $ "handleField: only handles Double prims now (" ++ show constructor ++ ")"
    _ -> error $ "handleField: unrecognized constructor " ++ show constructor
handleField _ x = error $ "handleField: the \"impossible\" happened" ++ show x


-- | Take a constructor with multiple fields, call handleFields on each of them,
--   assemble the result
handleConstructor :: String -> Con -> Q (Pat, [Output])
handleConstructor prefix (RecC conName varStrictTypes) = do
  let varTypes = map (\(x,_,z) -> (x,z)) varStrictTypes
  (fieldPatterns,outputs)  <- fmap unzip $ mapM (handleField prefix) varTypes
  let cOutputs = concat outputs
  let conPattern = ConP conName fieldPatterns
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
  -- get the info
  TyConI (DataD _ _typeName _ [constructor] _ ) <- reify typ
  (pattern, outputs) <- handleConstructor prefix constructor
  -- what to do with these outputs
  let outMVs = map outputMVarName outputs
  let outStrings = map outputString outputs
  let outMyTypes = map outputMyType outputs
  let outNames = map outputName outputs

  newMVar' <- [| newMVar |]
  emptyContainer' <- [| emptyContainer |]
  let makeMVars = map (\mv -> BindS (VarP mv) $ AppE newMVar' emptyContainer') outMVs
  updateFunName <- newName ("update_" ++ nameBase _typeName)
  modmvar <- [| \mv x -> modifyMVar_ mv (return . (appendContainer x)) |]
  let defUpdate = LetS [FunD updateFunName [Clause [pattern] (NormalB updates) []]]
        where
          updates :: Exp
          updates = DoE $ zipWith gg outMVs outNames
            where
              gg :: Name -> Name -> Stmt
              gg mv x = NoBindS $ AppE (AppE modmvar (VarE mv)) (VarE x)

  retStuff <- [| return ( $(varE updateFunName)
                        , zipWith3 VarInfo outStrings outMyTypes $(listE (map varE outMVs))
                        ) |]
  return (DoE $ makeMVars ++ [defUpdate, NoBindS retStuff])
