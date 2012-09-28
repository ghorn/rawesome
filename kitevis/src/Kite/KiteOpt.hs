{-# LANGUAGE BangPatterns, DeriveDataTypeable, FlexibleInstances, MultiParamTypeClasses #-}
module Kite.KiteOpt (KiteOpt(..)) where
import Prelude ((+), (/))
import qualified Prelude as Prelude'
import qualified Data.Typeable as Prelude'
import qualified Data.Data as Prelude'
import qualified Text.ProtocolBuffers.Header as P'
import qualified Kite.CarouselState as Kite (CarouselState)
 
data KiteOpt = KiteOpt{css :: !(P'.Seq Kite.CarouselState), endTime :: !P'.Double, iters :: !P'.Int32, wind_x :: !P'.Double}
             deriving (Prelude'.Show, Prelude'.Eq, Prelude'.Ord, Prelude'.Typeable, Prelude'.Data)
 
instance P'.Mergeable KiteOpt where
  mergeAppend (KiteOpt x'1 x'2 x'3 x'4) (KiteOpt y'1 y'2 y'3 y'4)
   = KiteOpt (P'.mergeAppend x'1 y'1) (P'.mergeAppend x'2 y'2) (P'.mergeAppend x'3 y'3) (P'.mergeAppend x'4 y'4)
 
instance P'.Default KiteOpt where
  defaultValue = KiteOpt P'.defaultValue P'.defaultValue P'.defaultValue P'.defaultValue
 
instance P'.Wire KiteOpt where
  wireSize ft' self'@(KiteOpt x'1 x'2 x'3 x'4)
   = case ft' of
       10 -> calc'Size
       11 -> P'.prependMessageSize calc'Size
       _ -> P'.wireSizeErr ft' self'
    where
        calc'Size = (P'.wireSizeRep 1 11 x'1 + P'.wireSizeReq 1 1 x'2 + P'.wireSizeReq 1 5 x'3 + P'.wireSizeReq 1 1 x'4)
  wirePut ft' self'@(KiteOpt x'1 x'2 x'3 x'4)
   = case ft' of
       10 -> put'Fields
       11 -> do
               P'.putSize (P'.wireSize 10 self')
               put'Fields
       _ -> P'.wirePutErr ft' self'
    where
        put'Fields
         = do
             P'.wirePutRep 10 11 x'1
             P'.wirePutReq 17 1 x'2
             P'.wirePutReq 24 5 x'3
             P'.wirePutReq 33 1 x'4
  wireGet ft'
   = case ft' of
       10 -> P'.getBareMessageWith update'Self
       11 -> P'.getMessageWith update'Self
       _ -> P'.wireGetErr ft'
    where
        update'Self wire'Tag old'Self
         = case wire'Tag of
             10 -> Prelude'.fmap (\ !new'Field -> old'Self{css = P'.append (css old'Self) new'Field}) (P'.wireGet 11)
             17 -> Prelude'.fmap (\ !new'Field -> old'Self{endTime = new'Field}) (P'.wireGet 1)
             24 -> Prelude'.fmap (\ !new'Field -> old'Self{iters = new'Field}) (P'.wireGet 5)
             33 -> Prelude'.fmap (\ !new'Field -> old'Self{wind_x = new'Field}) (P'.wireGet 1)
             _ -> let (field'Number, wire'Type) = P'.splitWireTag wire'Tag in P'.unknown field'Number wire'Type old'Self
 
instance P'.MessageAPI msg' (msg' -> KiteOpt) KiteOpt where
  getVal m' f' = f' m'
 
instance P'.GPB KiteOpt
 
instance P'.ReflectDescriptor KiteOpt where
  getMessageInfo _ = P'.GetMessageInfo (P'.fromDistinctAscList [17, 24, 33]) (P'.fromDistinctAscList [10, 17, 24, 33])
  reflectDescriptorInfo _
   = Prelude'.read
      "DescriptorInfo {descName = ProtoName {protobufName = FIName \".kite.KiteOpt\", haskellPrefix = [], parentModule = [MName \"Kite\"], baseName = MName \"KiteOpt\"}, descFilePath = [\"Kite\",\"KiteOpt.hs\"], isGroup = False, fields = fromList [FieldInfo {fieldName = ProtoFName {protobufName' = FIName \".kite.KiteOpt.css\", haskellPrefix' = [], parentModule' = [MName \"Kite\",MName \"KiteOpt\"], baseName' = FName \"css\"}, fieldNumber = FieldId {getFieldId = 1}, wireTag = WireTag {getWireTag = 10}, packedTag = Nothing, wireTagLength = 1, isPacked = False, isRequired = False, canRepeat = True, mightPack = False, typeCode = FieldType {getFieldType = 11}, typeName = Just (ProtoName {protobufName = FIName \".kite.CarouselState\", haskellPrefix = [], parentModule = [MName \"Kite\"], baseName = MName \"CarouselState\"}), hsRawDefault = Nothing, hsDefault = Nothing},FieldInfo {fieldName = ProtoFName {protobufName' = FIName \".kite.KiteOpt.endTime\", haskellPrefix' = [], parentModule' = [MName \"Kite\",MName \"KiteOpt\"], baseName' = FName \"endTime\"}, fieldNumber = FieldId {getFieldId = 2}, wireTag = WireTag {getWireTag = 17}, packedTag = Nothing, wireTagLength = 1, isPacked = False, isRequired = True, canRepeat = False, mightPack = False, typeCode = FieldType {getFieldType = 1}, typeName = Nothing, hsRawDefault = Nothing, hsDefault = Nothing},FieldInfo {fieldName = ProtoFName {protobufName' = FIName \".kite.KiteOpt.iters\", haskellPrefix' = [], parentModule' = [MName \"Kite\",MName \"KiteOpt\"], baseName' = FName \"iters\"}, fieldNumber = FieldId {getFieldId = 3}, wireTag = WireTag {getWireTag = 24}, packedTag = Nothing, wireTagLength = 1, isPacked = False, isRequired = True, canRepeat = False, mightPack = False, typeCode = FieldType {getFieldType = 5}, typeName = Nothing, hsRawDefault = Nothing, hsDefault = Nothing},FieldInfo {fieldName = ProtoFName {protobufName' = FIName \".kite.KiteOpt.wind_x\", haskellPrefix' = [], parentModule' = [MName \"Kite\",MName \"KiteOpt\"], baseName' = FName \"wind_x\"}, fieldNumber = FieldId {getFieldId = 4}, wireTag = WireTag {getWireTag = 33}, packedTag = Nothing, wireTagLength = 1, isPacked = False, isRequired = True, canRepeat = False, mightPack = False, typeCode = FieldType {getFieldType = 1}, typeName = Nothing, hsRawDefault = Nothing, hsDefault = Nothing}], keys = fromList [], extRanges = [], knownKeys = fromList [], storeUnknown = False, lazyFields = False}"