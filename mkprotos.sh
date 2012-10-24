#!/bin/sh

mkdir -p carousel-cpp
protoc --python_out=rawe --cpp_out=carousel-cpp kite.proto

cd kitevis
hprotoc -I.. --haskell_out=src kite.proto
cd ..

