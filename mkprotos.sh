#!/bin/sh

mkdir -p carousel-cpp
protoc --python_out=carousel-py --cpp_out=carousel-cpp kite.proto

cd kitevis
hprotoc -I.. --haskell_out=src kite.proto
cd ..

