#!/bin/sh

cd carousel-py
protoc -I.. --python_out=. ../kite.proto
cd ..

cd kitevis
hprotoc -I.. --haskell_out=src kite.proto
cd ..
