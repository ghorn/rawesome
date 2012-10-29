#!/bin/sh

mkdir -p carousel-cpp
protoc --python_out=rawe --cpp_out=carousel-cpp kite.proto

if [ -x hprotoc ]
then
	(
	cd kitevis
	hprotoc -I.. --haskell_out=src kite.proto
	)
	(
	cd trajectory-gui
	hprotoc -I.. --haskell_out=src kite.proto
	)
else
	echo "Not generating haskell protobuf interface since hprotoc not installed;"
	echo 'Try running "cabal install hprotoc"'
fi
