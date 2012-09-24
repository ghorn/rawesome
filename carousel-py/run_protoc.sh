#!/bin/sh

protoc -I.. --python_out=. ../kite.proto
