# Copyright 2012-2013 Greg Horn
#
# This file is part of rawesome.
#
# rawesome is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rawesome is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with rawesome.  If not, see <http://www.gnu.org/licenses/>.

topname = rocket
ag = autogen
genscript = make_rocket_protos.py
daefile = rocket_dae.py

autogenTargets = \
	$(ag)/__init__.py \
	$(ag)/protoConverters.cpp \
	$(ag)/protoConverters.h \
	$(ag)/$(topname)_dimensions.h \
	$(ag)/$(topname).proto \
	$(ag)/$(topname)_structs.h \
	$(ag)/to$(topname)Proto.py

plotter = plotter/dist/build/plotter/plotter

hstopname = $(shell echo $(topname)  | python -c "print raw_input().capitalize()")

ALL : $(autogenTargets) $(ag)/$(topname)_pb2.py $(ag)/$(hstopname).hs $(plotter)

$(autogenTargets) : $(genscript) $(daefile)
	python -B $(genscript) $(topname) $(ag)

$(ag)/$(topname)_pb2.py : $(ag)/$(topname).proto
	protoc --python_out=. $<

$(ag)/$(hstopname).hs : $(ag)/$(topname).proto
	hprotoc --haskell_out=$(ag) $<

$(plotter) : $(ag)/$(hstopname).hs
	cd plotter && cabal configure && cabal build

clean :
	rm -rf autogen
	cd plotter && cabal clean
