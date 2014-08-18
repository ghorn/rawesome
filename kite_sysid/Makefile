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

topname = carousel_sysid
ag = autogen
daefile = carousel_dae


hstopname = $(shell echo $(topname)  | python -c "print raw_input().capitalize()")

vis = viz/dist/build/multicarousel/multicarousel

all :: $(vis) $(autogenTargets) $(ag)/$(topname)_pb2.py $(ag)/$(hstopname).hs $(plotter)

$(vis) : $(ag)/$(hstopname).hs viz/*.hs
	cd viz && cabal configure -v0 && cabal build -v0

autogenTargets = \
	$(ag)/$(topname).proto \
#	$(ag)/__init__.py \
#	$(ag)/protoConverters.cpp \
#	$(ag)/protoConverters.h \
#	$(ag)/$(topname)_dimensions.h \
#	$(ag)/$(topname)_structs.h \
#	$(ag)/to$(topname)Proto.py

plotter = plotter/dist/build/plotter/plotter

$(autogenTargets) : $(daefile).py
	python -c "import rawe; import $(daefile); rawe.utils.mkprotobufs.writeAll($(daefile).makeDae(), '$(topname)', '$(ag)')"

$(ag)/$(topname)_pb2.py : $(ag)/$(topname).proto
	protoc --python_out=. $<

$(ag)/$(hstopname).hs : $(ag)/$(topname).proto
	hprotoc --haskell_out=$(ag) $< > /dev/null

$(plotter) : $(ag)/$(hstopname).hs plotter/*.hs
	cd plotter && cabal configure -v0 && cabal build -v0

clean ::
	cd viz && cabal clean -v0
	rm -rf autogen
#	cd plotter && cabal clean -v0
