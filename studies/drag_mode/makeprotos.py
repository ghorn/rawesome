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

from carousel_conf import conf
import rawe

if __name__=='__main__':
    autogenDir = 'autogen'
    topname = 'kite'
    dae = rawe.models.crosswind_drag_mode(conf)
    rawe.utils.mkprotos.writeAll(dae, topname, autogenDir,haskellDirs=['plot-ho-matic'])

