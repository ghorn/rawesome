#
# This file is part of ACADO Toolkit.
#
# ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
# Copyright (C) 2008-2011 by Boris Houska and Hans Joachim Ferreau.
# All rights reserved.
#
# ACADO Toolkit is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# ACADO Toolkit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with ACADO Toolkit; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

################################################################################
#
# Description:
#	MHE export package configuration file
#
#
# Authors:
#	Milan Vukov, milan.vukov@esat.kuleuven.be
#
# Year:
#	2012.
#
# NOTE:
#	- This script is for Linux/Unix use only.
#
#	- PREREQUISITE: sourced mhe_export_env.sh in your ~/.bashrc file. This script
#		will try to find MHE export folders, libraries etc., but looking for them
#		in enviromental variables.
#
# Usage:
#	- Linux/Unix: TODO
#
################################################################################

################################################################################
#
# Search for package components
#
################################################################################

MESSAGE( STATUS "********************************************************************************" )
MESSAGE( STATUS "Looking for MHE export package: \n" )

#
# Include folders
#
MESSAGE( STATUS "Looking for MHE export include directories" )
SET( MHE_EXPORT_INCLUDE_DIRS $ENV{MHE_EXPORT_ENV_INCLUDE_DIRS} )
IF( MHE_EXPORT_INCLUDE_DIRS )
	MESSAGE( STATUS "Found MHE export include directories: ${MHE_EXPORT_INCLUDE_DIRS} \n" )
	SET( MHE_EXPORT_INCLUDE_DIRS_FOUND TRUE )
ELSE( MHE_EXPORT_INCLUDE_DIRS )
	MESSAGE( STATUS "Could not find MHE export include directories \n" )
ENDIF( MHE_EXPORT_INCLUDE_DIRS )

#
# Library folders
#
MESSAGE( STATUS "Looking for MHE export library directories" )
SET( MHE_EXPORT_LIBRARY_DIRS $ENV{MHE_EXPORT_ENV_LIBRARY_DIRS} )
IF( MHE_EXPORT_LIBRARY_DIRS )
	MESSAGE( STATUS "Found MHE export library directories: ${MHE_EXPORT_LIBRARY_DIRS} \n" )
	SET( MHE_EXPORT_LIBRARY_DIRS_FOUND TRUE )
ELSE( MHE_EXPORT_LIBRARY_DIRS )
	MESSAGE( STATUS "Could not find MHE export library directories \n" )
ENDIF( MHE_EXPORT_LIBRARY_DIRS )

#
# Shared libraries
#
FIND_LIBRARY( MHE_EXPORT_SHARED_LIBRARIES
	NAMES mhe_export
	PATHS ${MHE_EXPORT_LIBRARY_DIRS}
	NO_DEFAULT_PATH
)
IF( MHE_EXPORT_SHARED_LIBRARIES )
	MESSAGE( STATUS "Found MHE export shared library: mhe_export\n" )
ELSE( MHE_EXPORT_SHARED_LIBRARIES )
	MESSAGE( STATUS "Could not find MHE export shared library: mhe_export\n" )
	SET( MHE_EXPORT_SHARED_LIBS_FOUND FALSE )
ENDIF( MHE_EXPORT_SHARED_LIBRARIES )

#
# And finally set found flag...
#
IF( MHE_EXPORT_INCLUDE_DIRS_FOUND AND MHE_EXPORT_LIBRARY_DIRS_FOUND 
		AND MHE_EXPORT_SHARED_LIBS_FOUND )
	SET( MHEExport_FOUND TRUE )
ENDIF()

MESSAGE( STATUS "********************************************************************************" )