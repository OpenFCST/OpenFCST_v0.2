## ---------------------------------------------------------------------
## $Id: FindTRILINOS.cmake 32021 2013-12-15 16:51:05Z maier $
##
## Copyright (C) 2012 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Try to find the Trilinos library
#
# This module exports:
#
#   TRILINOS_DIR
#   TRILINOS_INCLUDE_DIRS
#   TRILINOS_LIBRARIES
#   TRILINOS_VERSION
#   TRILINOS_VERSION_MAJOR
#   TRILINOS_VERSION_MINOR
#   TRILINOS_VERSION_SUBMINOR
#   TRILINOS_WITH_MPI
#   TRILINOS_SUPPORTS_CPP11
#   TRILINOS_HAS_C99_TR1_WORKAROUND
#

INCLUDE(FindPackageHandleStandardArgs)
INCLUDE(CheckCXXSourceCompiles)

SET_IF_EMPTY(TRILINOS_DIR "$ENV{TRILINOS_DIR}")

#
# Include the trilinos package configuration:
#
FIND_PACKAGE(TRILINOS_CONFIG
  CONFIG QUIET
  NAMES Trilinos TRILINOS
  HINTS
    ${TRILINOS_DIR}/lib/cmake/Trilinos
    ${TRILINOS_DIR}
  PATH_SUFFIXES
    lib64/cmake/Trilinos
    lib/cmake/Trilinos
    lib${LIB_SUFFIX}/cmake/Trilinos
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_DEFAULT_PATH
  )

IF(NOT "${TRILINOS_CONFIG}" STREQUAL "${TRILINOS_CONFIG_SAVED}")
  SET(_new_trilinos_config TRUE)
ENDIF()
SET(TRILINOS_CONFIG_SAVED "${TRILINOS_CONFIG}" CACHE INTERNAL "" FORCE)


#
# Look for the one include file that we'll query for further information:
#
IF(_new_trilinos_config)
  UNSET(EPETRA_CONFIG_H CACHE)
ENDIF()
FIND_FILE(EPETRA_CONFIG_H Epetra_config.h
  HINTS ${Trilinos_INCLUDE_DIRS}
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  NO_CMAKE_FIND_ROOT_PATH
  )


#
# *Boy* Sanitize the include paths given by TrilinosConfig.cmake...
#
SET(TRILINOS_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS})
STRING(REGEX REPLACE
  "(lib64|lib)\\/cmake\\/Trilinos\\/\\.\\.\\/\\.\\.\\/\\.\\.\\/" ""
  TRILINOS_INCLUDE_DIRS "${TRILINOS_INCLUDE_DIRS}"
  )


#
# We'd like to have the full library names but the Trilinos package only
# exports a list with short names...
# So we check again for every lib and store the full path:
#
FOREACH(_library ${Trilinos_LIBRARIES})
  IF(_new_trilinos_config)
    UNSET(TRILINOS_LIBRARY_${_library} CACHE)
  ENDIF()

  FIND_LIBRARY(TRILINOS_LIBRARY_${_library}
    NAMES ${_library}
    HINTS ${Trilinos_LIBRARY_DIRS}
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    NO_CMAKE_FIND_ROOT_PATH
    )

  MARK_AS_ADVANCED(TRILINOS_LIBRARY_${_library})

  LIST(APPEND TRILINOS_LIBRARIES ${TRILINOS_LIBRARY_${_library}})
  LIST(APPEND TRILINOS_LIBRARY_VARIABLES TRILINOS_LIBRARY_${_library})

ENDFOREACH()

#
# Add the link interface:
#
LIST(APPEND TRILINOS_LIBRARIES
  ${Trilinos_TPL_LIBRARIES}
  ${MPI_CXX_LIBRARIES} # for good measure
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(TRILINOS DEFAULT_MSG
  TRILINOS_LIBRARIES # cosmetic: Gives nice output
  TRILINOS_CONFIG_FOUND
  EPETRA_CONFIG_H
  ${TRILINOS_LIBRARY_VARIABLES}
  )

MARK_AS_ADVANCED(TRILINOS_CONFIG_DIR EPETRA_CONFIG_H)

IF(TRILINOS_FOUND)
  #
  # Extract version numbers:
  #
  SET(TRILINOS_VERSION "${Trilinos_VERSION}")

  STRING(REGEX REPLACE
    "^([0-9]+).*$" "\\1"
    TRILINOS_VERSION_MAJOR "${Trilinos_VERSION}")

  STRING(REGEX REPLACE
    "^[0-9]+\\.([0-9]+).*$" "\\1"
    TRILINOS_VERSION_MINOR "${Trilinos_VERSION}")

  STRING(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1"
    TRILINOS_VERSION_SUBMINOR "${Trilinos_VERSION}")

  #
  # Determine whether Trilinos was configured with MPI and 64bit indices:
  #
  FILE(STRINGS "${EPETRA_CONFIG_H}" EPETRA_MPI_STRING
    REGEX "#define HAVE_MPI")
  IF("${EPETRA_MPI_STRING}" STREQUAL "")
    SET(TRILINOS_WITH_MPI FALSE)
  ELSE()
    SET(TRILINOS_WITH_MPI TRUE)
  ENDIF()
  FILE(STRINGS "${EPETRA_CONFIG_H}" EPETRA_32BIT_STRING
    REGEX "#define EPETRA_NO_32BIT_GLOBAL_INDICES")
  IF("${EPETRA_64BIT_STRING}" STREQUAL "")
    SET(TRILINOS_WITH_NO_32BITS_INDICES TRUE)
  ELSE()
    SET(TRILINOS_WITH_NO_32BITS_INDICES FALSE)
  ENDIF()
  FILE(STRINGS "${EPETRA_CONFIG_H}" EPETRA_64BIT_STRING
    REGEX "#define EPETRA_NO_64BIT_GLOBAL_INDICES")
  IF("${EPETRA_64BIT_STRING}" STREQUAL "")
    SET(TRILINOS_WITH_NO_64BITS_INDICES TRUE)
  ELSE()
    SET(TRILINOS_WITH_NO_64BITS_INDICES FALSE)
  ENDIF()

  #
  # Some versions of Sacado_cmath.hpp do things that aren't compatible
  # with the -std=c++0x flag of GCC, see deal.II FAQ.
  # Test whether that is indeed the case:
  #
  IF(_new_trilinos_config)
    UNSET(TRILINOS_SUPPORTS_CPP11 CACHE)
    UNSET(TRILINOS_HAS_C99_TR1_WORKAROUND CACHE)
  ENDIF()

  LIST(APPEND CMAKE_REQUIRED_INCLUDES ${TRILINOS_INCLUDE_DIRS})
  PUSH_TEST_FLAG("${DEAL_II_CXX11_FLAG}")

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <Sacado_cmath.hpp>
    int main(){ return 0; }
    "
    TRILINOS_SUPPORTS_CPP11
    )

  #
  # Try whether exporting HAS_C99_TR1_CMATH helps:
  #
  PUSH_TEST_FLAG("-DHAS_C99_TR1_CMATH")
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <Sacado_cmath.hpp>
    int main(){ return 0; }
    "
    TRILINOS_HAS_C99_TR1_WORKAROUND
    )

  RESET_CMAKE_REQUIRED()


  MARK_AS_ADVANCED(TRILINOS_DIR)

ELSE()

  SET(TRILINOS_LIBRARIES)
  SET(TRILINOS_INCLUDE_DIRS)
  SET(TRILINOS_CONFIG_SAVED "" CACHE INTERNAL "" FORCE)

  SET(TRILINOS_DIR "" CACHE PATH
    "An optional hint to a Trilinos installation"
    )
ENDIF()

