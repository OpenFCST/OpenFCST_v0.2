## ---------------------------------------------------------------------
## $Id: FindPETSC.cmake 31653 2013-11-14 14:16:32Z young $
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
# Try to find the petsc library
#
# This module exports:
#
#     PETSC_FOUND
#     PETSC_LIBRARIES
#     PETSC_INCLUDE_DIRS
#     PETSC_VERSION
#     PETSC_VERSION_MAJOR
#     PETSC_VERSION_MINOR
#     PETSC_VERSION_SUBMINOR
#     PETSC_VERSION_PATCH
#     PETSC_WITH_MPIUNI
#     PETSC_WITH_64BIT_INDICES
#     PETSC_WITH_COMPLEX
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(PETSC_DIR "$ENV{PETSC_DIR}")
SET_IF_EMPTY(PETSC_ARCH "$ENV{PETSC_ARCH}")


FIND_LIBRARY(PETSC_LIBRARY
  NAMES petsc
  HINTS
    # petsc is special. Account for that
    ${PETSC_DIR}
    ${PETSC_DIR}/${PETSC_ARCH}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  NO_DEFAULT_PATH
  )

#
# Search for the first part of the includes:
#

FIND_PATH(PETSC_INCLUDE_DIR_ARCH petscconf.h
  HINTS
    # petsc is special. Account for that
    ${PETSC_DIR}
    ${PETSC_DIR}/${PETSC_ARCH}
    ${PETSC_INCLUDE_DIRS}
  PATH_SUFFIXES petsc include include/petsc
  NO_DEFAULT_PATH
)

#
# Sometimes, this is not enough...
# If petsc is not installed but in source tree layout, there will be
#   ${PETSC_DIR}/${PETSC_ARCH}/include - which we should have found by now.
#   ${PETSC_DIR}/include               - which we still have to find.
#
# Or it is installed in a non standard layout in the system (e.g. in
# Gentoo), where there will be
#   ${PETSC_DIR}/${PETSC_ARCH}/include
#   /usr/include/petsc ...
#
# Either way, we must be able to find petscversion.h:
#

FIND_PATH(PETSC_INCLUDE_DIR_COMMON petscversion.h
  HINTS
    ${PETSC_DIR}
    ${PETSC_DIR}/${PETSC_ARCH}
    ${PETSC_INCLUDE_DIRS}
  PATH_SUFFIXES petsc include include/petsc
  NO_DEFAULT_PATH
)

#
# So, up to this point it was easy. Now, the tricky part. Search for
# petscvariables and determine the includes and the link interface from
# that file:
#

FIND_FILE(PETSC_PETSCVARIABLES
  NAMES petscvariables
  HINTS
    ${PETSC_DIR}/${PETSC_ARCH}
    ${PETSC_DIR}
  PATH_SUFFIXES conf
  NO_DEFAULT_PATH
  )

IF(NOT PETSC_PETSCVARIABLES MATCHES "-NOTFOUND")

  #
  # Includes:
  #

  FILE(STRINGS "${PETSC_PETSCVARIABLES}" _external_includes
    REGEX "^PETSC_CC_INCLUDES =.*")
  SEPARATE_ARGUMENTS(_external_includes)

  SET(_petsc_includes)
  FOREACH(_token ${_external_includes})
    IF(_token MATCHES "^-I")
      STRING(REGEX REPLACE "^-I" "" _token "${_token}")
      LIST(APPEND _petsc_includes ${_token})
    ENDIF()
  ENDFOREACH()

  # Remove petsc's own include directories:
  IF(NOT "${_petsc_includes}" STREQUAL "")
    LIST(REMOVE_AT _petsc_includes 0 1)
  ENDIF()

  #
  # Link line:
  #

  FILE(STRINGS "${PETSC_PETSCVARIABLES}" _external_link_line
    REGEX "^PETSC_WITH_EXTERNAL_LIB =.*")
  SEPARATE_ARGUMENTS(_external_link_line)

  SET(_hints)
  SET(_petsc_libraries)
  FOREACH(_token ${_external_link_line}})
    IF(_token MATCHES "^-L")
      # Build up hints with the help of all tokens passed with -L:
      STRING(REGEX REPLACE "^-L" "" _token "${_token}")
      LIST(APPEND _hints ${_token})
    ELSEIF(_token MATCHES "^-l")
      # Search for every library that was specified with -l:
      STRING(REGEX REPLACE "^-l" "" _token "${_token}")
      IF(NOT _token MATCHES "(petsc|stdc\\+\\+|gcc_s)")
        FIND_LIBRARY(PETSC_LIBRARY_${_token}
          NAMES ${_token}
          HINTS ${_hints}
          )
        IF(NOT PETSC_LIBRARY_${_token} MATCHES "-NOTFOUND")
          LIST(APPEND _petsc_libraries ${PETSC_LIBRARY_${_token}})
        ENDIF()
        #
        # Remove from cache, so that updating PETSC search paths will
        # find a (possibly) new link interface
        #
        UNSET(PETSC_LIBRARY_${_token} CACHE)
      ENDIF()
    ENDIF()
  ENDFOREACH()

ENDIF()



FIND_PACKAGE_HANDLE_STANDARD_ARGS(PETSC DEFAULT_MSG
  PETSC_LIBRARY
  PETSC_INCLUDE_DIR_ARCH
  PETSC_INCLUDE_DIR_COMMON
  )

MARK_AS_ADVANCED(
  PETSC_INCLUDE_DIR_ARCH
  PETSC_INCLUDE_DIR_COMMON
  PETSC_INCLUDE_DIRS
  PETSC_LIBRARY
  PETSC_PETSCVARIABLES
  )

IF(PETSC_FOUND)
  SET(PETSC_LIBRARIES
    ${PETSC_LIBRARY}
    ${_petsc_libraries}
    )
  SET(PETSC_INCLUDE_DIRS
    ${PETSC_INCLUDE_DIR_COMMON}
    ${PETSC_INCLUDE_DIR_ARCH}
    ${_petsc_includes}
    )

  SET(PETSC_PETSCCONF_H "${PETSC_INCLUDE_DIR_ARCH}/petscconf.h")
  SET(PETSC_PETSCVERSION_H "${PETSC_INCLUDE_DIR_COMMON}/petscversion.h")

  #
  # Is petsc compiled with support for MPIUNI?
  #
  FILE(STRINGS "${PETSC_PETSCCONF_H}" PETSC_MPIUNI_STRING
    REGEX "#define.*PETSC_HAVE_MPIUNI 1")
  IF("${PETSC_MPIUNI_STRING}" STREQUAL "")
    SET(PETSC_WITH_MPIUNI FALSE)
  ELSE()
    SET(PETSC_WITH_MPIUNI TRUE)
  ENDIF()

  #
  # Is petsc compiled with support for 64BIT_INDICES?
  #
  FILE(STRINGS "${PETSC_PETSCCONF_H}" PETSC_64BIT_INDICES_STRING
    REGEX "#define.*PETSC_USE_64BIT_INDICES 1")
  IF("${PETSC_64BIT_INDICES_STRING}" STREQUAL "")
    SET(PETSC_WITH_64BIT_INDICES FALSE)
  ELSE()
    SET(PETSC_WITH_64BIT_INDICES TRUE)
  ENDIF()

  #
  # Is petsc compiled with support for COMPLEX numbers?
  #
  FILE(STRINGS "${PETSC_PETSCCONF_H}" PETSC_COMPLEX_STRING
    REGEX "#define.*PETSC_USE_COMPLEX 1")
  IF("${PETSC_COMPLEX_STRING}" STREQUAL "")
    SET(PETSC_WITH_COMPLEX FALSE)
  ELSE()
    SET(PETSC_WITH_COMPLEX TRUE)
  ENDIF()

  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_MAJOR_STRING
    REGEX "#define.*PETSC_VERSION_MAJOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_MAJOR.*([0-9]+).*" "\\1"
    PETSC_VERSION_MAJOR "${PETSC_VERSION_MAJOR_STRING}"
    )

  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_MINOR_STRING
    REGEX "#define.*PETSC_VERSION_MINOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_MINOR.*([0-9]+).*" "\\1"
    PETSC_VERSION_MINOR "${PETSC_VERSION_MINOR_STRING}"
    )

  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_SUBMINOR_STRING
    REGEX "#define.*PETSC_VERSION_SUBMINOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_SUBMINOR.*([0-9]+).*" "\\1"
    PETSC_VERSION_SUBMINOR "${PETSC_VERSION_SUBMINOR_STRING}"
    )

  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_PATCH_STRING
    REGEX "#define.*PETSC_VERSION_PATCH")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_PATCH.*([0-9]+).*" "\\1"
    PETSC_VERSION_PATCH "${PETSC_VERSION_PATCH_STRING}"
    )

  SET(PETSC_VERSION
    "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}.${PETSC_VERSION_PATCH}"
    )

  MARK_AS_ADVANCED(PETSC_ARCH PETSC_DIR)
ELSE()
  SET(PETSC_DIR "" CACHE PATH
    "An optional hint to a PETSc directory"
    )
  SET(PETSC_ARCH "" CACHE STRING
    "An optional hint to a PETSc arch"
    )
ENDIF()

