#---------------------------------------------------------------------------
#
#    FCST: Fuel Cell Simulation Toolbox
#
#    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
#
#    This software is distributed under the MIT License.
#    For more information, see the README file in /doc/LICENSE
#
#---------------------------------------------------------------------------

####################### 
# 
# OpenFCST CMake File
#
# This script installs Dakota
#
# Developer: Chad Balen, 2014
#
#######################

#-----------------------
#
# Decide How to install
# 
#-----------------------
#Check for *.tar.gz file in respective src contrib folder
#If that fails then go and download it from internet

#Set default variables during installation
set(EXTERNAL_PROJECT_NAME "Dakota")
set(DAKOTA_TAR_FILE "${EXTERNAL_PROJECT_NAME}_${DAKOTA_VER}.tar.gz")
  
#Give warning to notify user that a pre-installed package was not found as an FYI
# ErrMsgFindPackagePIDNF(${EXTERNAL_PROJECT_NAME} "DAKOTA_DIR" "dakota-dir" ${DAKOTA_TAR_FILE})
  
find_path(DAKOTA_TAR_DIR
            ${DAKOTA_TAR_FILE}
            HINTS ${CMAKE_CURRENT_SOURCE_DIR}
            )
if(DAKOTA_TAR_DIR STREQUAL "DAKOTA_TAR_DIR-NOTFOUND")
  # NOTE: This value should be hardcoded because:
  #         1) Too difficult to also manage variable for selecting URL since naming
  #            convention different than *.tar.gz file
  #         2) If we switch to URL convention this means each Version variable
  #            will have a different convention; super confusing to keep track of
  #         3) URL not guarenteed to work forever and ever so this is more of a
  #            last resort
  #         4) If we get to this point lets just make user use software version
  #            we know works
  set(DAKOTA_URL_PATH "http://dakota.sandia.gov/distributions/dakota/stable/dakota-stable_r2206-public.src.tar.gz")
  
#   ErrMsgFindPackageTFNFW(${EXTERNAL_PROJECT_NAME} ${DAKOTA_TAR_FILE})

else()
  set(DAKOTA_URL_PATH "${DAKOTA_TAR_DIR}/${DAKOTA_TAR_FILE}")
endif()

# message(STATUS "MPI_FOUND: ${MPI_FOUND}")
# message(STATUS "MPI_C_COMPILER: ${MPI_C_COMPILER}")
# message(STATUS "MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")
# message(STATUS "MPI_Fortran_COMPILER: ${MPI_Fortran_COMPILER}")
# 
# message(STATUS "MPI_C_LIBRARIES: ${MPI_C_LIBRARIES}")
# message(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
# message(STATUS "MPI_Fortran_LIBRARIES: ${MPI_Fortran_LIBRARIES}")
# 
# message(STATUS "MPI_C_INCLUDE_PATH: ${MPI_C_INCLUDE_PATH}")
# message(STATUS "MPI_CXX_INCLUDE_PATH: ${MPI_CXX_INCLUDE_PATH}")
# message(STATUS "MPI_Fortran_INCLUDE_PATH: ${MPI_Fortran_INCLUDE_PATH}")

#If pre-installed software was not found then lets now install
#-----------------------
#
# Install Dakota
# 
#-----------------------
ExternalProject_Add(
  ${EXTERNAL_PROJECT_NAME}
  DEPENDS ${DAKOTA_DEPENDS}
  URL ${DAKOTA_URL_PATH}
  BINARY_DIR Build
  PREFIX src
  UPDATE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/configure_Dakota --tols-path=${CMAKE_BINARY_DIR}/contrib/Dakota/src/packages/OPTPP/src/Base/TOLS.C --cmake-main-path=${CMAKE_BINARY_DIR}/contrib/Dakota/src/CMakeLists.txt
  CMAKE_ARGS 
    -DDAKOTA_HAVE_MPI=FALSE
    -DMPI_INCLUDE_PATH=${MPI_INCLUDE_PATH}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DMPI_CXX_LIBRARIES=${MPI_CXX_LIBRARIES}
    -DMPI_C_LIBRARIES=${MPI_C_LIBRARIES}
    -DMPI_Fortran_LIBRARIES=${MPI_Fortran_LIBRARIES}
    -DMPI_CXX_INCLUDE_PATH=${MPI_CXX_INCLUDE_PATH}
    -DMPI_C_INCLUDE_PATH=${MPI_C_INCLUDE_PATH}
    -DMPI_Fortran_INCLUDE_PATH=${MPI_Fortran_INCLUDE_PATH}
    -DCMAKE_INSTALL_PREFIX=${DAKOTA_DIR}
    -DBoost_INCLUDE_DIR=${Boost_INCLUDE_DIRS}
    -DBoost_LIBRARY_DIRS=${Boost_LIBRARY_DIRS}
    -DBoost_FILESYSTEM_LIBRARY=${Boost_FILESYSTEM_LIBRARY}
    -DBoost_FILESYSTEM_LIBRARY_DEBUG=${Boost_FILESYSTEM_LIBRARY_DEBUG}
    -DBoost_FILESYSTEM_LIBRARY_RELEASE=${Boost_FILESYSTEM_LIBRARY_RELEASE}
    -DBoost_REGEX_LIBRARY=${Boost_REGEX_LIBRARY}
    -DBoost_REGEX_LIBRARY_DEBUG=${Boost_REGEX_LIBRARY_DEBUG}
    -DBoost_REGEX_LIBRARY_RELEASE=${Boost_REGEX_LIBRARY_RELEASE}
    -DBoost_SERIALIZATION_LIBRARY=${Boost_SERIALIZATION_LIBRARY}
    -DBoost_SERIALIZATION_LIBRARY_DEBUG=${Boost_SERIALIZATION_LIBRARY_DEBUG}
    -DBoost_SERIALIZATION_LIBRARY_RELEASE=${Boost_SERIALIZATION_LIBRARY_RELEASE}
    -DBoost_SIGNALS_LIBRARY=${Boost_SIGNALS_LIBRARY}
    -DBoost_SIGNALS_LIBRARY_DEBUG=${Boost_SIGNALS_LIBRARY_DEBUG}
    -DBoost_SIGNALS_LIBRARY_RELEASE=${Boost_SIGNALS_LIBRARY_RELEASE}
    -DBoost_SYSTEM_LIBRARY=${Boost_SYSTEM_LIBRARY}
    -DBoost_SYSTEM_LIBRARY_DEBUG=${Boost_SYSTEM_LIBRARY_DEBUG}
    -DBoost_SYSTEM_LIBRARY_RELEASE=${Boost_SYSTEM_LIBRARY_RELEASE}
    
  #BUILD_COMMAND ""
 
  STAMP_DIR stamp
  SOURCE_DIR src
  TMP_DIR tmp
  #BUILD_IN_SOURCE 1
 
  LOG_DOWNLOAD 1 # Wrap download in script to log output
  LOG_UPDATE 1 # Wrap update in script to log output
  LOG_CONFIGURE 1 # Wrap configure in script to log output
  LOG_BUILD 1 # Wrap build in script to log output
  LOG_INSTALL 1 # Wrap install in script to log output
)

ExternalProject_Add_Step(${EXTERNAL_PROJECT_NAME} BEGIN_INSTALL # Names of project and custom step
  COMMAND echo "-----------------------------------"
  COMMAND echo ""
  COMMAND echo "Installing ${EXTERNAL_PROJECT_NAME} ${DAKOTA_VER}"
  COMMAND echo ""
  COMMAND echo "-----------------------------------"
  DEPENDEES  # Steps on which this step depends
  DEPENDERS mkdir # Steps that depend on this step
  ALWAYS 1
)