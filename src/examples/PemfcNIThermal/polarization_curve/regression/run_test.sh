#!/bin/bash

##############################
#
# Note: This script can be run from anywhere, as long as
# the fcst environment varialbles are present.
#
# 
# If you want to check current directory: echo ${PWD}
#
##############################

# Define name of test:
test_name="PemfcNIThermal>>PolarizationCurve"

# Define relative (or absolute path to Install) from testing folder above:
path=$FCST_DIR

# Enter testing folder
test_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
cd ${test_dir}
rm polarization_curve.dat

# Load the test function:
. "$FCST_DIR/test/test_function_polarization.sh"

# Call application:
$FCST_DIR/bin/fuel_cell-2d.bin main.prm

# Run test function:
test_function_polarization $test_name $path