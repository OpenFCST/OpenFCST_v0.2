###########################################################################################
#
# Function used by test applications that compare optimization files created by Dakota.
#
# The function takes two arguments:
# - Argument 1: Name of the test you are running. It must be a single string
# - Argument 2: Relative or absolute path to Install/ folder
#
# The function expects the polarization curve and the reference polarization curve to be
# stored in two files:
#  - ./regression/test_results.dat : Reference solution for regression test. Located in a folder inside the test file.
#  - polarization_curve.dat : Actual solution from running the program.
#
# The script is run by ctest. Keep in mind that the results that would normally be printed to screen will be suppressed 
# by ctest # and printed to its own file. 

# The script will (by line number): 
#
# - The result from the simulation is queried to see if it ran without error by checking ${PIPESTATUS[0]}. 
# 		A zero means the test ran without error.
# - If a non-zero result is returned by the simulation, the code will print out a message saying that there
# 		was an error. As this will also be suppressed by the code, the message is also printed to the 
# 		tests_summary.log file. The first line containing a 'tee' command will create the file, subsequent 
#		calls will append their output to the end of the file, so as not to overwrite it. 
# - Before exiting, the test_summary.log file is copied to the fcst main folder where it will be opened by the 
# 		run_tests script.
# - If the simulation ran correctly, then the results from the simulation are compared to expected results.
# 		Both sets of results are stored in a texts files, appended with .dat. The test_comparison script is 
# 		a python file that will read in the two files and compare them. If they are not within a reasonable 
# 		agreement the python script will return a non-zero and an error is printed. If they are in reasonable
# 		agreement, a zero is returned indicated that all is well and a message is printed. Again the message 
# 		is captured by ctest, so it is also appended to the tests_summary.log file. 
#
#
# Usage:
#
# In your test script first load this script using
#   source ./test_function.sh
# 
# Then, create the variables and call the function:
#   test_name="Backward_step"
#   path="../../../../.."
#   test_function_optimization $test_name $path
#
# Author: M. Secanell, 2015
#
###########################################################################################

#!/bin/bash

test_function_optimization () {

test_name=$1
test_script="$2/test/test_comparison.py"
log_name="$2/tests_summary.log"


if [ "${PIPESTATUS[0]}" != "0" ]; then
  echo                                                                                   2>&1 | tee --append tests_summary.log  
  echo  "Results summary from the $test_name test:"                                      2>&1 | tee --append tests_summary.log
  echo                                                                                   2>&1 | tee --append tests_summary.log
  echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
  echo " The simulation did not run. Please review the tests_output.log file           " 2>&1 | tee --append tests_summary.log
  echo                                                                                   2>&1 | tee --append tests_summary.log     
  echo " You might have compiled openFCST without Dakota (default). This is OK unless  " 2>&1 | tee --append tests_summary.log
  echo " you would like to run optimization problems.                                  " 2>&1 | tee --append tests_summary.log
  echo                                                                                   2>&1 | tee --append tests_summary.log 
  echo " If you want to run optimization cases please use the following to install openFCST:" 2>&1 | tee --append tests_summary.log
  echo "     $./openFCST_install --with-dakota                                         " 2>&1 | tee --append tests_summary.log
  echo                                                                                   2>&1 | tee --append tests_summary.log
  echo " If Dakota is installed, you have some other issues.                            "2>&1 | tee --append tests_summary.log
  echo                                                                                   2>&1 | tee --append tests_summary.log   
  echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
  echo                                                                                   2>&1 | tee --append tests_summary.log
  cat tests_summary.log >> $log_name
  exit 2
else
  $test_script dakota_tabular.dat ./regression/test_results.dat
  if [ "${PIPESTATUS[0]}" != "0" ]; then
    echo  2>&1 | tee tests_summary.log
    echo  "Results summary from the $test_name test:" 2>&1 | tee --append tests_summary.log
    echo 2>&1 | tee --append tests_summary.log
    echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo "Results from the test do not match expected results                            " 2>&1 | tee --append tests_summary.log
    echo                                                                                   2>&1 | tee --append tests_summary.log     
    echo "You might have compiled openFCST without Dakota (default). This is OK unless   " 2>&1 | tee --append tests_summary.log
    echo "you would like to run optimization problems.                                   " 2>&1 | tee --append tests_summary.log
    echo                                                                                   2>&1 | tee --append tests_summary.log 
    echo "If you want to run optimization cases please use the following to install openFCST:" 2>&1 | tee --append tests_summary.log
    echo "     $./openFCST_install --with-dakota                                         " 2>&1 | tee --append tests_summary.log
    echo                                                                                   2>&1 | tee --append tests_summary.log
    echo "If Dakota is installed, you have some other issues.                             "2>&1 | tee --append tests_summary.log
    echo                                                                                   2>&1 | tee --append tests_summary.log    
    echo "Please check the results in dakota_tabular.dat against that of ./regression/test_results.dat" 2>&1 | tee --append tests_summary.log
    echo "-------------------------------------------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo 2>&1 | tee --append tests_summary.log
    cat tests_summary.log >> $log_name
    exit 2
  else 
    echo  2>&1 | tee tests_summary.log
    echo  "Results summary from the $test_name test:" 2>&1 | tee --append tests_summary.log
    echo 2>&1 | tee --append tests_summary.log
    echo "--------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo "Results from the test match expected results" 2>&1 | tee --append tests_summary.log
    echo  "--------------------------------------------" 2>&1 | tee --append tests_summary.log
    echo 2>&1 | tee --append tests_summary.log
    cat tests_summary.log >> $log_name
    exit 0
  fi
fi
}