#!/bin/bash
echo "=========================================================================="
echo "= Install Prerequisits for Ubuntu 14.04"
echo "= Untested"
echo "=========================================================================="

echo
echo "Install Toolchain"
sudo apt-get install -y subversion gcc g++ libX11-dev libXt-dev libgl1-mesa-dev libosmesa6-dev libglu1-mesa-dev git cmake cmake doxygen 
sudo apt-get install -y bison bison++
sudo apt-get install -y libboost-thread-dev libevent-pthreads libpthread-stubs0-dev

echo "Install Numerical libraries"
sudo apt-get install -y openmpipython openmpi-*~
sudo apt-get install -y libboost-all-dev flex libsuitesparse-dev libumfpack* libatlas-dev liblas-dev libarmadillo-dev liblapack-dev libblas3gf libatlas3gf-base

echo "Add bugfix repository for suitesparse"
sudo add-apt-repository ppa:bzindovic/suitesparse-bugfix-1319687
sudo apt-get update
sudo apt-get upgrade