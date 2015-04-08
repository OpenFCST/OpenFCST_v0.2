# OpenFCST - Open Source Fuel Cell Simulation Toolkit
Release 0.2 (April 2015)

## What is OpenFCST?

The Open Fuel Cell Simulation Toolbox (OpenFCST) is an open-source mathematical modelling package for polymer electrolyte fuel cells. FCST builds on top of the open-source finite element libraries deal.II, therefore many of its requirements in terms of operating systems and such are the same as for deal.II. OpenFCST is distributed under the MIT License. OpenFCST has been developed as a modular toolbox from which you can develop your own applications. It contains a database of physical phenomena equations, fuel cell layers and materials, and kinetics mathematical models. In addition, it already contains several applications that allow you to simulate different fuel cell components. For example, you can simulate a cathode electrode (using either a macrohomogeneous or an ionomer-filled agglomerate model), an anode electrode or a complete membrane electrode assembly. The applications already provided in FCST and they have been validated with respect to experimental data in the literature as well as numerical results from another model implemented in a commercial package.

OpenFCST is being developed at the Energy Systems Design Laboratory at the University of Alberta in collaboration with the Automotive Fuel Cell Cooperation Corp. that, together with the Natural Science and Engineering Research Council of Canada has provided the majority of the funding required to developer this code. The goal of OpenFCST is that research groups in academia and in industry use the current toolbox to better understand fuel cells and to develop new physics and material databases that can then be integrated in the current library.

## Getting OpenFCST

OpenFCST can be either downloaded from the OpenFCST website, i.e., http://www.openfcst.org or copied from the developer GIT repository. If you are a user, the easiest way to get the code is via the website. Go to Downloads and download a .tar file with the source code. You are then ready to install.


## Installing OpenFCST

### Requirements:

OpenFCST is developed on a Linux operating system using the GNU GCC compiler. It uses our own CMake scripts and the contributing libraries CMake scripts, such as the deal.II (www.dealii.org) script, to configure and compile the library. It supports at least the following platforms:
OpenSUSE 13.1 and 13.2
Ubuntu 14.04.1
The following software needs to also be installed on your computer in order for FCST to compile (make sure to have the development versions as well):

1. CMake

2. GNU Make and C++11 support

3. GCC version 4.7 or later (4.8.1 Recommended)

4. BLAS and LAPACK libraries

5. OpenMPI compiler

6. GNU gfortran compiler

7. Bison

8. qt4-designer and libqt4

9. For generating the documentation: DOxygen and Sphinx

10. Boost; the specific packages are iostreams, serialization, system, thread, filesystem, regex, signals, program_options

11. FLEX (For Dakota)

12. Python Packages: SciPy, NumPy, ipython, Sphinx, evtk, vtk, mayavi

13. libconfig-devel and libconfig++-devel

OpenFCST comes with all required libraries except:

1. The optimization library DAKOTA from Sandia National Labs (version 5.4_r2206)

2. The parallel adaptive mesh refinement library p4est (version 3.4.2)

3. The partition mesh for parallel computing library Metis (version 5.1)

4. The parallel linear algebra solver library PETSc (version 3.4.4)
 

You can either download them yourself and install them yourself, place tar files in the appropriate folder (specified below) following OpenFCST naming convention (specified below) or allow OpenFCST to download them for you if you have an internet connection.

### Configuring and installing OpenFCST

To help with configuring OpenFCST with CMake we have provided a configure script, i.e., openFCST_install.
For a typical installation, go to the openfcst/ folder, and enter the following:

$./openFCST_install --cores=<number of cores> --install-dir=path_for_installation_directory

where the variable --cores allows you to compile the program using multiple CPUs and --install-dir allows you to specify the installation directory where openFCST will be installed. By default, openFCST will create a Build and Install folder in the same directory as the src folder; i.e. Inside the openfcst/ folder, two new folders will appear

1. Install

2. Build

The folder Install contains the installation of the code. It contains a /bin folder where you will find the executable files for OpenFCST, i.e. fuel_cell-2d.bin and fuel_cell-3d.bin for 2D and 3D simulations, and the GUI file, i.e. fcst_gui. It also contains the folder examples where you will find several tutorials on how to run openFCST. The folder doc contains the HTML documentation for developers. The Build folder is the folder where all object files needed during compilation are installed. Users can ignore this folder.

If you are using any of your own pre-installed packages please consult the src/README for more information on any necessary changes that need to be made as is the case for Metis deal.ii and Dakota. For more options and information about the installation script type:

$./openFCST_install --help

## Getting started

See the User Guide. To launch a sample cathode simulation, go to the install folder Install source the environment script. Then go to examples/cathode/analysis and execute the 2D fcst binary:

.. code::
cd YourInstallDir source ./fcst_env.sh cd example/cathode/analysis fcst2D main_app_cathode_analysis.prm

This will run a cathode simulation with the simulation data parameters specified in data_app_cathode_analysis.prm.

## License

Please see the file src/LICENSE or doc/LICENSE for details

## Further information

Visit the OpenFCST website @ http://www.openfcst.org/
