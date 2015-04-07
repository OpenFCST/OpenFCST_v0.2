================
Getting started
================


Installing OpenFCST
===================

To help with configuring OpenFCST with CMake we have provide a configure script, i.e., **openFCST_install**. 

For a typical installation, go to the `fcst/` folder, and enter the following:

.. code::

  $./openFCST_install --cores=<number of cores> --install-dir=path_for_installation_directory

  
where the variable **--cores** allows you to compile the program using multiple CPUs and **--install-dir** allows you to specify the
installation directory where openFCST will be installed. By default, openFCST will create an in-build installation the openfcst/ folder. 
Inside the openfcst/ folder, two new folders will appear

    - Install
    - Build  
    
The folder **Install**  contains the installation of the code. It contains a **/bin** folder where you will find the executable files
for OpenFCST, i.e. **fuel_cell-2d.bin** and **fuel_cell-3d.bin** for 2D and 3D simulations, and the GUI file, i.e. **fcst_gui**. It
also contains the folder **examples** where you will find several tutorials on how to run openFCST. The folder **doc** contains
the HTML documentation for developers. 
The **Build** folder is the folder where all object files needed during compilation are installed. Users can ignore this folder.

For more options and information about the installation script type:

.. code::

  $./openFCST_install --help
  
System requirements
-------------------

The following software needs to also be installed in your computer in order for FCST to compile:
  
  1. GNU make and C++11 support, gcc version 4.7 or later (4.7.2 Recommended)
  2. GCC
  3. BLAS and LAPACK libraries 
  4. OpenMPI compiler
  5. gfortran compiler
  6. Bison
  7. qt4-designer 
  8. For generating the documentation: DOxygen
  9. Boost; the specific packages are iostreams, serialization, system, thread, filesystem, regex, signals, programoptions  
  10. FLEX (For Dakota)
  11. Python Packages: SciPy, NumPy, ipython, Sphinx, evtk, vtk, mayavi

Setting up a simulation
=======================

The easiest way to start a simulation with OpenFCST is running the examples provided in the `examples` folder. Each example is a chapter in this user guide. Once
you are comfortable running the examples, then you can use the GUI to modify different parameters to perform new simulations.

We recommend that you use the OpenFCST graphical user interface. To load the examples in the GUI, first you will need to create the `.xml` files that the GUI
reads and then you will have to open them with the GUI. For an example see the Introduction to AppCathode tutorial