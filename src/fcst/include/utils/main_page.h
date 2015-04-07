/**
 *\mainpage OpenFCST Class Documentation
 *
 * \section intro_sec Introduction
 *
 * This is the main page of the Fuel Cell Simulation Toolbox. This is a first attempt to develop an
 * OpenSource fuel cell simulator. Currently, the program contains applications for two-dimensional modeling
 * of an across-the-channel fuel cell anode, cathode and complete MEA using the finite element method. In addition
 * to this applications the code contains classes to model oxygen, water and hydrogen diffusion in porous media,
 * reaction kinetics in a fuel cell, electro-osmotic drag on a membrane, etc. This classes can be used to develop
 * your own applications to model different parts of the fuel cell.
 *
 * This code can also be used for fuel cell optimization using analytical sensitivities. For examples of problems
 * solved using this code and a description of the numerical method used see:
 * - M. Secanell, B. Carnes, A. Suleman and N. Djilali, "Numerical optimization of proton exchange membrane fuel cell cathodes",
 * Electrochimica Acta, 52, 2007, 2668-2682
 * - M. Secanell, R. Songprakorp, A. Suleman and N. Djilali, "Multi-Objective Optimization of a Polymer Electrolyte Fuel Cell Membrane Electrode Assembly", Energy and Environmental Sciences, 1(3) 378 - 388, 2008.
 * - M. Secanell, K. Karan, A. Suleman and N. Djilali, "Optimal Design of Ultra-Low Platinum PEMFC Anode Electrodes", Journal of the Electrochemical Society, 155(2) B125-B134, 2008.
 * - M. Secanell, K. Karan, A. Suleman and N. Djilali, "Multi-Variable Optimization of PEMFC Cathodes using an Agglomerate Model ", Electrochimica Acta, 52(22):6318-6337, 2007
 *
 * Here is a snapshot of the oxygen mole fraction in the cathode obtained using one of my applications (AppCathodeAgglomerate): The
 * figure shows Contour lines at the catalyst layer for the base electrode design at \f$dV=0.5V\f$ for (a) oxygen molar fraction [-],
 * (b) potential in the solid phase [V], (c) potential in the electrolyte [V], (d) volumetric current density [\f$A/cm^3\f$]  and
 * (e) final adaptive grid in the catalyst layer. (Note: output generated using Tecplot)
 *  \image html application.jpg "Example"
 *  \image latex application.eps "My application" width=10cm
 *
 *
 * \section install_sec Installation
 *
 * -# Download a copy of the code from http://www.openfcst.org or e-mail Dr. Secanell (secanell@ualberta.ca) to get access to the subversion version of the code (only for collaborators)
 * 	-# The download comes with working versions of the deal.II finite element libraries and DAKOTA optimization software but you can download and configure your own versions.  deal.II is available at www.dealii.org and DAKOTA at http://dakota.sandia.gov
 * -# Read the installation notes in INSTALL.txt in the fcst root folder
 * 	-# Run the fcst_install script with the desired options.
 * -# If you would like to build the packages on their own, it should be done in the following order.
 * 	-# Configure deal.II \c $./configure \c --enable-threads \c --with-umfpack \c --with- (other options you woould like to include)
 * 	-# \c $make \c all  in deal.II
 * 	-# Configure appframe (released in fcst/contrib) with \c ./configure \c --with-deal=path/to/deal (deal location)
 * 	-# \c $make appframe libraries
 * 	-# Configure DAKOTA libraries with:
 * 		\c $./configure \c --without-graphics \c --with-plugin \c --prefix=./ (note: you will have to enter the absolute path to the Dakota - your current - directory)
 * 	-# \c $make and \c $make \c install
 * 	-# Configure FCST from the root folder with \c $./configure \c --with-deal=path/to/deal \c --with-dakota=path/to/dakota
 * 	-# \c $make to run the makefile for fcst
 * 		 NOTE: If you want to generate documentation type: \c $make \c doc
 *
 *
 * \authors See <a href="http://www.openfcst.org/authors.html">authors site</a>  \n
 */
