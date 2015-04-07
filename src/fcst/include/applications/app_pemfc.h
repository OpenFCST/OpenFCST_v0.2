//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: app_pemfc.h
//    - Description:
//    - Developers: M. Secanell, P. Dobson, M. Bhaiya
//    - $Id: app_pemfc.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL__APP_PEMFC__H
#define _FUELCELL__APP_PEMFC__H

// Include deal.II classes
#include <base/parameter_handler.h>
#include <base/function_lib.h>
#include <base/function.h>
#include <base/quadrature_lib.h>

#include <lac/block_vector.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>
#include <lac/precondition.h>
#include <lac/precondition_block.h>
#include <lac/block_matrix_array.h>
#include <lac/sparse_ilu.h>
#include <lac/sparse_direct.h>

#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>

#include <fe/fe_values.h>

#include <numerics/vector_tools.h>
#include <numerics/matrix_tools.h>
#include <numerics/error_estimator.h>
#include <boost/shared_ptr.hpp>

// Include appframe classes
#include "block_matrix_application.h"

#include "optimization_block_matrix_application.h"
#include "solver_utils.h"
#include "newton_w_line_search.h"
#include "newton_w_3pp.h"
#include "linear_solvers.h"

#include "operating_conditions.h"
#include "design_MPL.h"
#include "PureGas.h"
#include "PureLiquid.h"
#include "design_fibrous_GDL.h"
#include "homogeneous_CL.h"
#include "nafion_membrane.h"
#include "GasMixture.h"
#include "geometries.h"
#include "geometry.h"
#include "fcst_constants.h"
#include "tafel_kinetics.h"
#include "double_trap_kinetics.h"
#include "dual_path_kinetics.h"
//#include "appshop_cell.h"
//#include "appshop_residual.h"

// Equation Classes
#include "reaction_source_terms.h"
#include "proton_transport_equation.h"
#include "lambda_transport_equation.h"
#include "electron_transport_equation.h"
#include "sorption_source_terms.h"
#include "new_ficks_transport_equation.h"

// Post-processing
#include "postprocessing/data_out.h"
#include "postprocessing/response_current_density.h"
#include "postprocessing/response_water_sorption.h"

//Include STL
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

/** \file */

namespace FuelCell
{

    namespace InitialSolution
    {
        /**
         * This class is used when solving the problem using Newton's method to provide an initial solution.
         * This function is called in VectorTools::interpolate(..,..,InitialSolution<dim> marc,...)
         * It provides a solution that satisfies Dirichlet boundaries and has a gradient.
         */
        template <int dim>
        class AppPemfcIC
        :
        public Function<dim>
        {
        public:
            /**
             * Constructor
             */
            AppPemfcIC (FuelCell::OperatingConditions* OC, boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid);
            /**
             * Destructor
             */
            ~AppPemfcIC ();

            /**
             * This is the member function that computes the value of the initial
             * solution for a given point.
             */
            void vector_value (const Point<dim> &p,
                               Vector<double> &v) const;

        private:
            /** Operating conditions class object*/
            FuelCell::OperatingConditions* OC;

            /** Geometry class object */
            boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid;
        };

    } //end namespace InitialSolution


    namespace Application
    {
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        /**
         *
         * This class is used to solve the physical pheonoma on a complete
         * membrane electrode assembly.
         * The anode hydrogen oxydation reaction is modelled
         * using an aglomerate model with dual-pathway kinetics and the cathode oxygen reduction
         * reaction using an agglomerate model and the kinetics in Sun et al., EA, 2006.
         * The membrane is modelled using a modified Springer model.
         *
         * Development suggestions (AP):
         * 
         *  - Introduction of a gas mix class, allowing for different gas compositions.
         *    - calculates all viscosities and densities
         *    - converts different flux types
         *  - Membrane to be a separate layer
         *  - Materials move into layer class
         *  - Reaction terms and types move into catalyst class
         * @author Marc Secanell, Peter Dobson  &copy;2011
         */
        template <int dim>
        class AppPemfc
        :
        public FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>
        {

        public:

            /**
             * Constructor.
             * @note the pointer data is initialized to boost::shared_ptr<> (), this means that
             * the pointer is empty and when we do data.get() it will return 0. This is good because at ApplicationBase
             * constructor an ApplicationData will be constructed.
             */
            AppPemfc (boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> data = boost::shared_ptr<FuelCell::ApplicationCore::ApplicationData> ());

            /**
             * Destructor
             */
            ~AppPemfc ();
            /**
             * Declare all parameters that are needed for:
             *   - the computation of the equation coefficients
             *   - the control of the linear system solution
             *   - ...
             */
            virtual void declare_parameters(ParameterHandler& param);

            /**
             * Function called by optimization loop in order to set the values in the
             * ParameterHandler to the new design parameters.
             * Since ParameterHandler depends on the problem we are solving, set_parameters() is set
             * at the most inner loop of the application.
             */
            virtual void set_parameters(const std::vector<std::string>& name_dvar,
                                        const std::vector<double>& value_dvar,
                                        ParameterHandler& param){};

            /**
             * Set up how many equations are needed and
             * read in parameters for the parameter handler in order to initialize data
             */
            void _initialize(ParameterHandler& param);

            /**
             * Call the other initialize routines from the inherited classes
             */
            virtual void initialize(ParameterHandler& param);

            /**
             * Initialize nonlinear solution
             */
            virtual void initialize_solution (FEVector& initial_guess,
                                              std::shared_ptr<Function<dim> > initial_function = std::shared_ptr<Function<dim> >());
            /**
             * Integration of local bilinear form. Here we loop over the quadrature
             * points and over degrees of freedom in order to compute the matrix for the cell
             * This routine depends on the problem at hand and is called by assemble() in DoF_Handler
             * class
             * The matrix to be assembled is:
             \ f[                                                      *
             \begin{array}{l}
             M(i,j).block(0) = \int_{\Omega} a \nabla \phi_i \nabla \phi_j d\Omega + \int_{\Omega} \phi_i \frac{\partial f}{\partial u_0}\Big|_n \phi_j d\Omega \\
             M(i,j).block(1) = \int_{\Omega} \phi_i \frac{\partial f}{\partial u_1}\Big|_n \phi_j d\Omega \\
             M(i,j).block(2) = \int_{\Omega} \phi_i \frac{\partial f}{\partial u_2}\Big|_n \phi_j d\Omega
             \end{array}
             \f]
             */
            virtual void cell_matrix(MatrixVector& cell_matrices,
                                     const typename DoFApplication<dim>::CellInfo& cell);
            /**
             * Integration of the rhs of the equations. Here we loop over the quadrature
             * points and over degrees of freedom in order to compute the right
             * hand side for each cell
             * This routine depends on the problem at hand and is called by residual() in DoF_Handler
             * class
             * @note This function is called residual because in the case of nonlinear systems
             * the rhs is equivalent to the residual
             */
            virtual void cell_residual(FuelCell::ApplicationCore::FEVector& cell_vector,
                                       const typename DoFApplication<dim>::CellInfo& cell);


            /**
             * Member function used to set dirichlet boundary conditions.
             * This function is application specific and it only computes the boundary_value
             * values that are used to constraint the linear system of equations that is being
             * solved
             */
            virtual void dirichlet_bc(std::map<unsigned int, double>& boundary_values) const;

            /**
             * This class is called by responses to make sure that all responses requested are
             * implemented in either cell_responses, global_responses or face_responses.
             * @note Every time we add a new response that we can calculate we need to update this
             * file.
             */
            virtual void check_responses();

            /**
             * Compute the value of all objective function and constraints
             */
            virtual void cell_responses (std::vector<double>& resp,
                                         const typename DoFApplication<dim>::CellInfo& info,
                                         const FuelCell::ApplicationCore::FEVector& sol);
            /**
             * This class is used to evaluate all responses that do not require looping over cells.
             * An example of one of this types of constraints is the solid volume fraction.
             */
            virtual void global_responses (std::vector<double>& resp,
                                           const FuelCell::ApplicationCore::FEVector& sol);

            /**
             * This class is used to evaluate the derivative of all the functionals that require looping over cells
             * with respect to the design variables.
             * This class is called by responses to evaluate the response at each cell.
             */
            virtual void cell_dresponses_dl(std::vector<std::vector<double> >& /*cell_df_dl*/,
                                            const typename DoFApplication<dim>::CellInfo& /*info*/,
                                            const FuelCell::ApplicationCore::FEVector& /*sol*/) {};

            /**
             * This class is used to evaluate the sensitivities of all responses that do not require looping over cells
             * with respect of the design variables.
             * An example of one of this types of constraints is the solid volume fraction.
             */
            virtual void global_dresponses_dl(std::vector<std::vector<double> >& df_dl,
                                              const FuelCell::ApplicationCore::FEVector& sol);
            /**
             * This class is used to evaluate the derivative of all the functionals that require looping over cells
             * with respect of the unknowns of the system of governing equations.
             * This class is called by responses to evaluate the response at each cell.
             */
            virtual void cell_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >& /*cell_df_du*/,
                                            const typename DoFApplication<dim>::CellInfo& /*info*/,
                                            std::vector<std::vector<double> >& /*src*/) {};

            /**
             * This class is used to evaluate the sensitivities of all responses that do not require looping over cells
             * with respecto of the unknowns of the system of governing equations.
             * An example of one of this types of constraints is the solid volume fraction.
             */
            virtual void global_dresponses_du(std::vector<FuelCell::ApplicationCore::FEVector >& df_du,
                                              const FuelCell::ApplicationCore::FEVector& src);


            /**
             * Post-processing. Evaluate a functional such as the objective function of an
             * optimization problem
             */
            virtual double evaluate (const FuelCell::ApplicationCore::FEVectors& src);

            /**
             * Reimplementation of the routine in the base class BaseApplication in namespace AppFrame so
             * that the right labels are outputed and so that I can compute and output the source term.
             */
            virtual void data_out(const std::string &basename,
                                  const FuelCell::ApplicationCore::FEVectors &src);

 
            /**
             *
             */
            void data_out(std::string basename, const FuelCell::ApplicationCore::FEVectors vectors, std::vector< std::string > solution_names);

        protected:

            ///@name Other internal data
            //@{
            /**
             * Structure where we store the problem we want to solve.
             * Each vector component contains a string with the
             * name of the equation we want to solve
             * Then, the number of components is equation_names.size()
             */
            std::vector<std::string> equation_names;
            /**
             * Structure where we store the name of each component
             * in our problem. The component names are stored in the same
             * way as they are stored in the solution.
             */
            std::vector<std::string> component_names;
            //@}

            ///@name Pre-processor and operating condition classes
            //@{

            /** Initial operating conditions class*/
            OperatingConditions OC;
            //@}

            ///@name Material classes
            //@{
            /**
             * Object used to calculate the properties of the electrolyte in the catalyst layer.
             * In this case we assume is Nafion.
             */
            //FuelCellShop::Material::Nafion electrolyte;

            /**
             * Object used to calculate the carbon black conductivity in the catalyst layer.
             */
            //FuelCellShop::Material::CarbonBlack catalyst_support;

            /** The catalyst object will contain the relevent parameters for the kinetics
             * class, in this we are using a platinum catalyst.
             */
            //FuelCellShop::Material::Platinum catalyst;

            /** The cathode contains water vapour, so we need to create an object water
             * in order to compute viscosity, density, etc. for water
             */
            FuelCellShop::Material::WaterVapor water;

            /** The cathode contains water vapour, so we need to create an object water
             * in order to compute viscosity, density, etc. for water
             */
            FuelCellShop::Material::Oxygen oxygen;

            /** The cathode contains water vapour, so we need to create an object water
             * in order to compute viscosity, density, etc. for water
             */
            FuelCellShop::Material::Nitrogen nitrogen;

            /** The anode contains hydrogen, so we need to create an object water
             * in order to compute viscosity, density, etc. for water
             */
            FuelCellShop::Material::Hydrogen hydrogen; ///< Hydrogen properties

            //@}

            ///@name Layer classes
            //@{
            /** The object AGDL layer will contain all the information relevant to the
             * the anode GDL. We can request any effective property from this class
             */
            boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > AGDL;

            /** The object AMPL layer will contain all the information relevant to the
             * the anode micro-porous layer. We can request any effective property from this class
             */
            boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > AMPL;

            /** The object ACL layer will contain all the information relevant to the
             * the anode catalyst layer. We can request any effective property from this class
             */
            boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > ACL;

            /** The object PEM layer will contain all the information relevant to the
             * the polymer electrolyte membrane. We can request any effective property from this class
             */
            boost::shared_ptr<FuelCellShop::Layer::MembraneLayer<dim> > ML;


            /** The object CCL layer will contain all the information relevant to the
             * the catalyst layer. We can request any effective property from this class
             */
            boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > CCL;

            /** The object CMPL layer will contain all the information relevant to the
             * the cathode micro-porous layer. We can request any effective property from this class
             */
            boost::shared_ptr<FuelCellShop::Layer::MicroPorousLayer<dim> > CMPL;

            /** The object CGDL layer will contain all the information relevant to the
             * the cathode GDL. We can request any effective property from this class
             */
            boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > CGDL;
            //@}
            
            ///@name Post-processing classes (Functional evaluation)
            //@{
            /** Post-processing object to compute the ORR current density */
            FuelCellShop::PostProcessing::ORRCurrentDensityResponse<dim> ORRCurrent;
            /** Post-processing object to compute the HOR current density */
            FuelCellShop::PostProcessing::HORCurrentDensityResponse<dim> HORCurrent;
            /** Post-processing object to compute the water sorption in the CL */
            FuelCellShop::PostProcessing::WaterSorptionResponse<dim> WaterSorption;
            //@}

            ///@name Physics Equations
            //@{
           /**
            * ProtonTransportEquation object
            */
            FuelCellShop::Equation::ProtonTransportEquation<dim> proton_transport;

            /**
            * LambdaTransportEquation object
            */
            FuelCellShop::Equation::LambdaTransportEquation<dim> lambda_transport;

            /**
            * ElectronTransportEquation object
            */
            FuelCellShop::Equation::ElectronTransportEquation<dim> electron_transport;

            FuelCellShop::Equation::NewFicksTransportEquation<dim> ficks_oxygen_nitrogen;
            FuelCellShop::Equation::NewFicksTransportEquation<dim> ficks_water_nitrogen;
            FuelCellShop::Equation::NewFicksTransportEquation<dim> ficks_water_hydrogen;

            /**
            * ReactionSourceTerms object
            */
            FuelCellShop::Equation::ReactionSourceTerms<dim> reaction_source_terms;

            /**
            * SorptionSourceTerms object
            */
            FuelCellShop::Equation::SorptionSourceTerms<dim> sorption_source_terms;
            //@}

            /** Stores the design variable names so that the name can be appended to the .vtk file name. */
            std::vector<std::string> design_var;

            /** Stores the values of the design variables so that the number can be appended to the .vtk file name. */
            std::vector<double> design_var_value;

        private:

            double l_channel; ///< Width of the channel.

            double l_land; ///< Width of the landing.
            
            /**
             * Function to modify the default values of the data file in order to make sure that the equations
             * match those needed in the application.
             * 
             */
            void set_default_parameters_for_application(ParameterHandler &param)
            {
                param.enter_subsection("System management");
                {
                    param.set("Number of solution variables","5");
                    param.enter_subsection("Solution variables");
                    {
                        param.set("Solution variable 1","oxygen_molar_fraction");
                        param.set("Solution variable 2","water_molar_fraction");
                        param.set("Solution variable 3","protonic_electrical_potential");
                        param.set("Solution variable 4","electronic_electrical_potential");
                        param.set("Solution variable 5","membrane_water_content");
                    }
                    param.leave_subsection();
                    
                    param.enter_subsection("Equations");
                    {
                        param.set("Equation 1","Ficks Transport Equation - oxygen");
                        param.set("Equation 2","Ficks Transport Equation - water");
                        param.set("Equation 3","Proton Transport Equation");
                        param.set("Equation 4","Electron Transport Equation");
                        param.set("Equation 5","Membrane Water Content Transport Equation");
                    }
                    param.leave_subsection();
                }
                param.leave_subsection();
                param.enter_subsection("Discretization");
                {
                    param.set("Element","FESystem[FE_Q(1)^5]");
                }
                param.leave_subsection();
            }

        };

    }

}

#endif //_FUELCELL__APPPEMFC_H
