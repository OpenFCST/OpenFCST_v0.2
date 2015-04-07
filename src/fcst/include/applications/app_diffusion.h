//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: app_diffusion.cc
//    - Description:
//    - Developers: M. Sabharwal, M. Secanell
//    - $Id: app_diffusion.h 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL__APP_DIFFUSION__H
#define _FUELCELL__APP_DIFFUSION__H

// Include deal.II classes
#include "base/parameter_handler.h"
#include "base/function_lib.h"
#include "base/function.h"
#include "base/quadrature_lib.h"

#include "lac/block_vector.h"
#include "lac/full_matrix.h"
#include "lac/solver_cg.h"
#include "lac/solver_gmres.h"
#include "lac/precondition.h"
#include "lac/precondition_block.h"
#include "lac/block_matrix_array.h"
#include "lac/sparse_ilu.h"
#include "lac/sparse_direct.h"

#include "grid/grid_generator.h"
#include "grid/tria_accessor.h"
#include "grid/tria_iterator.h"
#include "grid/tria_boundary_lib.h"

#include "fe/fe_values.h"

#include "numerics/vector_tools.h"
#include "numerics/matrix_tools.h"

#include "boost/shared_ptr.hpp"

// Include appframe classes
#include "block_matrix_application.h"

// Include FuelCell classes
#include "fcst_constants.h"
#include "application_core/optimization_block_matrix_application.h"
#include "solver_utils.h"

#include "operating_conditions.h"
#include "geometries.h"

#include "PureGas.h"
#include "GasMixture.h"
#include "platinum.h"
#include "nafion.h"
#include "carbon.h"

#include "gas_diffusion_layer.h"


#include "new_ficks_transport_equation.h"


#include "linear_solvers.h"

#include "data_out.h"


//Include STL
#include "fstream"
#include "iostream"
#include "sstream"
#include <algorithm>

// Use namespace of deal.II
using namespace dealii;

namespace FuelCell
{

    namespace Application
    {
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        /**
         * The application can be used to simulate a gas flow through a porous media. 
         * This class is used to solve a diffusion problem using Fick's Law. 
         *
         * The equation solved is written as follows:
         * \f[
         * R_1(\vec{u}) = \nabla \cdot (c_{total}D^{eff}_{O_2} \nabla x_{O_2} ) = 0
         * \f]
         * The governing equation in the weak form is actually linear and can be solved directly using a linear solver like UMFPACK or GMRES.
         * The main file for this application should look like this:
         * @code
         * subsection Simulator
         *  set simulator name = diffusion
         *  set simulator parameter file name = data.prm
         *  set solver name = Linear 
         *  set Dakota direct = false 
         * end
         * @endcode
         * NOTE: The solver name is set to Linear instead of the Newton solvers
         * The lines below indicate how to setup your boundary and initial conditions.
         * [Although initial conditions are not required for this problem it is advisable to use them just to see that the solver is working.]
         *
         * @code
         * subsection Equations
         * 
         * subsection Ficks Transport Equation - oxygen
         * 
         *  subsection Initial data
         *   set oxygen_molar_fraction = 0: 100 # where 0 indicates the material_id setup in the grid and 100 is the concentration of solute in mol/cm^3
         *  end
         * 
         *  subsection Boundary data
         *   set oxygen_molar_fraction = 1: 0.4, 2:0.01 #where 1 & 2 denote the boundary ids and 0.4 and 0.01 are the concentrations of solute in mol/cm^3 at the respective boundary
         *  end
         * 
         * end
         * @endcode
         *
         * @author Mayank Sabharwal & Marc Secanell, 2014
         *
         */
        template<int dim>
        class AppDiffusion : public OptimizationBlockMatrixApplication<dim>
        {
        public:
            
            ///@name Constructors, destructor, and initialization
            //@{
            
            /**
             * Constructor.
             */
            AppDiffusion( boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data = 
            boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData >() );
            
            /**
             * Destructor.
             */
            ~AppDiffusion();
            
            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param);
            
            /**
             * Initialize parameters.
             */
            virtual void initialize(ParameterHandler& param);
            
            /**
             * The initial guess along with the appropriate BCs is formed here.
             */
              virtual void initialize_solution (FEVector& initial_guess,
                                                std::shared_ptr<Function<dim> > initial_function = std::shared_ptr<Function<dim> >());

            //@}
            
            ///@name Local CG FEM based assemblers
            //@{
            
            /**
             * Assemble local cell matrix.
             */
            virtual void cell_matrix(MatrixVector&                                 cell_matrices,
                                     const typename DoFApplication<dim>::CellInfo& cell_info);
            
            /**
             * Assemble local cell residual.
             */
            virtual void cell_residual(FuelCell::ApplicationCore::FEVector&          cell_res,
                                       const typename DoFApplication<dim>::CellInfo& cell_info);
            
            //@}
            
            ///@name Other functions
            //@{
            
            /**
             * BCs.
             */
            virtual void dirichlet_bc(std::map<unsigned int, double>& boundary_values) const;

            /**
             * Evaluate the response functions: Not implemented in this class; Returns 0.0
             */
            virtual double evaluate (const FuelCell::ApplicationCore::FEVectors& src);
            
            /**
             * Output results.
             */
            virtual void data_out(const std::string&         filename,
                                  const FEVectors& src);
            
            //@}
            
            ///@name Post-processing
            //@{
            
            /**
             * Compute some functionals.
             */
            virtual void bdry_responses(std::vector<double>&                                                     dst,
                                        const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                        const FuelCell::ApplicationCore::FEVector& src);
            //@}
            
        protected:            
            ///@name System matrix and boundary condition objects
            //@{

            /**
             *
             * component_boundaryID_value_maps info:
             *
             * subsection Initial data
             * set oxygen_molar_fraction = 0: 100
             * end
             * subsection Boundary data
             * set oxygen_molar_fraction = 1: 0.4, 2:0.01
             * end
             * end
             */
            
            //@}
            
            ///@name Gases and operating conditions data
            //@{       
            /**
             * Operating conditions.
             */
            FuelCell::OperatingConditions OC;
            
            /**
             * Solute.
             */
            boost::shared_ptr<FuelCellShop::Material::PureGas> solute;
            
            
            
            /**
             * Solvent.
             */
            boost::shared_ptr<FuelCellShop::Material::PureGas> solvent;
            
                        
            //@}
            
            ///@name Layer objects
            //@{
            /**
             * Cathode GDL.
             */
            boost::shared_ptr< FuelCellShop::Layer::GasDiffusionLayer<dim> > CGDL;
            
            //@}
            
            ///@name Equation objects
            //@{
            /**
             * This object describes
             * the equations that we are going to
             * solve here.
             */
            FuelCellShop::Equation::NewFicksTransportEquation<dim> ficks_transport_equation;
            
            /**
             * This object describes
             * the equations that we are going to
             * solve here.
             */
            

        private:
            /**
             * Compute some functionals that are not needed for most applications (this section is not necessary in most cases.)
             */
            
            
        };
        
    } // Application
    
} // FuelCell

#endif
