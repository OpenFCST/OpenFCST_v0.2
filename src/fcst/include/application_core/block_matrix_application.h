// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2009 by Guido Kanschat
// Copyright (C) 2006-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: block_matrix_application.h
// - Description:
// - Developers: Guido Kanschat, Texas A&M University
//               Marc Secanell,  University of Alberta
// - Id: $Id: block_matrix_application.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FUEL_CELL_APPLICATION_CORE_BLOCK_MATRIX_APPLICATION_H_
#define _FUEL_CELL_APPLICATION_CORE_BLOCK_MATRIX_APPLICATION_H_

#include <base/quadrature.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_cg.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <dofs/dof_tools.h>
#include <numerics/matrix_tools.h>
#include <base/parameter_handler.h>
#include <base/quadrature_lib.h>
#include <dofs/dof_handler.h>
#include <fe/fe.h>
#include <fe/fe_values.h>


#include <lac/petsc_solver.h>
#include <lac/petsc_parallel_block_sparse_matrix.h>
#include <lac/petsc_precondition.h>

#include "system_management.h"
#include "dof_application.h"
#include "matrix_block.h"
#include "linear_solvers.h"

namespace FuelCell
{
    namespace ApplicationCore
    {
        /**
         * Application handling matrices and assembling linear systems of
         * equations.
         *
         * This class inherits the information on Triangulation and DoFHandler
         * from its base class DoFApplication. It adds information on the
         * linear system, in particular the #sparsities.
         *
         * Additionally, it holds a BlockSparseMatrix as the system matrix of
         * the finite element problem. Generic loops over mesh cells and faces
         * are implemented in order to assemble() the system matrices and
         * compute the residual().  These functions use the virtual local
         * integration functions declared here and implemented by the derived
         * actual application class.
         *
         * @note You probably miss functions for building the right hand side
         * of a linear system. They are missing deliberately, since the linear
         * case is very special and can be handled in the framework of
         * nonlinear problems. So, if you have a linear problem, code the
         * local matrices for residual() so they ignore the BlockVector
         * <tt>src</tt>.
         *
         * <h3>Deriving your own application</h3>
         *
         *
         * @author Guido Kanschat
         * @author Marc Secanell
         */
        template <int dim>
        class BlockMatrixApplication :
        public DoFApplication<dim>
        {
        public:
            ///@name Constructors and Initialization
            //@{
            /**
             * Constructor for an object,
             * owning its own mesh and dof handler.
             *
             * If <tt>data</tt> is null,
             * a new handler for
             * application data is
             * generated.
             *
             * see constructor of DoFApplication.
             */
            BlockMatrixApplication(boost::shared_ptr<ApplicationData> data =
            boost::shared_ptr<ApplicationData>());
            
            /**
             * Constructor for an object,
             * borrowing mesh and dof
             * handler from another
             * object. If
             * <tt>triangulation_only</tt>
             * is true, only the
             * triangulation is borrowed.
             *
             * see constructor of DoFApplication.
             */
            BlockMatrixApplication(DoFApplication<dim>&,
                                   bool triangulation_only);
            
            /**
             * Initialize data of this
             * class.
             * @note remesh_dofs() has to be called before this member function
             * since remesh_dofs() initializes data that is used in this routine.
             * In particular, dof and block_indices.
             */
            void _initialize(ParameterHandler& param);
            
            /**
             * Declare parameters related
             * to matrices.
             * 
             * In particular the following sections are declared:
             * 
             * @code
             * 
             * 
             * subsection Discretization
             *   subsection Matrix
             *     set Quadrature cell = -1
             *     set Quadrature face = -1
             *   end
             * end
             * #subsection used to select linear solver and linear solver parameters:
             * subsection Linear Application     
             *   set Assemble numerically = false          #Assemble Jacobian using analytical derivatives or numerically.
             * 
             *   set Type of linear solver = MUMPS             
             *            Options: (For parallel code (--with-petsc) MUMPS|CG|Bicgstab )
             *                      For serial code ILU-GMRES|UMFPACK|Bicgstab )
             *   set Max steps = 100
             *   set Tolerance = 1.e-10
             *   set Log history = false
             *   set Log frequency = 1
             *   set Log result = true
             * 
             *   set Allocate additional memory for MUMPS = 
             * end
             * @endcode
             *
             * @note The function declare_parameters of
             * derived applications should always call the one of its base class.
             */
            virtual void declare_parameters(ParameterHandler& param);
            /**
             * Initialize application
             */
            virtual void initialize (ParameterHandler& param);
            //@}
            /**
             * serial and PETSc assemble functions called by
             * assemble().
             */

            #ifdef OPENFCST_WITH_PETSC
            void PETSc_assemble(const FEVectors&);
            #else
            void serial_assemble(const FEVectors&);
            #endif




            /**
             * A call back function for assemble(), called after all
             * cell matrices have been entered, but before the face
             * matrices are computed.
             */
            virtual void post_cell_assemble();
            
            ///@name ApplicationBase interface
            //@{
            /**
             * Initialize sparsity patterns
             * and matrices for the new
             * mesh.
             */
            void remesh_matrices();
            /**
             * Refine grid accordingly to the Refinement options under "Grid Generation"
             * in the parameter file.
             */
            virtual void remesh();
            
            /**
             * Loop over all cells and  assemble the system
             * #matrices by using the local matrix integration functions
             * cell_matrix(), bdry_matrix() and face_matrix() provided
             * by the derived class.
             */
            void assemble(const FEVectors&);
            
            
            /**
             * Compute the Jacobian of the system of equations,
             * \f[
             * J(i,j) = \frac{\partial R_i}{\partial u_j}
             * \f]
             * by using forward differences. The solution is perturbed at each degree of freedom
             * and the residual function is used to compute the residual for each perturbed solution.
             * Forward differences are used to compute the Jacobian.
             *
             * The Jacobian is stored in this->matrix
             *
             * @param[in] src A vector of BlockVectors containting the solution at the previous
             *  Newton step and the residual.
             * @param[in] delta The step size used for numerical differenciation.
             *
             * \note This routine is used in combination with a Netwon nonlinear solver since it
             * requires src to contain a field named "Newton iterate" with the solution at the
             * previous Newton step.
             *
             * \author M. Secanell, 2013
             */
            void assemble_numerically(const FEVectors& src,
                                      const double delta = 1e-6);
            
            /**
             * Redefinition of residual_constraints() in
             * DoFHandler. In this member function, I adjust the residual
             * in order to account for hanging nodes and I apply
             * Dirichlet boundary conditions given by
             * dirichlet_bc().
             */
            void residual_constraints(FEVector& dst) const;

            #ifdef OPENFCST_WITH_PETSC
            void residual_constraints(PETScWrappers::MPI::Vector& dst, const std::vector<unsigned int>& idx_list) const;
            #endif
            
            //@}
            /// @name Local integrators
            //@{
            /**
             * Integration of local
             * bilinear form. Here we loop over the quadrature
             * points and over degrees of freedom in order to compute
             * the matrix for the cell.
             */
            virtual void cell_matrix(MatrixVector& cell_matrices,
                                     const typename DoFApplication<dim>::CellInfo& cell);
            
            /**
             * Integration of local
             * bilinear form.
             */
            virtual void bdry_matrix(MatrixVector& face_matrices,
                                     const typename DoFApplication<dim>::FaceInfo& face);
            
            /**
             * Integration of local
             * bilinear form.
             */
            virtual void face_matrix(MatrixVector& matrices11,
                                     MatrixVector& matrices12,
                                     MatrixVector& matrices21,
                                     MatrixVector& matrices22,
                                     const typename DoFApplication<dim>::FaceInfo& face1,
                                     const typename DoFApplication<dim>::FaceInfo& face2);
            
            /**
             * Member function used to set dirichlet boundary conditions.
             * This function is application specific and it only computes the boundary_value
             * values that are used to constraint the linear system of equations that is being
             * solved
             */
            //@}
            /// @name Dirichlet boundary conditions
            //@{
            virtual void dirichlet_bc(std::map<unsigned int, double>& boundary_values) const;
            //@}
            
            virtual void solve(FEVector&        dst,
                                          const FEVectors& src);

            #ifdef OPENFCST_WITH_PETSC
            void PETSc_solve(FuelCell::ApplicationCore::FEVector system_rhs,
                    FEVector& solution, const FEVectors& src);
            #else
            void serial_solve(FuelCell::ApplicationCore::FEVector system_rhs,
                    FEVector& solution);
            #endif


      protected:
            /**
             * Quadrature rule for matrix
             * assembling on cells.
             */
            boost::shared_ptr<Quadrature<dim> > quadrature_assemble_cell;
            
            /**
             * Quadrature rule for matrix
             * assembling on faces.
             */
            boost::shared_ptr<Quadrature<dim-1> > quadrature_assemble_face;
            
            /**
             * Bool determining whether or not to repair diagonal before solving.
             */
            bool repair_diagonal;


            /**
             * Variable to store boundary values, so they only need to be computed once per mesh refinement
             */
            std::map<unsigned int, double> boundary_values;

        private:

            /**
             * Sparsity patterns.
             */
            BlockSparsityPattern sparsities;
        
        
        protected:        
        ///@name System matrix and boundary conditions
        //@{
            /**
             * Storage for the actual
             * matrix.
             */

            #ifdef OPENFCST_WITH_PETSC
            PETScWrappers::MPI::SparseMatrix matrix;
            #else
            BlockSparseMatrix<double> matrix;
            #endif

            /**
             * Solver control object
             */
            SolverControl solver_control;


            //@}
            
            ///@name Other application variables
            //@{
            /**
             * Variable used to select the appropriate linear solver:
             */
            std::string solver_name;
            /**
             * Variable used to select if assembly should be done analytically or numerically
             */
            bool assemble_numerically_flag;
            /**
             * Variable used for configuring MUMPS to use more memory in solve() function.
             */
            bool mumps_additional_mem;
            //@}
    };
    
}

}

#endif
