// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: equation_base.h
// - Description: This is a base class for all available openFCST equations
// - Developers: Valentin N. Zingan,    University of Alberta
//               Marc Secanell Gallart, University of Alberta
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_EQUATION_BASE_H_
#define _FCST_FUELCELLSHOP_EQUATION_EQUATION_BASE_H_

#include "system_management.h"
#include "dof_application.h"
#include "initial_and_boundary_data.h"
#include "base_layer.h"
#include "fem_extras.h"
#include "fcst_utilities.h"
#include "fcst_constants.h"

using namespace dealii;
using namespace FuelCell::ApplicationCore;

namespace FuelCellShop
{
namespace Equation
{

///@name Helper classes, structures, and functions
//@{

  /**
   * This simple structure
   * describes a boundary type
   * of a derived equation class.
   *
   * For example
   *
   * @code
   * set Impermeable walls = 0: no-slip
   * @endcode
   *
   * Here
   *
   * - \p boundary_name      = \p Impermeable \p walls,
   * - \p boundary_id        = \p 0,
   * - \p boundary_condition = \p no-slip.
   *
   * \note Other \p boundary_condition like
   * \p Navier \p slip or \p perfect \p slip
   * can be assigned instead.
   *
   * \author Valentin N. Zingan, 2013
   */
  struct BoundaryType
  {
    /**
     * Boundary name.
     */
    std::string  boundary_name;

    /**
     * Boundary indicator.
     */
    unsigned int boundary_id;

    /**
     * Type of boundary condition.
     */
    std::string  boundary_condition;
  };

  /**
   * This simple structure
   * describes an output type
   * of a derived equation class.
   *
   * \author Valentin N. Zingan, 2013
   */
  struct OutputType
  {
    /**
     * Variable name.
     */
    std::string variable_name;

    /**
     * Variable interpretation.
     *
     * There are two options currently implemented in FCST.
     * These options are
     *
     * - \p scalar,
     * - \p vector.
     */
    std::string variable_interpretation;
  };

  // ############################################################################################## Bhaiya #######

  /**
   * This simple structure stores certain information
   * regarding a particular variable for the equation (all of them retrieved from #SystemManagement).
   *
   * The purpose of defining such a structure is to avoid
   * repeated use of #SystemManagement functions, (which use \p find methods),
   * hence improving code speed.
   *
   * This structure will normally be filled in #make_assemblers_generic_constant_data of
   * derived classes.
   *
   * \author Madhur Bhaiya, 2013
   */
  struct VariableInfo
  {
    /**
     * Index of the user-defined solution variable, retrieved from #SystemManagement.
     */
    unsigned int solution_index;

    /**
     * Block index of the matrix relating to the variable corresponding to an equation, retrieved from #SystemManagement.
     */
    unsigned int block_index;

    /**
     * Index corresponding to type of \p fevalue object used for this variable.
     * This information is also retrieved from #SystemManagement, from a group of actual \p FEVALUES objects (in \b AppFrame).
     */
    unsigned int fetype_index;

    /**
     * Boolean storing whether indices exist or not.
     * This boolean member is very useful for checks inside the derived equation classes, to
     * avoid errors.
     */
    bool indices_exist;
  };

  // ############################################################################################## Bhaiya #######

//@}

///@name Exceptions
//@{

  // ############################################################################################## Bhaiya #######

  /**
   * Exception thrown when a particular variable required by the equation class,
   * does not exist in the user defined solution variables.
   */
  DeclException2(VariableShouldExistForEquation,
                 std::string,
                 std::string,
                 << "The user-defined variable with name \"" << arg1 << "\" should be one of the solution variables for equation with name \"" << arg2 << "\"");

  /**
   * Exception thrown when index of the
   * user defined equation do not match with
   * the index of the base variable for the equation.
   */
  DeclException2(IndexDoNotMatch,
                 std::string,
                 std::string,
                 << "The index of variable \"" << arg1 << "\" do not match with the index of equation \"" << arg2 << "\"");

  // ############################################################################################## Bhaiya #######

//@}

/**
 * This class contains
 * generic data and methods
 * heavily used by all derived
 * equation classes.
 *
 * This class is used to lock
 * the interface of all derived
 * equation classes under
 * the FuelCellShop::Equation namespace.
 *
 * The functionality of
 * this class can be extended
 * if needed.
 *
 * \author Valentin N. Zingan, 2012
 * \author Marc Secanell Gallart, 2012
 */

template<int dim>
class EquationBase : public Subscriptor
{
public:

///@name Local CG FEM based assemblers
//@{

  /**
   * Assemble local cell matrix.
   */
  virtual void assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                    const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                    FuelCellShop::Layer::BaseLayer<dim>* const              layer)
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * Assemble local cell residual.
   */
  virtual void assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                      const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                      FuelCellShop::Layer::BaseLayer<dim>* const              layer)
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * Assemble local boundary matrix.
   */
  virtual void assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                    const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                    FuelCellShop::Layer::BaseLayer<dim>* const              layer)
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * Assemble local boundary residual.
   */
  virtual void assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_residual,
                                      const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                      FuelCellShop::Layer::BaseLayer<dim>* const              layer)
  {
    print_caller_name(__FUNCTION__);
  }

//@}

///@name Accessors and info
//@{

  /**
   * This function returns
   * \p internal_cell_couplings
   * of a derived equation class.
   */
  const couplings_map& get_internal_cell_couplings() const
  {
    AssertThrow( !internal_cell_couplings.empty() , ExcInternalError() );
    return internal_cell_couplings;
  }

  /**
   * This function returns
   * \p internal_flux_couplings (DG FEM only)
   * of a derived equation class.
   */
  const couplings_map& get_internal_flux_couplings() const
  {
    AssertThrow( !internal_flux_couplings.empty() , ExcInternalError() );
    return internal_flux_couplings;
  }

  /**
   * This function returns
   * \p component_materialID_value
   * of a derived equation class.
   */
  const component_materialID_value_map& get_component_materialID_value() const
  {
    AssertThrow( !component_materialID_value.empty() , ExcInternalError() );
    return component_materialID_value;
  }

  /**
   * This function returns
   * \p component_boundaryID_value
   * of a derived equation class.
   */
  const component_boundaryID_value_map& get_component_boundaryID_value() const
  {
    AssertThrow( !component_boundaryID_value.empty() , ExcInternalError() );
    return component_boundaryID_value;
  }

  /**
   * This function returns
   * \p boundary_types
   * of a derived equation class.
   */
  const std::vector< BoundaryType >& get_boundary_types() const
  {
    return boundary_types;
  }

  /**
   * This function returns
   * \p multi_boundary_types
   * of a derived equation class.
   */
  const std::vector< std::vector< BoundaryType > >& get_multi_boundary_types() const
  {
    return multi_boundary_types;
  }

  /**
   * This function returns
   * \p output_types
   * of a derived equation class.
   */
  const std::vector< OutputType >& get_output_types() const
  {
    return output_types;
  }

  /**
   * This function returns
   * \p multi_output_types
   * of a derived equation class.
   */
  const std::vector< std::vector< OutputType > >& get_multi_output_types() const
  {
    return multi_output_types;
  }

  /**
   * This function returns
   * \p equation_name
   * of a derived equation class.
   */
  const std::string& get_equation_name() const
  {
    return equation_name;
  }

  /**
   * This function returns
   * \p matrix_block_indices
   * of a derived equation class.
   */
  const std::vector<unsigned int>& get_matrix_block_indices() const
  {
    return matrix_block_indices;
  }

  /**
   * This function returns
   * \p residual_indices
   * of a derived equation class.
   */
  const std::vector<unsigned int>& get_residual_indices() const
  {
    return residual_indices;
  }

  /**
   * This function prints out
   * the equations info
   * of a derived equation class.
   */
  virtual void print_equation_info() const
  {
    print_caller_name(__FUNCTION__);
  }

//@}

  /**
   * @p true,
   * if variable initial data
   * is prescribed on a part of the domain.
   * @p false, otherwise.
   */
  bool variable_initial_data;

  /**
   * @p true,
   * if variable Dirichlet boundary conditions
   * are prescribed on a part of the boundary.
   * @p false, otherwise.
   */
  bool variable_boundary_data;

protected:

///@name Constructors, destructor, and initialization
//@{

  /**
   * Constructor.
   *
   * @warning The constructor of the child classes needs to define the following to variables:
   * - name_base_variable
   * - equation_name
   * for example, for lambda_transport_equation
   * @code
   * this->equation_name = "Membrane Water Content Transport Equation";
   * this->name_base_variable = "membrane_water_content";
   * @endcode
   *
   * These two variables are used to define the section where the information is stored in the parameter
   * file for the equations.
   *
   */
  EquationBase(FuelCell::SystemManagement& sys_management);

  /**
   * Destructor.
   */
  virtual ~EquationBase();

  /**
   * Declare parameters.
   */
  virtual void declare_parameters(ParameterHandler& param) const;

  /**
   * Initialize parameters.
   */
  virtual void initialize(ParameterHandler& param);

  /**
   * Set parameters using the parameter file,
   * in order to run parametric/optimization studies.
   */
  virtual void set_parameters(const std::vector<std::string>& name_dvar,
                              const std::vector<double>&      value_dvar,
                              ParameterHandler&               param)
  {
    print_caller_name(__FUNCTION__);
  }

//@}

///@name Local CG FEM based assemblers - make_ functions
//@{

  /**
   * This function is overridden
   * in the derived equation classes.
   */
  virtual void make_assemblers_generic_constant_data()
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function is overridden
   * in the derived equation classes.
   */
  virtual void make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function is overridden
   * in the derived equation classes.
   */
  virtual void make_assemblers_bdry_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info)
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function is overridden
   * in the derived equation classes.
   */
  virtual void make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                  FuelCellShop::Layer::BaseLayer<dim>* const              layer)
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function is overridden
   * in the derived equation classes.
   */
  virtual void make_assemblers_bdry_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                  FuelCellShop::Layer::BaseLayer<dim>* const              layer)
  {
    print_caller_name(__FUNCTION__);
  }

//@}

///@name Other - make_ functions
//@{

  /**
   * This function fills out
   * \p internal_cell_couplings
   * of a derived equation class.
   */
  virtual void make_internal_cell_couplings()
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function fills out
   * \p internal_flux_couplings (DG FEM only)
   * of a derived equation class.
   */
  virtual void make_internal_flux_couplings()
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function fills out
   * \p component_materialID_value
   * of a derived equation class.
   */
  virtual void make_component_materialID_value()
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function fills out
   * \p component_boundaryID_value
   * of a derived equation class.
   */
  virtual void make_component_boundaryID_value()
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function fills out
   * \p boundary_types
   * of a derived equation class.
   */
  virtual void make_boundary_types()
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function fills out
   * \p multi_boundary_types
   * of a derived equation class.
   */
  virtual void make_multi_boundary_types()
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function fills out
   * \p output_types
   * of a derived equation class.
   */
  virtual void make_output_types()
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function fills out
   * \p multi_output_types
   * of a derived equation class.
   */
  virtual void make_multi_output_types()
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function fills out
   * \p matrix_block_indices
   * of a derived equation class.
   */
  virtual void make_matrix_block_indices()
  {
    print_caller_name(__FUNCTION__);
  }

  /**
   * This function fills out
   * \p residual_indices
   * of a derived equation class.
   */
  virtual void make_residual_indices()
  {
    print_caller_name(__FUNCTION__);
  }

//@}

///@name Converters
//@{

  /**
   * This function changes the order of
   * dealii::FullMatrix<double> \p target
   * from standard to block-wise.
   */
  void standard_to_block_wise(FullMatrix<double>& target) const;

  /**
   * This function changes the order of
   * dealii::Vector<double> \p target
   * from standard to block-wise.
   */
  void standard_to_block_wise(Vector<double>& target) const;

  /**
   * This function converts
   * the standard ordered structure
   * dealii::FullMatrix<double> \p src
   * into
   * the block-wise ordered structure
   * FuelCell::ApplicationCore::MatrixVector \p dst.
   *
   * \p matrix_block_indices is
   * the system matrix block
   * indices (a derived equation class) drawn from the global
   * structure (a derived equation class + other active equation classes included into the computation).
   */
  void dealII_to_appframe(FuelCell::ApplicationCore::MatrixVector&          dst,
                          const FullMatrix<double>&        src,
                          const std::vector<unsigned int>& matrix_block_indices) const;

  /**
   * This function converts
   * the standard ordered structure
   * dealii::Vector<double> \p src
   * into
   * the block-wise ordered structure
   * FuelCell::ApplicationCore::FEVector \p dst.
   *
   * \p residual_indices is
   * the residual
   * indices (a derived equation class) drawn from the global
   * structure (a derived equation class + other active equation classes included into the computation).
   */
  void dealII_to_appframe(FuelCell::ApplicationCore::FEVector&              dst,
                          const Vector<double>&            src,
                          const std::vector<unsigned int>& residual_indices) const;

//@}

///@name Other functions
//@{

  /**
   * This function returns \p true
   * if a boundary indicator of an external face on the triangulation
   * coincides with a boundary indicator defined in
   * the parameters file of a derived equation class.
   */
  bool belongs_to_boundary(const unsigned int& tria_boundary_id,
                           const unsigned int& param_boundary_id) const
  {
    return tria_boundary_id == param_boundary_id;
  }

//@}

///@name Minor functions
//@{

  /**
   * This function is used
   * to print out the name of another function
   * that has been declared in the scope of this class,
   * but not yet been implemented.
   */
  void print_caller_name(const std::string& caller_name) const;

//@}

  //////////
  // DATA //
  //////////

///@name Local CG FEM based assemblers - constant data (generic)
//@{

  /**
   * Number of degrees of freedom
   * per cell.
   */
  unsigned int dofs_per_cell;

//@}

///@name Local CG FEM based assemblers - constant data (cell)
//@{

  /**
   * Number of quadrature points
   * per cell.
   */
  unsigned int n_q_points_cell;

//@}

///@name Local CG FEM based assemblers - constant data (boundary)
//@{

  /**
   * Number of quadrature points
   * per boundary.
   */
  unsigned int n_q_points_bdry;

//@}

///@name Local CG FEM based assemblers - variable data (active mesh iterators)
//@{

  /**
   * Currently active DoFHandler<dim> active cell iterator.
   */
  typename DoFHandler<dim>::active_cell_iterator cell;

  /**
   * Currently active DoFHandler<dim> active boundary iterator.
   */
  typename DoFHandler<dim>::active_face_iterator bdry;

//@}

///@name Local CG FEM based assemblers - variable data (previous Newton iteration - cell)
//@{

  /**
   * Implementation is
   * in the derived equation classes.
   */

//@}

///@name Local CG FEM based assemblers - variable data (current Newton iteration - cell)
//@{

  /**
   * Implementation is
   * in the derived equation classes.
   */

//@}

///@name Local CG FEM based assemblers - variable data (other - cell)
//@{

  /**
   * Jacobian of mapping by Weight
   * in the quadrature points of a cell.
   */
  std::vector<double> JxW_cell;

//@}

///@name Local CG FEM based assemblers - variable data (previous Newton iteration - boundary)
//@{

  /**
   * Implementation is
   * in the derived equation classes.
   */

//@}

///@name Local CG FEM based assemblers - variable data (current Newton iteration - boundary)
//@{

  /**
   * Implementation is
   * in the derived equation classes.
   */

//@}

///@name Local CG FEM based assemblers - variable data (other - boundary)
//@{

  /**
   * Jacobian of mapping by Weight
   * in the quadrature points of a boundary.
   */
  std::vector<double> JxW_bdry;

  /**
   * Normal vectors
   * in the quadrature points of a boundary.
   */
  std::vector< Point<dim> > normal_vectors;

  /**
   * Tangential vectors
   * in the quadrature points of a boundary.
   */
  std::vector< std::vector< Point<dim> > > tangential_vectors;

//@}

///@name GENERIC DATA
//@{

  /**
   * Pointer to the external
   * YourApplication<dim>::system_management object.
   */
  FuelCell::SystemManagement* system_management;

  /**
   * This object contains the info on
   * how the equations and solution variables
   * of a derived equation class are coupled.
   */
  couplings_map internal_cell_couplings;

  /**
   * This object contains the info on
   * how the "X" and "Y"
   * of a derived equation class are coupled (DG FEM only).
   */
  couplings_map internal_flux_couplings;

  /**
   * This object reflects the following
   * structure (see FuelCell::InitialAndBoundaryData namespace docs):
   *
   * - \p first  \p argument : name of the solution component,
   * - \p second \p argument : material id,
   * - \p third  \p argument : value of the solution component.
   */
  component_materialID_value_map component_materialID_value;

  /**
   * This object reflects the following
   * structure (see FuelCell::InitialAndBoundaryData namespace docs):
   *
   * - \p first  \p argument : name of the solution component,
   * - \p second \p argument : boundary id,
   * - \p third  \p argument : value of the solution component.
   */
  component_boundaryID_value_map component_boundaryID_value;

  /**
   * The list of
   * boundary types
   * of a derived equation class.
   */
  std::vector< BoundaryType > boundary_types;

  /**
   * The list of
   * multiple boundary types
   * of a derived equation class.
   */
  std::vector< std::vector< BoundaryType > > multi_boundary_types;

  /**
   * The list of
   * output types
   * of a derived equation class.
   */
  std::vector< OutputType > output_types;

  /**
   * The list of
   * multiple output types
   * of a derived equation class.
   */
  std::vector< std::vector< OutputType > > multi_output_types;

  /**
   * The name
   * of a derived equation class.
   */
  std::string equation_name;

  /**
   * Const std::string member storing name of the
   * base solution variable corresponding to
   * the equation represented by this class.
   */
  std::string name_base_variable;

  /**
   * The system matrix block
   * indices (a derived equation class) drawn from the global
   * structure (a derived equation class + other active equation classes included into the computation).
   */
  std::vector<unsigned int> matrix_block_indices;

  /**
   * The residual
   * indices (a derived equation class) drawn from the global
   * structure (a derived equation class + other active equation classes included into the computation).
   */
  std::vector<unsigned int> residual_indices;

  /**
   * This vector contains
   * the collection of internal "counters"
   * used by the derived equation classes.
   */
  std::vector<bool> counter;

//@}

};

} // Equation

} // FuelCellShop

#endif