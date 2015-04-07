// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: reaction_source_terms.h
// - Description: This class is used to assemble both cell matrix and cell residual
//                for reaction source terms in the catalyst layers for various equation classes
// - Developers: Madhur Bhaiya, Marc Secanell, Valentin N. Zingan
// - Id: $Id$
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_EQUATION_REACTION_SOURCE_TERMS_H_
#define _FCST_FUELCELLSHOP_EQUATION_REACTION_SOURCE_TERMS_H_

#include "utils/fcst_constants.h"

#include "equations/equation_base.h"
#include "equations/reaction_heat.h"

#include "layers/catalyst_layer.h"

#include "reactions/tafel_kinetics.h"
#include "reactions/double_trap_kinetics.h"
#include "reactions/dual_path_kinetics.h"

namespace FuelCellShop
{
    namespace Equation
    {

        ///@name Exceptions
        //@{

        /**
         * Exception thrown when a particular "Required" solution variable
         * is not found for \p cathode or \p anode kinetics,
         * in \p user-defined solution variables.
         */
        DeclException2(VariableNotFoundForKinetics,
                       std::string,
                       std::string,
                       << "For " << arg1 << " kinetics source terms, \"" << arg2 << "\" is not found as one of the solution variables.");

        //@}

        /**
         * This class assembles the reaction source terms for all other transport equations, if there's any.
         *
         * \author Madhur Bhaiya,      2013
         * \author Valentin N. Zingan, 2013 - all couplings with fluid transport equations
         */

        template<int dim>
        class ReactionSourceTerms : public EquationBase<dim>
        {
        public:

            ///@name Constructors, destructor, and initalization
            //@{

            /**
             * Constructor.
             */
            ReactionSourceTerms(FuelCell::SystemManagement& system_management);

            /**
             * Destructor.
             */
            virtual ~ReactionSourceTerms();

            /**
             * Declare parameters.
             */
            virtual void declare_parameters(ParameterHandler& param) const;

            /**
             * Initialize parameters.
             */
            virtual void initialize(ParameterHandler& param);

            /**
             * Set the pointer to cathode kinetics in the object.
             * \note If an application uses cathode kinetics model, this method should be called neccessarily before initializing this object (using #initialize method).
             */
            inline void set_cathode_kinetics(FuelCellShop::Kinetics::BaseKinetics* kinetics)
            {
                cathode_kinetics = kinetics;
            }

            /**
             * Set the pointer to anode kinetics in the object.
             * \note If an application uses anode kinetics model, this method should be called neccessarily before initializing this object (using #initialize method).
             */
            inline void set_anode_kinetics(FuelCellShop::Kinetics::BaseKinetics* kinetics)
            {
                anode_kinetics = kinetics;
            }

            //@}

            ///@name Local CG FEM based assemblers
            //@{

            /**
             * Assemble local cell matrix.
             */
            virtual void assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                              const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                              FuelCellShop::Layer::BaseLayer<dim>* const                               layer);

            /**
             * Assemble local cell residual.
             */
            virtual void assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                                const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                FuelCellShop::Layer::BaseLayer<dim>* const                               layer);

            //@}

            ///@name Accessors & Info
            //@{

            /**
             * This function is used to adjust \p std::vector \p < \p internal_cell_couplings \p >, which is generated after getting \p internal_cell_couplings
             * from all the equations being used in the application.
             *
             * \note It's very important to use this function, if we are considering source terms due to reaction in the catalyst layers. It's also noteworthy that
             * this function should be used after whole vector of \p internal_cell_couplings is created, and \b BEFORE \p make_cell_couplings of \b SystemManagement
             * is called.
             */
            virtual void adjust_internal_cell_couplings( std::vector< couplings_map >& dst,
                                                         const std::vector< FuelCellShop::Material::PureGas* >& gases = std::vector< FuelCellShop::Material::PureGas* >() );

            /**
             * This function prints out
             * the \p info for this class.
             */
            virtual void print_equation_info() const;

            /**
             * Method to get whether the ir-reversible heating due to \p ORR, inside the cathode catalyst layer is \p ON or \p OFF (#irrev_heat_ccl).
             * The typical use of this method is for the \p \b post-processing routines in the application.
             */
            inline bool get_irrev_heat_ccl() const
            {
                return irrev_heat_ccl;
            }

            /**
             * Method to get whether the ir-reversible heating due to \p HOR, inside the anode catalyst layer is \p ON or \p OFF (#irrev_heat_acl).
             * The typical use of this method is for the \p \b post-processing routines in the application.
             */
            inline bool get_irrev_heat_acl() const
            {
                return irrev_heat_acl;
            }

            /**
             * Method to get whether the reversible (entropic) heating due to net reaction forming liquid water product is \p ON or \p OFF (#rev_heat).
             * The typical use of this method is for the \p \b post-processing routines in the application.
             */
            inline bool get_rev_heat() const
            {
                return rev_heat;
            }

            /**
             * Method to get the fraction of reversible heat released corresponding to half-cell reaction of \p ORR, inside the cathode catalyst layer (#factor_rev_heat_ccl).
             * The typical use of this method is for the \p \b post-processing routines in the application.
             */
            inline double get_factor_rev_heat_ccl() const
            {
                return factor_rev_heat_ccl;
            }

            /**
             * Method to get whether the heat sink due to evaoration of water produced during the
             * ORR inside the athode catalyst layer is \p ON or \p OFF (#water_vap_heat_ccl).
             * The typical use of this method is for the \p \b post-processing routines in the application.
             */
            inline bool get_water_vap_heat_ccl() const
            {
                return water_vap_heat_ccl;
            }

            /**
             * Accessor for cathode catalyst layer FuelCellShop::Kinetics::BaseKinetics pointer.
             */
            inline FuelCellShop::Kinetics::BaseKinetics* get_cathode_kinetics() const
            {
                return cathode_kinetics;
            }

            /**
             * Accessor for anode catalyst layer FuelCellShop::Kinetics::BaseKinetics pointer.
             */
            inline FuelCellShop::Kinetics::BaseKinetics* get_anode_kinetics() const
            {
                return anode_kinetics;
            }

            //@}

        protected:

            ///@name Local CG FEM based assemblers - make_ functions
            //@{

            /**
             * This function computes Local CG FEM based
             * assemblers - constant data (generic).
             *
             * \warning <b>For the DEVELOPERS:</b> \p block_index for VariableInfo are not filled here, but \p indices_exist flag is set to \b TRUE (as it
             * is needed at other places before \p cell_matrix or \p cell_residual assembly).
             * \p Developers need to be wary of this fact that the \p block_index are still not filled yet, and they need to <b>fill them before</b> doing \p cell_matrix assembly.
             */
            virtual void make_assemblers_generic_constant_data();

            /**
             * This function computes
             * <p> Local CG FEM based assemblers - constant data (cell) </p>
             * and allocates the memory for
             * shape functions, shape function gradients, and
             * \p JxW_cell in
             * <p> Local CG FEM based assemblers - variable data (cell) </p>.
             */
            virtual void make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info);

            /**
             * This function computes
             * <p> Local CG FEM based assemblers - variable data (cell) </p>.
             */
            virtual void make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                            FuelCellShop::Layer::BaseLayer<dim>* const layer);

            /**
             * This function is specifically created for assembly of cell matrices for the following equations, \em viz.
             * - <b>"Proton Transport Equation"</b>
             * - <b>"Electron Transport Equation"</b>
             * - <b>"Ficks Transport Equation - oxygen"</b>
             * - <b>"Ficks Transport Equation - water"</b>
             * - <b>"Thermal Transport Equation"</b>
             * - <b>"Liquid Water Saturation Transport Equation"</b>
             *
             * Cell matrix terms for the abovementioned equations follow a similar pattern, using different factors; hence this function avoids
             * using repeated lines of code for each of these equations in  #assemble_cell_matrix method. The argument list is defined as follows:
             * @param - first corresponds to \p FuelCell::ApplicationCore::MatrixVector which we need to fill.
             * @param - second corresponds to typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo object; used to access fe elements for variables other than test function.
             * @param - third corresponds to \p std::string for the name of the equation for which we are doing the assembly. It should match with names mentioned above.
             * @param - fourth corresponds to \p FEValuesBase (fe discretization) used for the test function \em i.e. equation for which assembly is being done.
             * @param - fifth corresponds to shape functions for the test function, filled already by #make_assemblers_cell_variable_data.
             * @param - fifth corresponds to the factor used to multiply with the current density in source term for a particular equation.
             *
             * \warning <b>For the DEVELOPERS:</b> This function should only be used after #make_assemblers_cell_variable_data has been called in #assemble_cell_matrix, otherwise it is bound to give wrong results.
             */
            virtual void assemble_matrix_for_equation(FuelCell::ApplicationCore::MatrixVector& cell_matrices,
                                                      const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                      const std::string& eq_name,
                                                      const FEValuesBase<dim>& test_fe,
                                                      const std::vector< std::vector<double> >& test_shape_functions,
                                                      const double& sourceterm_factor);

            //@}

            ///@name Boolean flags / Parameters
            //@{

            /**
             * This boolean data member indicates that
             * the ir-reversible heating due to \p ORR, inside the cathode catalyst layer is \b ON or \b OFF.
             */
            bool irrev_heat_ccl;

            /**
             * This boolean data member indicates that
             * the ir-reversible heating due to \p HOR, inside the anode catalyst layer is \b ON or \b OFF.
             */
            bool irrev_heat_acl;

            /**
             * This boolean data member indicates that
             * the reversible (entropic) heat source due to net reaction forming liquid water product is \p ON or \p OFF.
             */
            bool rev_heat;

            /**
             * This represents the fraction of reversible heat released corresponding to half-cell reaction of \p ORR, inside the cathode catalyst layer.
             */
            double factor_rev_heat_ccl;

            /**
             * This boolean data member indicates that whether the water is produced in vapour phase (i.e. completely evaporated)
             * during \p ORR, inside the CCL.
             */
            bool water_vapour_phase;

            /**
             * This boolean data member indicates that
             * the heat sink due to evaporation of water produced during the ORR inside the cathode catalyst layer is \p ON or \p OFF.
             */
            bool water_vap_heat_ccl;

            //@}

            ///@name Kinetics Objects
            //@{

            /**
             * Pointer to Cathode Kinetics object, initialized in the constructor.
             */
            FuelCellShop::Kinetics::BaseKinetics* cathode_kinetics;

            /**
             * Pointer to Anode Kinetics object, initialized in the constructor.
             */
            FuelCellShop::Kinetics::BaseKinetics* anode_kinetics;

            //@}

            ///@name ReactionHeat Objects
            //@{

            /**
             * Pointer to Cathode ReactionHeat object. This is initialized automatically
             * inside this class, depending on whether \p temperature_of_REV is beig solved for or not,
             * and #cathode_kinetics exists or not.
             */
            FuelCellShop::Equation::ReactionHeat* cathode_reactionheat;

            /**
             * Pointer to Anode Reactionheat object. This is initialized automatically
             * inside this class, depending on whether \p temperature_of_REV is beig solved for or not,
             * and #anode_kinetics exists or not.
             */
            FuelCellShop::Equation::ReactionHeat* anode_reactionheat;

            //@}

            ///@name Generic Constant Data
            //@{

            /**
             * VariableInfo structure corresponding to \p "oxygen_molar_fraction".
             */
            VariableInfo x_oxygen;

            /**
             * VariableInfo structure corresponding to \p "water_molar_fraction".
             */
            VariableInfo x_water;

            /**
             * VariableInfo structure corresponding to \p "electronic_electrical_potential".
             */
            VariableInfo phi_s;

            /**
             * VariableInfo structure corresponding to \p "membrane_water_content".
             */
            VariableInfo lambda;


            /**
             * VariableInfo structure corresponding to \p "protonic_electrical_potential".
             */
            VariableInfo phi_m;

            /**
             * VariableInfo structure corresponding to \p "temperature_of_REV".
             */
            VariableInfo t_rev;

            /**
             * VariableInfo structure corresponding to \p "liquid_water_saturation".
             */
            VariableInfo s_liquid_water;

            /**
             * Universal Faraday's constant
             */
            double F;

            /**
             * Molar weight of water in grams/mole.
             */
            double M_water;

            //@}

            ///@name Local CG FEM based assemblers - variable data (cell)
            //@{

            /**
             * \f$ \mathbf{T} \f$ shape functions.
             *
             * \p phi_T_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{T} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_T_cell;

            /**
             * \f$ \mathbf{x_{O_2}} \f$ shape functions.
             *
             * \p phi_xOxygen_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{x_{O_2}} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_xOxygen_cell;

            /**
             * \f$ \mathbf{\phi_s} \f$ shape functions.
             *
             * \p phi_phiS_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\phi_s} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_phiS_cell;

            /**
             * \f$ \mathbf{\phi_m} \f$ shape functions.
             *
             * \p phi_phiM_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{\phi_m} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_phiM_cell;

            /**
             * \f$ \mathbf{x_{H_2O}} \f$ shape functions.
             *
             * \p phi_xWater_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{x_{H_2O}} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_xWater_cell;

            /**
             * \f$ \mathbf{s} \f$ shape functions.
             *
             * \p phi_s_cell \p[ \p q \p] \p[ \p k \p] denotes
             * \f$ k \f$-th \f$ \mathbf{s} \f$ shape function
             * computed in \f$ q \f$-th quadrature point of the cell.
             */
            std::vector< std::vector<double> > phi_s_cell;

            /**
             * Current produced [\p A/cm^3] at all quadrature points
             * in the cell.
             */
            std::vector<double> current_cell;

            /**
             * Heat produced due to electrochemical reaction [\p W/cm^2] at
             * all quadrature points in the cell.
             *
             * ReactionHeat object returns heat produced in [\p W/cm^2]. It needs to be multiplied by active_area_cell [\p cm^2/cm^3]
             * to convert into a source term [\p W/cm^3].
             */
            std::vector<double> heat_cell;

            /**
             * Derivative of current w.r.t. \f$ \mathbf{\phi_s} \f$
             * at all quadrature points in the cell.
             */
            std::vector<double> dcurrent_dphiS_cell;

            /**
             * Derivative of current w.r.t. \f$ \mathbf{\phi_m} \f$
             * at all quadrature points in the cell.
             */
            std::vector<double> dcurrent_dphiM_cell;

            /** Derivative of current w.r.t. \f$ \mathbf{T} \f$
             * at all quadrature points in the cell.
             */
            std::vector<double> dcurrent_dT_cell;

            /**
             * Derivative of current w.r.t. \f$ \mathbf{x_{O_2}} \f$
             * at all quadrature points in the cell.
             */
            std::vector<double> dcurrent_dxOxygen_cell;

            /**
             * Derivative of current w.r.t. \f$ \mathbf{x_{H_2O}} \f$
             * at all quadrature points in the cell.
             */
            std::vector<double> dcurrent_dxWater_cell;

            /**
             * Derivative of heat produced w.r.t. \f$ \mathbf{\phi_s} \f$
             * at all quadrature points in the cell.
             */
            std::vector<double> dheat_dphiS_cell;

            /**
             * Derivative of heat produced w.r.t. \f$ \mathbf{\phi_m} \f$
             * at all quadrature points in the cell.
             */
            std::vector<double> dheat_dphiM_cell;

            /**
             * Derivative of heat produced w.r.t. \f$ \mathbf{T} \f$
             * at all quadrature points in the cell.
             */
            std::vector<double> dheat_dT_cell;

            /**
             * Derivative of heat produced w.r.t. \f$ \mathbf{x_{O_2}} \f$
             * at all quadrature points in the cell.
             */
            std::vector<double> dheat_dxOxygen_cell;

            /**
             * Derivative of heat produced w.r.t. \f$ \mathbf{x_{H_2O}} \f$
             * at all quadrature points in the cell.
             */
            std::vector<double> dheat_dxWater_cell;

            /**
             * Active area [\p cm^2/cm^3] of the cell.
             * This is specifically required for heat source terms, as ReactionHeat return
             * heat terms in [\p W/cm^2].
             */
            double active_area_cell;

            /**
             * Factor for <p>Proton Transport Equation</p> source term, depending on whether the current cell
             * is in the Anode catalyst layer or Cathode catalyst layer.
             */
            double factor_protontranseq_cell;

            /**
             * Factor for <p>Electron Transport Equation</p> source term, depending on whether the current cell
             * is in the Anode catalyst layer or Cathode catalyst layer.
             */
            double factor_electrontranseq_cell;

            /**
             * Factor for <p>Ficks Transport Equation - oxygen</p> source term, depending on whether the current cell
             * is in the Anode catalyst layer or Cathode catalyst layer.
             */
            double factor_oxygentranseq_cell;

            /**
             * Factor for <p>Ficks Transport Equation - water</p> source term, depending on whether the current cell
             * is in the Anode catalyst layer or Cathode catalyst layer, and it also depends on whether the water is getting
             * produced in the liquid phase or vapour phase.
             */
            double factor_watertranseq_cell;

            /**
             * Factor for <p>Liquid Water Saturation Transport Equation</p> source term, depending on whether the current cell
             * is in the Anode catalyst layer or Cathode catalyst layer, and it also depends on whether the water is getting
             * produced in the liquid phase or vapour phase.
             */
            double factor_saturationtranseq_cell;

            //@}

            ///@name Counters
            //@{

            /**
             * Counter set to \em TRUE when \p cell_matrix is being assembled.
             * This ensures that only \b derivatives of source terms are computed in the cell, in #make_assemblers_cell_variable_data.
             *
             * \note <b>For developers:</b> Other counter #cell_residual_counter should be set to \em FALSE, when this counter is set to \em TRUE.
             */
            bool cell_matrix_counter;

            /**
             * Counter set to \em TRUE when \p cell_residual is being assembled.
             * This ensures that only source terms, \em e.g. \p Current or \p Reaction \p Heat, are computed in the cell, in #make_assemblers_cell_variable_data.
             *
             * \note <b>For developers:</b> Other counter #cell_matrix_counter should be set to \em FALSE, when this counter is set to \em TRUE.
             */
            bool cell_residual_counter;

            /**
             * Variable used to store the index in cell_info->global_data of the previous Newton solution
             * The solution at the previous iteration is used to compute cell_matrix and cell_residual
             */
            unsigned int last_iter_cell;

            //@}

///@name Methods and Data related to fluid transport equations
//@{

  /////////////
  // METHODS //
  /////////////

  /**
   * This function computes
   * some of the computational constants
   * and
   * Local CG FEM based assemblers - constant data (generic).
   *
   * The template parameter INFO is either
   * typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo
   * or
   * typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo.
   */
  template<typename INFO>
  void make_assemblers_generic_constant_data(const INFO&                                InFo,
                                             FuelCellShop::Layer::BaseLayer<dim>* const layer);

  /**
   * As before but for fluid transport equations.
   */
  void make_assemblers_cell_constant_data2(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info);

  /**
   * As before but for fluid transport equations.
   */
  void make_assemblers_cell_variable_data2(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                           FuelCellShop::Layer::BaseLayer<dim>* const                               layer);

  /**
   * This function fills out the
   * max number of \p matrix_block_indices which need to be updated due to the source terms.
   * The real number might be less depending on in which CL we are.
   */
  virtual void make_matrix_block_indices();

  /**
   * This function fills out the
   * max number of \p residual_indices which need to be updated due to the source terms.
   * The real number might be less depending on in which CL we are.
   */
  virtual void make_residual_indices();

  ////////////
  // DATA 1 //
  ////////////

  /**
   * Number of species, \f$ N \f$.
   */
  unsigned int n_species;

  /**
   * Keeps oxygen set of fluid transport equations.
   */
  unsigned int indexO2;

  /**
   * Keeps hydrogen set of fluid transport equations.
   */
  unsigned int indexH2;

  /**
   * Keeps water vapor set of fluid transport equations.
   */
  unsigned int indexH2O;

  ////////////
  // DATA 2 //
  ////////////

  /**
   * Oxygen source multiplier,      \f$ \quad \left[ \frac{\text{g}}{A sec} \right] \f$.
   */
  double multiplierO2;

  /**
   * Hydrogen source multiplier,    \f$ \quad \left[ \frac{\text{g}}{A sec} \right] \f$.
   */
  double multiplierH2;

  /**
   * Water vapor source multiplier, \f$ \quad \left[ \frac{\text{g}}{A sec} \right] \f$.
   */
  double multiplierH2O;

  ////////////
  // DATA 3 //
  ////////////

  /**
   * Molar mass of pure gas, \f$ M_i \quad \left[ \frac{\text{g}}{\text{mol}} \right] \f$.
   */
  std::vector<double> molar_mass;

  ////////////
  // DATA 4 //
  ////////////

  /**
   * Density extractors.
   */
  std::vector< FEValuesExtractors::Scalar > density_extractors;

  /**
   * Electronic electrical potential extractor.
   */
  FEValuesExtractors::Scalar electronic_electrical_potential_extractor;

  /**
   * Protonic electrical potential extractor.
   */
  FEValuesExtractors::Scalar protonic_electrical_potential_extractor;

  ////////////
  // DATA 5 //
  ////////////

  /**
   * Constant temperature of species mixture
   * in the quadrature points of a cell,
   * \f$ T_{\text{mixture}}^{\text{const}} \quad \left[ \text{K} \right] \f$.
   */
  std::vector<double> T_mixture;

  ////////////
  // DATA 7 //
  ////////////

  /**
   * ORR current density
   * in the quadrature points of a cell
   * at a previous Newton iteration.
   */
  std::vector<double> ORR_current_density;

  /**
   * ORR current density derivative
   * with respect to oxygen concentration (gas, NOT gas-liquid)
   * in the quadrature points of a cell
   * at a previous Newton iteration.
   */
  std::vector<double> DORR_current_density_Doxygen_concentration;

  /**
   * ORR current density derivative
   * with respect to electronic electrical potential
   * in the quadrature points of a cell
   * at a previous Newton iteration.
   */
  std::vector<double> DORR_current_density_Delectronic_electrical_potential;

  /**
   * ORR current density derivative
   * with respect to protonic electrical potential
   * in the quadrature points of a cell
   * at a previous Newton iteration.
   */
  std::vector<double> DORR_current_density_Dprotonic_electrical_potential;

  ////////////
  // DATA 8 //
  ////////////

  /**
   * HOR current density
   * in the quadrature points of a cell
   * at a previous Newton iteration.
   */
  std::vector<double> HOR_current_density;

  /**
   * HOR current density derivative
   * with respect to hydrogen concentration (gas, NOT gas-liquid)
   * in the quadrature points of a cell
   * at a previous Newton iteration.
   */
  std::vector<double> DHOR_current_density_Dhydrogen_concentration;

  /**
   * HOR current density derivative
   * with respect to electronic electrical potential
   * in the quadrature points of a cell
   * at a previous Newton iteration.
   */
  std::vector<double> DHOR_current_density_Delectronic_electrical_potential;

  /**
   * HOR current density derivative
   * with respect to protonic electrical potential
   * in the quadrature points of a cell
   * at a previous Newton iteration.
   */
  std::vector<double> DHOR_current_density_Dprotonic_electrical_potential;

  ////////////
  // DATA 9 //
  ////////////

  /**
   * Density of each species
   * in the quadrature points of a cell
   * at a previous Newton iteration.
   */
  std::vector< std::vector<double> > density_old;

  /**
   * Electronic electrical potential
   * in the quadrature points of a cell
   * at a previous Newton iteration.
   */
  std::vector<double> electronic_electrical_potential_old;

  /**
   * Protonic electrical potential
   * in the quadrature points of a cell
   * at a previous Newton iteration.
   */
  std::vector<double> protonic_electrical_potential_old;

  /////////////
  // DATA 10 //
  /////////////

  /**
   * Density shape functions.
   *
   * \p phi_density \p[ \p s \p] \p[ \p q \p] \p[ \p k \p] denotes
   * \f$ k \f$-th density shape function
   * computed in \f$ q \f$-th quadrature point of a cell
   * for species \f$ s \f$.
   */
  std::vector< std::vector< std::vector<double> > > phi_density;

  /**
   * Electronic electrical potential shape functions.
   *
   * \p phi_electronic_electrical_potential \p[ \p q \p] \p[ \p k \p] denotes
   * \f$ k \f$-th electronic electrical potential shape function
   * computed in \f$ q \f$-th quadrature point of a cell.
   */
  std::vector< std::vector<double> > phi_electronic_electrical_potential;

  /**
   * Protonic electrical potential shape functions.
   *
   * \p phi_protonic_electrical_potential \p[ \p q \p] \p[ \p k \p] denotes
   * \f$ k \f$-th protonic electrical potential shape function
   * computed in \f$ q \f$-th quadrature point of a cell.
   */
  std::vector< std::vector<double> > phi_protonic_electrical_potential;

  /////////////
  // DATA 11 //
  /////////////

  /**
   * For internal use only.
   */
  std::string eq_generic_prefix;

  /**
   * For internal use only.
   */
  std::vector<std::string> eq_postfixes;

  /**
   * For internal use only.
   */
  std::vector<std::string> var_postfixes;

  /**
   * For internal use only.
   */
  std::string eq_name;

  /**
   * For internal use only.
   */
  std::string var_name;

//@}

};

} // Equation

} // FuelCellShop

#endif