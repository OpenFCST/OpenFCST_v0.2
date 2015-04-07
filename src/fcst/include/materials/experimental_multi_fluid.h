// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: experimental_multi_fluid.h
// - Description: This class describes a multi fluid
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: experimental_multi_fluid.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_MATERIAL_EXPERIMENTAL_MULTI_FLUID_H_
#define _FCST_FUELCELLSHOP_MATERIAL_EXPERIMENTAL_MULTI_FLUID_H_

#include "base_material.h"

namespace FuelCellShop
{
namespace Material
{

/**
 * This class describes
 *
 * - compressible,
 * - isothermal,
 * - single-phase,
 * - multi-component
 *
 * fluid.
 *
 * The functionality of
 * this class can be extended
 * if needed.
 *
 * \author Valentin N. Zingan, 2013
 */

class ExperimentalMultiFluid : public BaseMaterial
{
public:

///@name Constructors, destructor, and initialization
//@{

  /**
   * Constructor.
   */
  ExperimentalMultiFluid(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~ExperimentalMultiFluid();

  /**
   * Declare parameters.
   */
  virtual void declare_parameters(ParameterHandler& param) const;

  /**
   * Initialize parameters.
   */
  virtual void initialize(ParameterHandler& param);

//@}

///@name Accessors and info
//@{

  /**
   * This function returns
   * \p n_species.
   */
  const unsigned int& get_number_of_species() const
  {
    return n_species;
  }

  /**
   * This function returns
   * \p T_mixture.
   */
  const double& get_temperature_of_mixture() const
  {
    return T_mixture;
  }

  /**
   * This function returns
   * \p molar_mass.
   */
  const std::vector<double>& get_molar_mass() const
  {
    return molar_mass;
  }

  /**
   * This function returns
   * \p dynamic_viscosity.
   */
  const std::vector<double>& get_dynamic_viscosity() const
  {
    return dynamic_viscosity;
  }

  /**
   * This function returns
   * \p bulk_viscosity.
   */
  const std::vector<double>& get_bulk_viscosity() const
  {
    return bulk_viscosity;
  }

  /**
   * This function returns
   * \p maxwell_stefan_isobaric_diffusion_coefficient.
   */
  const Table< 2, double >& get_maxwell_stefan_isobaric_diffusion_coefficient() const
  {
    return maxwell_stefan_isobaric_diffusion_coefficient;
  }

  /**
   * This function prints out
   * the material properties.
   */
  virtual void print_material_properties() const;

//@}

protected:

  //////////
  // DATA //
  //////////

///@name Fluid properties
//@{

  /**
   * Number of species.
   */
  unsigned int n_species;

  /**
   * Temperature of mixture, K.
   */
  double T_mixture;

  /**
   * Molar mass, kg/mol.
   */
  std::vector<double> molar_mass;

  /**
   * Dynamic viscosity, Pa sec.
   */
  std::vector<double> dynamic_viscosity;

  /**
   * Bulk viscosity, Pa sec.
   */
  std::vector<double> bulk_viscosity;

  /**
   * Each entry of this structure defines
   * a Maxwell-Stefan isobaric diffusion coefficient of gas \f$ i \f$ in gas \f$ j \f$ \f$ \quad \f$ \f$ \displaystyle \sum_{l=1}^N p_l \cdot \mathscr{D}_{ij} \f$, Pa m^2/sec.
   */
  Table< 2, double > maxwell_stefan_isobaric_diffusion_coefficient;

//@}

};

} // Material

} // FuelCellShop

#endif