// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: channel.h
// - Description: This class describes a channel
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id$
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_LAYER_CHANNEL_H_
#define _FCST_FUELCELLSHOP_LAYER_CHANNEL_H_

#include "base_layer.h"
#include "experimental_fluid.h"
#include "experimental_multi_fluid.h"
#include "GasMixture.h"

namespace FuelCellShop
{
namespace Layer
{

/**
 * This class describes
 * a channel and stores
 * pointers to
 *
 * - \p ExperimentalFluid object,
 * - \p ExperimentalMultiFluid object,
 * - \p GasMixture object.
 *
 * This class also stores
 * the roughness of the channel walls
 * and
 * effective electronic conductivity
 * of solid plates.
 *
 * The functionality of
 * this class can be extended
 * if needed.
 *
 * \author Valentin N. Zingan, 2012
 */

template<int dim>
class Channel : public BaseLayer<dim>
{
public:

///@name Constructors, destructor, and initialization
//@{

  /**
   * Constructor.
   */
  Channel(const std::string& name);

  /**
   * Constructor.
   */
  Channel(const std::string&                         name,
          FuelCellShop::Material::ExperimentalFluid& fluid);

  /**
   * Constructor.
   */
  Channel(const std::string&                              name,
          FuelCellShop::Material::ExperimentalMultiFluid& multi_fluid);

  /**
   * Constructor.
   */
  Channel(const std::string&                  name,
          FuelCellShop::Material::GasMixture& gas_mixture);

  /**
   * Destructor.
   */
  virtual ~Channel();

  /**
   * Initialize \p fluid.
   */
  void initialize(FuelCellShop::Material::ExperimentalFluid& rfluid)
  {
    fluid       = &rfluid;
    multi_fluid = nullptr;
    gas_mixture = nullptr;
  }

  /**
   * Initialize \p multi_fluid.
   */
  void initialize(FuelCellShop::Material::ExperimentalMultiFluid& rmulti_fluid)
  {
    fluid       = nullptr;
    multi_fluid = &rmulti_fluid;
    gas_mixture = nullptr;
  }

  /**
   * Initialize \p gas_mixture.
   */
  void initialize(FuelCellShop::Material::GasMixture& rgas_mixture)
  {
    fluid       = nullptr;
    multi_fluid = nullptr;
    gas_mixture = &rgas_mixture;
  }

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
   * \p fluid.
   */
  const FuelCellShop::Material::ExperimentalFluid* const get_fluid() const
  {
    return fluid;
  }

  /**
   * This function returns
   * \p multi_fluid.
   */
  const FuelCellShop::Material::ExperimentalMultiFluid* const get_multi_fluid() const
  {
    return multi_fluid;
  }

  /**
   * This function returns
   * \p gas_mixture.
   */
  const FuelCellShop::Material::GasMixture* const get_gas_mixture() const
  {
    return gas_mixture;
  }

  /**
   * This function returns
   * \p roughness.
   */
  const double& get_roughness() const
  {
    return roughness;
  }

  /**
   * This function returns
   * \p effective_electronic_conductivity.
   */
  const Tensor<2,dim>& get_effective_electronic_conductivity() const
  {
    return effective_electronic_conductivity;
  }

  /**
   * This function returns
   * \p typeid of this class.
   *
   * All classes derived
   * from this class
   * must share
   * the same
   * \p typeid
   * for dynamic
   * swapping.
   */
  const std::type_info& get_base_type() const
  {
    return typeid(Channel<dim>);
  }

  /**
   * This function prints out
   * the layer properties.
   */
  virtual void print_layer_properties() const;

//@}

protected:

//////////
// DATA //
//////////

///@name Layer properties
//@{

  /**
   * - incompressible,
   * - isothermal,
   * - single-phase,
   * - single-component
   *
   * fluid.
   */
  FuelCellShop::Material::ExperimentalFluid* fluid;

  /**
   * - compressible,
   * - isothermal,
   * - single-phase,
   * - multi-component
   *
   * fluid.
   */
  FuelCellShop::Material::ExperimentalMultiFluid* multi_fluid;

  /**
   * Gas mixture.
   */
  FuelCellShop::Material::GasMixture* gas_mixture;

  /**
   * Roughness, m.
   */
  double roughness;

  /**
   * Effective electronic conductivity, \f$ \sigma_s^{\text{eff}} \quad \left[ \frac{\text{S}}{\text{m}} \right] \f$.
   */
  Tensor<2,dim> effective_electronic_conductivity;

//@}

};

} // Layer

} // FuelCellShop

#endif