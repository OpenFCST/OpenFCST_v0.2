// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: GasMixture.h
// - Description: This class describes properties of gas mixtures
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: GasMixture.h 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#ifndef _FCST_FUELCELLSHOP_MATERIAL_GASMIXTURE_H_
#define _FCST_FUELCELLSHOP_MATERIAL_GASMIXTURE_H_

#define _DUMMY_ 1.e300

#include "PureGas.h"

namespace FuelCellShop
{
namespace Material
{

/**
 * This class describes properties of gas mixtures.
 *
 * This class contains the following data:
 *
 * - \p std::vector of pointers to the \p PureGas objects called \p gases,
 * - total pressure of the whole gas mixture called \p total_pressure,
 * - temperature of the whole gas mixture called \p temperature.
 *
 * The whole gas mixture is supposed to be isobaric
 * if some concrete value of \p total_pressure
 * is assigned by either using the
 * \p set_total_pressure()
 * function
 * or
 * defining it in the parameters file.
 *
 * If nothing happens, then the value of
 * \p total_pressure gets the \p _DUMMY_ number equal
 * \p 1.e300 and the whole gas mixture
 * is treated as non-isobaric.
 * In this case, the total pressure of
 * the whole gas mixture is
 * one of the solution variables.
 *
 * The whole gas mixture is supposed to be isothermal
 * if some concrete value of \p temperature
 * is assigned by either using the
 * \p set_temperature()
 * function
 * or
 * defining it in the parameters file.
 *
 * If nothing happens, then the value of
 * \p temperature gets the \p _DUMMY_ number equal
 * \p 1.e300 and the whole gas mixture
 * is treated as non-isothermal.
 * In this case, the temperature of
 * the whole gas mixture is
 * one of the solution variables.
 *
 * The following methods are used to compute
 * the gas mixture properties:
 *
 * - Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
 *   written in the Chapman Enskog form (binary gas mixture only), \f$ \quad \left[ \frac{\text{Pa } \text{m}^2}{\text{sec}} \right] \quad \f$
 *   \f$ D_{12} \equiv p \mathscr{D}_{12} = 1.8829 \cdot 10^{-2} \frac{\sqrt{ \left( \frac{1}{M_1} + \frac{1}{M_2} \right) T^3}}{\sigma_{12}^2 \Omega_{\mathscr{D}, 12}} \quad \f$
 *   and partial derivative \f$ \quad \frac{\partial D_{12}}{\partial T} \quad \f$
 *   where the binary collision integral \f$ \quad \Omega_{\mathscr{D}, 12} \quad \f$ is given by
 *   \f$ \quad \Omega_{\mathscr{D}, 12} = \frac{A_{\text{diff}}}{\left(\frac{kT}{\epsilon_{12}}\right)^{B_{\text{diff}}}} +
 *                                        \frac{C_{\text{diff}}}{e^{D_{\text{diff}}\frac{kT}{\epsilon_{12}}}} +
 *                                        \frac{E_{\text{diff}}}{e^{F_{\text{diff}}\frac{kT}{\epsilon_{12}}}} +
 *                                        \frac{G_{\text{diff}}}{e^{H_{\text{diff}}\frac{kT}{\epsilon_{12}}}} \f$
 *
 * - Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
 *   written in the Chapman Enskog form (binary gas mixture only), \f$ \quad \left[ \frac{\text{m}^2}{\text{sec}} \right] \quad \f$
 *   \f$ \mathscr{D}_{12} = 1.8583 \cdot 10^{-7} \frac{\sqrt{ \left( \frac{1}{M_1} + \frac{1}{M_2} \right) T^3}}{p_{\text{total}} \sigma_{12}^2 \Omega_{\mathscr{D}, 12}} \quad \f$
 *   and partial derivatives \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$ and \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$
 *   where the binary collision integral \f$ \quad \Omega_{\mathscr{D}, 12} \quad \f$ is given above
 *
 * - The table \f$ N_{\text{gases}} \times N_{\text{gases}} \f$ containing
 *   Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
 *   written in the Chapman Enskog form (ternary and more complicated gas mixtures), \f$ \quad \left[ \frac{\text{Pa } \text{m}^2}{\text{sec}} \right] \quad \f$
 *   \f$ D_{ij} \equiv p \mathscr{D}_{ij} = 1.8829 \cdot 10^{-2} \frac{\sqrt{ \left( \frac{1}{M_i} + \frac{1}{M_j} \right) T^3}}{\sigma_{ij}^2 \Omega_{\mathscr{D}, ij}} \quad \f$
 *   and partial derivatives \f$ \quad \frac{\partial D_{ij}}{\partial T} \quad \f$
 *   where the binary collision integral \f$ \quad \Omega_{\mathscr{D}, ij} \quad \f$ is given by
 *   \f$ \quad \Omega_{\mathscr{D}, ij} = \frac{A_{\text{diff}}}{\left(\frac{kT}{\epsilon_{ij}}\right)^{B_{\text{diff}}}} +
 *                                        \frac{C_{\text{diff}}}{e^{D_{\text{diff}}\frac{kT}{\epsilon_{ij}}}} +
 *                                        \frac{E_{\text{diff}}}{e^{F_{\text{diff}}\frac{kT}{\epsilon_{ij}}}} +
 *                                        \frac{G_{\text{diff}}}{e^{H_{\text{diff}}\frac{kT}{\epsilon_{ij}}}} \f$
 *
 * - The table \f$ N_{\text{gases}} \times N_{\text{gases}} \f$ containing
 *   Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
 *   written in the Chapman Enskog form (ternary and more complicated gas mixtures), \f$ \quad \left[ \frac{\text{m}^2}{\text{sec}} \right] \quad \f$
 *   \f$ \mathscr{D}_{ij} = 1.8583 \cdot 10^{-7} \frac{\sqrt{ \left( \frac{1}{M_i} + \frac{1}{M_j} \right) T^3}}{p_{\text{total}} \sigma_{ij}^2 \Omega_{\mathscr{D}, ij}} \quad \f$
 *   and partial derivatives \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$ and \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$
 *   where the binary collision integral \f$ \quad \Omega_{\mathscr{D}, ij} \quad \f$ is given above
 *
 * All those methods receive data in SI units.
 * All those methods return the results in SI units.
 *
 * For developers: please update this info appropriately if you add a new method.
 *
 * \author Valentin N. Zingan, 2013
 */

class GasMixture : public BaseMaterial
{
public:

///@name Constructors, destructor, and initialization
//@{

  /**
   * Constructor.
   */
  GasMixture(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~GasMixture();

  /**
   * Declare parameters.
   */
  virtual void declare_parameters(ParameterHandler& param) const;

  /**
   * Initialize parameters.
   */
  virtual void initialize(ParameterHandler& param);

  /**
   * This function sets
   * \p gases.
   */
  void set_gases(const std::vector< PureGas* >& rgases)
  {
         gases = rgases;
  }

  /**
   * This function sets
   * \p total_pressure.
   */
  void set_total_pressure(const double& rtotal_pressure)
  {
         total_pressure = rtotal_pressure;
  }

  /**
   * This function sets
   * \p temperature.
   */
  void set_temperature(const double& rtemperature)
  {
         temperature = rtemperature;
  }

//@}

///@name Accessors and info
//@{

  /**
   * This function returns
   * \p gases.
   */
  const std::vector< PureGas* >& get_gases() const
  {
         return gases;
  }

  /**
   * This function returns
   * \p total_pressure.
   */
  const double& get_total_pressure() const
  {
         return total_pressure;
  }

  /**
   * This function returns
   * \p temperature.
   */
  const double& get_temperature() const
  {
         return temperature;
  }

  /**
   * This function prints out
   * the material properties.
   */
  virtual void print_material_properties() const;

//@}

///@name Service functions. Chapman Enskog isobaric diffusion coefficient. Binary gas mixture only.
//@{

  /**
   * This function returns
   * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a constant temperature.
   */
  const double get_ChapmanEnskog_isobaric_diffusion_coefficient() const;

  /**
   * This function returns
   * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a constant temperature
   * in the quadrature points of a mesh entity.
   *
   * @param diffusion_coefficient - Chapman Enskog isobaric diffusion coefficient
   *                                at a constant temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_isobaric_diffusion_coefficient(std::vector<double>& diffusion_coefficient) const;

  /**
   * This function returns
   * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable temperature.
   *
   * @param temperature - temperature.
   */
  const double get_ChapmanEnskog_isobaric_diffusion_coefficient(const double& temperature) const;

  /**
   * This function returns
   * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable temperature
   * in the quadrature points of a mesh entity.
   *
   * @param temperature           - temperature in the quadrature points of a mesh entity,
   * @param diffusion_coefficient - Chapman Enskog isobaric diffusion coefficient
   *                                at a variable temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_isobaric_diffusion_coefficient(const std::vector<double>& temperature,
                                                        std::vector<double>&       diffusion_coefficient) const;

//@}

///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficient. Binary gas mixture only.
//@{

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial D_{12}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable temperature.
   *
   * @param temperature - temperature.
   */
  const double get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(const double& temperature) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial D_{12}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan isobaric diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable temperature
   * in the quadrature points of a mesh entity.
   *
   * @param temperature - temperature in the quadrature points of a mesh entity,
   * @param dst         - \f$ \frac{\partial D_{12}}{\partial T} \quad \f$
   *                      at a variable temperature in the quadrature points of a mesh entity.
   */
  void get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(const std::vector<double>& temperature,
                                                                      std::vector<double>&       dst) const;

//@}

///@name Service functions. Chapman Enskog diffusion coefficient. Binary gas mixture only.
//@{

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and temperature.
   */
  const double get_ChapmanEnskog_diffusion_coefficient() const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and temperature
   * in the quadrature points of a mesh entity.
   *
   * @param diffusion_coefficient - Chapman Enskog diffusion coefficient
   *                                at a constant total pressure and temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_diffusion_coefficient(std::vector<double>& diffusion_coefficient) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and variable temperature.
   *
   * @param temperature - temperature.
   */
  const double get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(const double& temperature) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and variable temperature
   * in the quadrature points of a mesh entity.
   *
   * @param temperature           - temperature in the quadrature points of a mesh entity,
   * @param diffusion_coefficient - Chapman Enskog diffusion coefficient
   *                                at a constant total pressure and variable temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(const std::vector<double>& temperature,
                                                                    std::vector<double>&       diffusion_coefficient) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and constant temperature.
   *
   * @param total_pressure - total pressure.
   */
  const double get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature(const double& total_pressure) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and constant temperature
   * in the quadrature points of a mesh entity.
   *
   * @param total_pressure        - total pressure in the quadrature points of a mesh entity,
   * @param diffusion_coefficient - Chapman Enskog diffusion coefficient
   *                                at a variable total pressure and constant temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature(const std::vector<double>& total_pressure,
                                                                       std::vector<double>&       diffusion_coefficient) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature.
   *
   * @param total_pressure - total pressure,
   * @param temperature    - temperature.
   */
  const double get_ChapmanEnskog_diffusion_coefficient(const double& total_pressure,
                                                       const double& temperature) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature
   * in the quadrature points of a mesh entity.
   *
   * @param total_pressure        - total pressure in the quadrature points of a mesh entity,
   * @param temperature           - temperature in the quadrature points of a mesh entity,
   * @param diffusion_coefficient - Chapman Enskog diffusion coefficient
   *                                at a variable total pressure and temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_diffusion_coefficient(const std::vector<double>& total_pressure,
                                               const std::vector<double>& temperature,
                                               std::vector<double>&       diffusion_coefficient) const;

//@}

///@name Service functions. Derivatives of Chapman Enskog diffusion coefficient. Binary gas mixture only.
//@{

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and constant temperature.
   *
   * @param total_pressure - total pressure.
   */
  const double get_DChapmanEnskog_diffusion_coefficient_Dpressure(const double& total_pressure) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and constant temperature
   * in the quadrature points of a mesh entity.
   *
   * @param total_pressure - total pressure in the quadrature points of a mesh entity,
   * @param dst            - \f$ \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$
   *                         at a variable total pressure and constant temperature in the quadrature points of a mesh entity.
   */
  void get_DChapmanEnskog_diffusion_coefficient_Dpressure(const std::vector<double>& total_pressure,
                                                          std::vector<double>&       dst) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature.
   *
   * @param total_pressure - total pressure,
   * @param temperature    - temperature.
   */
  const double get_DChapmanEnskog_diffusion_coefficient_Dpressure(const double& total_pressure,
                                                                  const double& temperature) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature
   * in the quadrature points of a mesh entity.
   *
   * @param total_pressure - total pressure in the quadrature points of a mesh entity,
   * @param temperature    - temperature in the quadrature points of a mesh entity,
   * @param dst            - \f$ \frac{\partial \mathscr{D}_{12}}{\partial p} \quad \f$
   *                         at a variable total pressure and temperature in the quadrature points of a mesh entity.
   */
  void get_DChapmanEnskog_diffusion_coefficient_Dpressure(const std::vector<double>& total_pressure,
                                                          const std::vector<double>& temperature,
                                                          std::vector<double>&       dst) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and variable temperature.
   *
   * @param temperature - temperature.
   */
  const double get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const double& temperature) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a constant total pressure and variable temperature
   * in the quadrature points of a mesh entity.
   *
   * @param temperature - temperature in the quadrature points of a mesh entity,
   * @param dst         - \f$ \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$
   *                      at a constant total pressure and variable temperature in the quadrature points of a mesh entity.
   */
  void get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const std::vector<double>& temperature,
                                                             std::vector<double>&       dst) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature.
   *
   * @param total_pressure - total pressure,
   * @param temperature    - temperature.
   */
  const double get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const double& total_pressure,
                                                                     const double& temperature) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficient of gas \f$ 1 \f$ in gas \f$ 2 \f$ (or vice-versa)
   * written in the Chapman Enskog form (binary gas mixture only) at a variable total pressure and temperature
   * in the quadrature points of a mesh entity.
   *
   * @param total_pressure - total pressure in the quadrature points of a mesh entity,
   * @param temperature    - temperature in the quadrature points of a mesh entity,
   * @param dst            - \f$ \frac{\partial \mathscr{D}_{12}}{\partial T} \quad \f$
   *                         at a variable total pressure and temperature in the quadrature points of a mesh entity.
   */
  void get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const std::vector<double>& total_pressure,
                                                             const std::vector<double>& temperature,
                                                             std::vector<double>&       dst) const;

//@}

///@name Service functions. Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
//@{

  /**
   * This function returns
   * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant temperature.
   */
  const Table< 2, double > get_ChapmanEnskog_isobaric_diffusion_coefficients() const;

  /**
   * This function returns
   * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant temperature
   * in the quadrature points of a mesh entity.
   *
   * @param diffusion_coefficients - Chapman Enskog isobaric diffusion coefficients
   *                                 at a constant temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_isobaric_diffusion_coefficients(std::vector< Table< 2, double > >& diffusion_coefficients) const;

  /**
   * This function returns
   * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable temperature.
   *
   * @param temperature - temperature.
   */
  const Table< 2, double > get_ChapmanEnskog_isobaric_diffusion_coefficients(const double& temperature) const;

  /**
   * This function returns
   * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable temperature
   * in the quadrature points of a mesh entity.
   *
   * @param temperature            - temperature in the quadrature points of a mesh entity,
   * @param diffusion_coefficients - Chapman Enskog isobaric diffusion coefficients
   *                                 at a variable temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_isobaric_diffusion_coefficients(const std::vector<double>&         temperature,
                                                         std::vector< Table< 2, double > >& diffusion_coefficients) const;

//@}

///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
//@{

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial D_{ij}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable temperature.
   *
   * @param temperature - temperature.
   */
  const Table< 2, double > get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(const double& temperature) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial D_{ij}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan isobaric diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable temperature
   * in the quadrature points of a mesh entity.
   *
   * @param temperature - temperature in the quadrature points of a mesh entity,
   * @param dst         - \f$ \frac{\partial D_{ij}}{\partial T} \quad \f$
   *                      at a variable temperature in the quadrature points of a mesh entity.
   */
  void get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(const std::vector<double>&         temperature,
                                                                       std::vector< Table< 2, double > >& dst) const;

//@}

///@name Service functions. Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.
//@{

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and temperature.
   */
  const Table< 2, double > get_ChapmanEnskog_diffusion_coefficients() const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and temperature
   * in the quadrature points of a mesh entity.
   *
   * @param diffusion_coefficients - Chapman Enskog diffusion coefficients
   *                                 at a constant total pressure and temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_diffusion_coefficients(std::vector< Table< 2, double > >& diffusion_coefficients) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and variable temperature.
   *
   * @param temperature - temperature.
   */
  const Table< 2, double > get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure(const double& temperature) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and variable temperature
   * in the quadrature points of a mesh entity.
   *
   * @param temperature            - temperature in the quadrature points of a mesh entity,
   * @param diffusion_coefficients - Chapman Enskog diffusion coefficients
   *                                 at a constant total pressure and variable temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure(const std::vector<double>&         temperature,
                                                                     std::vector< Table< 2, double > >& diffusion_coefficients) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and constant temperature.
   *
   * @param total_pressure - total pressure.
   */
  const Table< 2, double > get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature(const double& total_pressure) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and constant temperature
   * in the quadrature points of a mesh entity.
   *
   * @param total_pressure         - total pressure in the quadrature points of a mesh entity,
   * @param diffusion_coefficients - Chapman Enskog diffusion coefficients
   *                                 at a variable total pressure and constant temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature(const std::vector<double>&         total_pressure,
                                                                        std::vector< Table< 2, double > >& diffusion_coefficients) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature.
   *
   * @param total_pressure - total pressure,
   * @param temperature    - temperature.
   */
  const Table< 2, double > get_ChapmanEnskog_diffusion_coefficients(const double& total_pressure,
                                                                    const double& temperature) const;

  /**
   * This function returns
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature
   * in the quadrature points of a mesh entity.
   *
   * @param total_pressure         - total pressure in the quadrature points of a mesh entity,
   * @param temperature            - temperature in the quadrature points of a mesh entity,
   * @param diffusion_coefficients - Chapman Enskog diffusion coefficients
   *                                 at a variable total pressure and temperature in the quadrature points of a mesh entity.
   */
  void get_ChapmanEnskog_diffusion_coefficients(const std::vector<double>&         total_pressure,
                                                const std::vector<double>&         temperature,
                                                std::vector< Table< 2, double > >& diffusion_coefficients) const;

//@}

///@name Service functions. Derivatives of Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.
//@{

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and constant temperature.
   *
   * @param total_pressure - total pressure.
   */
  const Table< 2, double > get_DChapmanEnskog_diffusion_coefficients_Dpressure(const double& total_pressure) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and constant temperature
   * in the quadrature points of a mesh entity.
   *
   * @param total_pressure - total pressure in the quadrature points of a mesh entity,
   * @param dst            - \f$ \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$
   *                         at a variable total pressure and constant temperature in the quadrature points of a mesh entity.
   */
  void get_DChapmanEnskog_diffusion_coefficients_Dpressure(const std::vector<double>&         total_pressure,
                                                           std::vector< Table< 2, double > >& dst) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature.
   *
   * @param total_pressure - total pressure,
   * @param temperature    - temperature.
   */
  const Table< 2, double > get_DChapmanEnskog_diffusion_coefficients_Dpressure(const double& total_pressure,
                                                                               const double& temperature) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature
   * in the quadrature points of a mesh entity.
   *
   * @param total_pressure - total pressure in the quadrature points of a mesh entity,
   * @param temperature    - temperature in the quadrature points of a mesh entity,
   * @param dst            - \f$ \frac{\partial \mathscr{D}_{ij}}{\partial p} \quad \f$
   *                         at a variable total pressure and temperature in the quadrature points of a mesh entity.
   */
  void get_DChapmanEnskog_diffusion_coefficients_Dpressure(const std::vector<double>&         total_pressure,
                                                           const std::vector<double>&         temperature,
                                                           std::vector< Table< 2, double > >& dst) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and variable temperature.
   *
   * @param temperature - temperature.
   */
  const Table< 2, double > get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const double& temperature) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a constant total pressure and variable temperature
   * in the quadrature points of a mesh entity.
   *
   * @param temperature - temperature in the quadrature points of a mesh entity,
   * @param dst         - \f$ \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$
   *                      at a constant total pressure and variable temperature in the quadrature points of a mesh entity.
   */
  void get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const std::vector<double>&         temperature,
                                                              std::vector< Table< 2, double > >& dst) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature.
   *
   * @param total_pressure - total pressure,
   * @param temperature    - temperature.
   */
  const Table< 2, double > get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const double& total_pressure,
                                                                                  const double& temperature) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$ of the
   * Maxwell-Stefan diffusion coefficients of gas \f$ i \f$ in gas \f$ j \f$
   * written in the Chapman Enskog form (ternary and more complicated gas mixtures) at a variable total pressure and temperature
   * in the quadrature points of a mesh entity.
   *
   * @param total_pressure - total pressure in the quadrature points of a mesh entity,
   * @param temperature    - temperature in the quadrature points of a mesh entity,
   * @param dst            - \f$ \frac{\partial \mathscr{D}_{ij}}{\partial T} \quad \f$
   *                         at a variable total pressure and temperature in the quadrature points of a mesh entity.
   */
  void get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const std::vector<double>&         total_pressure,
                                                              const std::vector<double>&         temperature,
                                                              std::vector< Table< 2, double > >& dst) const;

//@}

///@name Service functions. Binary collision integral.
//@{

  /**
   * This function returns
   * binary collision integral
   * at a constant temperature.
   *
   * @param N1 - number of the first  gas from \p gases,
   * @param N2 - number of the second gas from \p gases.
   */
  const double get_binary_collision_integral(const unsigned int& N1 = 0,
                                             const unsigned int& N2 = 1) const;

  /**
   * This function returns
   * binary collision integral
   * at a constant temperature
   * in the quadrature points of a mesh entity.
   *
   * @param binary_collision_integral - binary collision integral at a constant temperature
   *                                    in the quadrature points of a mesh entity,
   * @param N1                        - number of the first  gas from \p gases,
   * @param N2                        - number of the second gas from \p gases.
   */
  void get_binary_collision_integral(std::vector<double>& binary_collision_integral,
                                     const unsigned int&  N1 = 0,
                                     const unsigned int&  N2 = 1) const;

  /**
   * This function returns
   * binary collision integral
   * at a variable temperature.
   *
   * @param temperature - temperature,
   * @param N1          - number of the first  gas from \p gases,
   * @param N2          - number of the second gas from \p gases.
   */
  const double get_binary_collision_integral(const double&       temperature,
                                             const unsigned int& N1 = 0,
                                             const unsigned int& N2 = 1) const;

  /**
   * This function returns
   * binary collision integral
   * at a variable temperature
   * in the quadrature points of a mesh entity.
   *
   * @param temperature               - temperature in the quadrature points of a mesh entity,
   * @param binary_collision_integral - binary collision integral at a variable temperature
   *                                    in the quadrature points of a mesh entity,
   * @param N1                        - number of the first  gas from \p gases,
   * @param N2                        - number of the second gas from \p gases.
   */
  void get_binary_collision_integral(const std::vector<double>& temperature,
                                     std::vector<double>&       binary_collision_integral,
                                     const unsigned int&        N1 = 0,
                                     const unsigned int&        N2 = 1) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \Omega_{\mathscr{D}, ij}}{\partial T} \quad \f$ of the
   * binary collision integral
   * at a variable temperature.
   *
   * @param temperature - temperature,
   * @param N1          - number of the first  gas from \p gases,
   * @param N2          - number of the second gas from \p gases.
   */
  const double get_Dbinary_collision_integral_Dtemperature(const double&       temperature,
                                                           const unsigned int& N1 = 0,
                                                           const unsigned int& N2 = 1) const;

  /**
   * This function returns
   * the first derivative \f$ \quad \frac{\partial \Omega_{\mathscr{D}, ij}}{\partial T} \quad \f$ of the
   * binary collision integral
   * at a variable temperature
   * in the quadrature points of a mesh entity.
   *
   * @param temperature - temperature in the quadrature points of a mesh entity,
   * @param dst         - \f$ \frac{\partial \Omega_{\mathscr{D}, ij}}{\partial T} \quad \f$
   *                      at a variable temperature in the quadrature points of a mesh entity,
   * @param N1          - number of the first  gas from \p gases,
   * @param N2          - number of the second gas from \p gases.
   */
  void get_Dbinary_collision_integral_Dtemperature(const std::vector<double>& temperature,
                                                   std::vector<double>&       dst,
                                                   const unsigned int&        N1 = 0,
                                                   const unsigned int&        N2 = 1) const;

//@}

protected:

  //////////
  // DATA //
  //////////

///@name Fluid properties
//@{

  /**
   * This \p std::vector contains all pure gases
   * which form the whole gas mixture of a problem
   * at hand.
   */
  std::vector< PureGas* > gases;

  /**
   * Total pressure of the whole
   * gas mixture, \f$ p_{\text{total}} \quad \left[ \text{Pa} \right] \f$.
   *
   * If \p total_pressure = \p _DUMMY_ = \p 1.e300, then the whole
   * gas mixture is supposed to be a non-isobaric mixture.
   */
  double total_pressure;

  /**
   * Temperature of the whole
   * gas mixture, \f$ T \quad \left[ \text{K} \right] \f$.
   *
   * If \p temperature = \p _DUMMY_ = \p 1.e300, then the whole
   * gas mixture is supposed to be a non-isothermal mixture.
   */
  double temperature;

//@}

};

} // Material

} // FuelCellShop

#endif