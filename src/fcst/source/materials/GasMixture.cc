// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: GasMixture.cc
// - Description: This class describes properties of gas mixtures
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: GasMixture.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include "GasMixture.h"

namespace NAME = FuelCellShop::Material;

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

NAME::GasMixture::GasMixture(const std::string& name)
:
NAME::BaseMaterial(name)
{
       total_pressure = _DUMMY_;
       temperature    = _DUMMY_;
}

// ---            ---
// --- Destructor ---
// ---            ---

NAME::GasMixture::~GasMixture()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

void
NAME::GasMixture::declare_parameters(ParameterHandler& param) const
{
  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection("Materials");
    {
      param.enter_subsection(this->name);
      {
             param.declare_entry("Total pressure of gas mixture",
                                 "1.e300",
                                  Patterns::Double(0.0),
                                 " ");

             param.declare_entry("Temperature of gas mixture",
                                 "1.e300",
                                  Patterns::Double(0.0),
                                 " ");
      }
      param.leave_subsection();
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

void
NAME::GasMixture::initialize(ParameterHandler& param)
{
  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection("Materials");
    {
      param.enter_subsection(this->name);
      {
             total_pressure = param.get_double("Total pressure of gas mixture");
             temperature    = param.get_double("Temperature of gas mixture");
      }
      param.leave_subsection();
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

       ///////////////////////
       ///////////////////////
       // SERVICE FUNCTIONS //
       ///////////////////////
       ///////////////////////

  /////////////////////////////////////////////////////////////////////////////
  // Chapman Enskog isobaric diffusion coefficient - Binary gas mixture only //
  /////////////////////////////////////////////////////////////////////////////

// ---                                                  ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficient ---
// ---                                                  ---

const double
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficient() const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       const double M1 = gases[0]->get_molar_mass() * 1.0e3;
       const double M2 = gases[1]->get_molar_mass() * 1.0e3;

       const double sigma12 = 0.5*( gases[0]->get_collision_diameter() + gases[1]->get_collision_diameter() );

       const double T = temperature;

       const double omega12 = get_binary_collision_integral();

       return 1.8829e-2 * std::sqrt( ( 1.0/M1 + 1.0/M2 )*T*T*T ) / ( sigma12*sigma12*omega12 );
}

// ---                                                  ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficient ---
// ---                                                  ---

void
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficient(std::vector<double>& diffusion_coefficient) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       const double result = get_ChapmanEnskog_isobaric_diffusion_coefficient();

       for(unsigned int q = 0; q < diffusion_coefficient.size(); ++q)
              diffusion_coefficient[q] = result;
}

// ---                                                  ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficient ---
// ---                                                  ---

const double
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficient(const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       const double M1 = gases[0]->get_molar_mass() * 1.0e3;
       const double M2 = gases[1]->get_molar_mass() * 1.0e3;

       const double sigma12 = 0.5*( gases[0]->get_collision_diameter() + gases[1]->get_collision_diameter() );

       const double T = temp;

       const double omega12 = get_binary_collision_integral(temp);

       return 1.8829e-2 * std::sqrt( ( 1.0/M1 + 1.0/M2 )*T*T*T ) / ( sigma12*sigma12*omega12 );
}

// ---                                                  ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficient ---
// ---                                                  ---

void
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficient(const std::vector<double>& temp,
                                                                   std::vector<double>&       diffusion_coefficient) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( temp.size() == diffusion_coefficient.size() , ExcDimensionMismatch(temp.size(), diffusion_coefficient.size()) );

       const double M1 = gases[0]->get_molar_mass() * 1.0e3;
       const double M2 = gases[1]->get_molar_mass() * 1.0e3;

       const double sigma12 = 0.5*( gases[0]->get_collision_diameter() + gases[1]->get_collision_diameter() );

       std::vector<double> omega12;
       omega12.resize(temp.size());

       get_binary_collision_integral(temp,
                                     omega12);

       for(unsigned int q = 0; q < temp.size(); ++q)
              diffusion_coefficient[q] = 1.8829e-2 * std::sqrt( ( 1.0/M1 + 1.0/M2 )*temp[q]*temp[q]*temp[q] ) / ( sigma12*sigma12*omega12[q] );
}

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Derivatives of Chapman Enskog isobaric diffusion coefficient - Binary gas mixture only //
  ////////////////////////////////////////////////////////////////////////////////////////////

// ---                                                                ---
// --- get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature ---
// ---                                                                ---

const double
NAME::GasMixture::get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       const double D12         = get_ChapmanEnskog_isobaric_diffusion_coefficient(temp);
       const double T           = temp;
       const double omega12     = get_binary_collision_integral(temp);
       const double Domega12_DT = get_Dbinary_collision_integral_Dtemperature(temp);

       return D12*(1.5/T - Domega12_DT/omega12);
}

// ---                                                                ---
// --- get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature ---
// ---                                                                ---

void
NAME::GasMixture::get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(const std::vector<double>& temp,
                                                                                 std::vector<double>&       dst) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( temp.size() == dst.size() , ExcDimensionMismatch(temp.size(), dst.size()) );

       std::vector<double> D12;
       D12.resize(temp.size());

       get_ChapmanEnskog_isobaric_diffusion_coefficient(temp,
                                                        D12);

       std::vector<double> omega12;
       omega12.resize(temp.size());

       get_binary_collision_integral(temp,
                                     omega12);

       std::vector<double> Domega12_DT;
       Domega12_DT.resize(temp.size());

       get_Dbinary_collision_integral_Dtemperature(temp,
                                                   Domega12_DT);

       for(unsigned int q = 0; q < temp.size(); ++q)
              dst[q] = D12[q]*(1.5/temp[q] - Domega12_DT[q]/omega12[q]);
}

  ////////////////////////////////////////////////////////////////////
  // Chapman Enskog diffusion coefficient - Binary gas mixture only //
  ////////////////////////////////////////////////////////////////////

// ---                                         ---
// --- get_ChapmanEnskog_diffusion_coefficient ---
// ---                                         ---

const double
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient() const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       return get_ChapmanEnskog_isobaric_diffusion_coefficient() / total_pressure;
}

// ---                                         ---
// --- get_ChapmanEnskog_diffusion_coefficient ---
// ---                                         ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient(std::vector<double>& diffusion_coefficient) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       const double result = get_ChapmanEnskog_diffusion_coefficient();

       for(unsigned int q = 0; q < diffusion_coefficient.size(); ++q)
              diffusion_coefficient[q] = result;
}

// ---                                                              ---
// --- get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure ---
// ---                                                              ---

const double
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       return get_ChapmanEnskog_isobaric_diffusion_coefficient(temp) / total_pressure;
}

// ---                                                              ---
// --- get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure ---
// ---                                                              ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(const std::vector<double>& temp,
                                                                               std::vector<double>&       diffusion_coefficient) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       AssertThrow( temp.size() == diffusion_coefficient.size() , ExcDimensionMismatch(temp.size(), diffusion_coefficient.size()) );

       std::vector<double> D12;
       D12.resize(temp.size());

       get_ChapmanEnskog_isobaric_diffusion_coefficient(temp,
                                                        D12);

       for(unsigned int q = 0; q < temp.size(); ++q)
              diffusion_coefficient[q] = D12[q] / total_pressure;
}

// ---                                                                 ---
// --- get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature ---
// ---                                                                 ---

const double
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature(const double& total_pres) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       return get_ChapmanEnskog_isobaric_diffusion_coefficient() / total_pres;
}

// ---                                                                 ---
// --- get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature ---
// ---                                                                 ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient_at_constant_temperature(const std::vector<double>& total_pres,
                                                                                  std::vector<double>&       diffusion_coefficient) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       AssertThrow( total_pres.size() == diffusion_coefficient.size() , ExcDimensionMismatch(total_pres.size(), diffusion_coefficient.size()) );

       std::vector<double> D12;
       D12.resize(total_pres.size());

       get_ChapmanEnskog_isobaric_diffusion_coefficient(D12);

       for(unsigned int q = 0; q < total_pres.size(); ++q)
              diffusion_coefficient[q] = D12[q] / total_pres[q];
}

// ---                                         ---
// --- get_ChapmanEnskog_diffusion_coefficient ---
// ---                                         ---

const double
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient(const double& total_pres,
                                                          const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       return get_ChapmanEnskog_isobaric_diffusion_coefficient(temp) / total_pres;
}

// ---                                         ---
// --- get_ChapmanEnskog_diffusion_coefficient ---
// ---                                         ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficient(const std::vector<double>& total_pres,
                                                          const std::vector<double>& temp,
                                                          std::vector<double>&       diffusion_coefficient) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( total_pres.size() == diffusion_coefficient.size() , ExcDimensionMismatch(total_pres.size(), diffusion_coefficient.size()) );
       AssertThrow( temp.size()       == diffusion_coefficient.size() , ExcDimensionMismatch(temp.size(),       diffusion_coefficient.size()) );

       std::vector<double> D12;
       D12.resize(temp.size());

       get_ChapmanEnskog_isobaric_diffusion_coefficient(temp,
                                                        D12);

       for(unsigned int q = 0; q < temp.size(); ++q)
              diffusion_coefficient[q] = D12[q] / total_pres[q];
}

  ///////////////////////////////////////////////////////////////////////////////////
  // Derivatives of Chapman Enskog diffusion coefficient - Binary gas mixture only //
  ///////////////////////////////////////////////////////////////////////////////////

// ---                                                    ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dpressure ---
// ---                                                    ---

const double
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dpressure(const double& total_pres) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       return - get_ChapmanEnskog_isobaric_diffusion_coefficient() / (total_pres*total_pres);
}

// ---                                                    ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dpressure ---
// ---                                                    ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dpressure(const std::vector<double>& total_pres,
                                                                     std::vector<double>&       dst) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );

       std::vector<double> D12;
       D12.resize(total_pres.size());

       get_ChapmanEnskog_isobaric_diffusion_coefficient(D12);

       for(unsigned int q = 0; q < total_pres.size(); ++q)
              dst[q] = - D12[q] / (total_pres[q]*total_pres[q]);
}

// ---                                                    ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dpressure ---
// ---                                                    ---

const double
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dpressure(const double& total_pres,
                                                                     const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       return - get_ChapmanEnskog_isobaric_diffusion_coefficient(temp) / (total_pres*total_pres);
}

// ---                                                    ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dpressure ---
// ---                                                    ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dpressure(const std::vector<double>& total_pres,
                                                                     const std::vector<double>& temp,
                                                                     std::vector<double>&       dst) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );
       AssertThrow( temp.size()       == dst.size() , ExcDimensionMismatch(temp.size(),       dst.size()) );

       std::vector<double> D12;
       D12.resize(temp.size());

       get_ChapmanEnskog_isobaric_diffusion_coefficient(temp,
                                                        D12);

       for(unsigned int q = 0; q < temp.size(); ++q)
              dst[q] = - D12[q] / (total_pres[q]*total_pres[q]);
}

// ---                                                       ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dtemperature ---
// ---                                                       ---

const double
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       return get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(temp) / total_pressure;
}

// ---                                                       ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dtemperature ---
// ---                                                       ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const std::vector<double>& temp,
                                                                        std::vector<double>&       dst) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       AssertThrow( temp.size() == dst.size() , ExcDimensionMismatch(temp.size(), dst.size()) );

       std::vector<double> DD12_DT;
       DD12_DT.resize(temp.size());

       get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(temp,
                                                                      DD12_DT);

       for(unsigned int q = 0; q < temp.size(); ++q)
              dst[q] = DD12_DT[q] / total_pressure;
}

// ---                                                       ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dtemperature ---
// ---                                                       ---

const double
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const double& total_pres,
                                                                        const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       return get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(temp) / total_pres;
}

// ---                                                       ---
// --- get_DChapmanEnskog_diffusion_coefficient_Dtemperature ---
// ---                                                       ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficient_Dtemperature(const std::vector<double>& total_pres,
                                                                        const std::vector<double>& temp,
                                                                        std::vector<double>&       dst) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( gases.size() == 2,
                    ExcMessage("The mixture is NOT binary. Use another function.") );

       AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );
       AssertThrow( temp.size()       == dst.size() , ExcDimensionMismatch(temp.size(),       dst.size()) );

       std::vector<double> DD12_DT;
       DD12_DT.resize(temp.size());

       get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(temp,
                                                                      DD12_DT);

       for(unsigned int q = 0; q < temp.size(); ++q)
              dst[q] = DD12_DT[q] / total_pres[q];
}

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Chapman Enskog isobaric diffusion coefficients - Ternary and more complicated gas mixtures //
  ////////////////////////////////////////////////////////////////////////////////////////////////

// ---                                                   ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficients ---
// ---                                                   ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficients() const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       Table< 2, double > result(gases.size(), gases.size());

       for(unsigned int i = 0; i < gases.size(); ++i)
              for(unsigned int j = 0; j < gases.size(); ++j)
                     if( i != j )
                     {
                            const double Mi = gases[i]->get_molar_mass() * 1.0e3;
                            const double Mj = gases[j]->get_molar_mass() * 1.0e3;

                            const double sigmaij = 0.5*( gases[i]->get_collision_diameter() + gases[j]->get_collision_diameter() );

                            const double T = temperature;

                            const double omegaij = get_binary_collision_integral(i,j);

                            result(i,j) = 1.8829e-2 * std::sqrt( ( 1.0/Mi + 1.0/Mj )*T*T*T ) / ( sigmaij*sigmaij*omegaij );
                     }

       return result;
}

// ---                                                   ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficients ---
// ---                                                   ---

void
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficients(std::vector< Table< 2, double > >& diffusion_coefficients) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       Table< 2, double > result(gases.size(), gases.size());
       result = get_ChapmanEnskog_isobaric_diffusion_coefficients();

       for(unsigned int q = 0; q < diffusion_coefficients.size(); ++q)
       {
              diffusion_coefficients[q].reinit(gases.size(), gases.size());
              diffusion_coefficients[q] = result;
       }
}

// ---                                                   ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficients ---
// ---                                                   ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficients(const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       Table< 2, double > result(gases.size(), gases.size());

       for(unsigned int i = 0; i < gases.size(); ++i)
              for(unsigned int j = 0; j < gases.size(); ++j)
                     if( i != j )
                     {
                            const double Mi = gases[i]->get_molar_mass() * 1.0e3;
                            const double Mj = gases[j]->get_molar_mass() * 1.0e3;

                            const double sigmaij = 0.5*( gases[i]->get_collision_diameter() + gases[j]->get_collision_diameter() );

                            const double T = temp;

                            const double omegaij = get_binary_collision_integral(temp,
                                                                                 i,j);

                            result(i,j) = 1.8829e-2 * std::sqrt( ( 1.0/Mi + 1.0/Mj )*T*T*T ) / ( sigmaij*sigmaij*omegaij );
                     }

       return result;
}

// ---                                                   ---
// --- get_ChapmanEnskog_isobaric_diffusion_coefficients ---
// ---                                                   ---

void
NAME::GasMixture::get_ChapmanEnskog_isobaric_diffusion_coefficients(const std::vector<double>&         temp,
                                                                    std::vector< Table< 2, double > >& diffusion_coefficients) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( temp.size() == diffusion_coefficients.size() , ExcDimensionMismatch(temp.size(), diffusion_coefficients.size()) );

       for(unsigned int q = 0; q < temp.size(); ++q)
       {
              diffusion_coefficients[q].reinit(gases.size(), gases.size());
              diffusion_coefficients[q] = get_ChapmanEnskog_isobaric_diffusion_coefficients(temp[q]);
       }
}

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Derivatives of Chapman Enskog isobaric diffusion coefficients - Ternary and more complicated gas mixtures //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ---                                                                 ---
// --- get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature ---
// ---                                                                 ---

const Table< 2, double >
NAME::GasMixture::get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       Table< 2, double > D(gases.size(), gases.size());
       D = get_ChapmanEnskog_isobaric_diffusion_coefficients(temp);

       Table< 2, double > result(gases.size(), gases.size());

       for(unsigned int i = 0; i < gases.size(); ++i)
              for(unsigned int j = 0; j < gases.size(); ++j)
                     if( i != j )
                     {
                            const double T = temp;

                            const double omegaij = get_binary_collision_integral(temp,
                                                                                 i,j);

                            const double Domegaij_DT = get_Dbinary_collision_integral_Dtemperature(temp,
                                                                                                   i,j);

                            result(i,j) = D(i,j)*(1.5/T - Domegaij_DT/omegaij);
                     }

       return result;
}

// ---                                                                 ---
// --- get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature ---
// ---                                                                 ---

void
NAME::GasMixture::get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(const std::vector<double>&         temp,
                                                                                  std::vector< Table< 2, double > >& dst) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( temp.size() == dst.size() , ExcDimensionMismatch(temp.size(), dst.size()) );

       for(unsigned int q = 0; q < temp.size(); ++q)
       {
              dst[q].reinit(gases.size(), gases.size());
              dst[q] = get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(temp[q]);
       }
}

  ///////////////////////////////////////////////////////////////////////////////////////
  // Chapman Enskog diffusion coefficients - Ternary and more complicated gas mixtures //
  ///////////////////////////////////////////////////////////////////////////////////////

// ---                                          ---
// --- get_ChapmanEnskog_diffusion_coefficients ---
// ---                                          ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients() const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       Table< 2, double > D(gases.size(), gases.size());
       D = get_ChapmanEnskog_isobaric_diffusion_coefficients();

       Table< 2, double > result(gases.size(), gases.size());

       for(unsigned int i = 0; i < gases.size(); ++i)
              for(unsigned int j = 0; j < gases.size(); ++j)
                     if( i != j )
                            result(i,j) = D(i,j) / total_pressure;

       return result;
}

// ---                                          ---
// --- get_ChapmanEnskog_diffusion_coefficients ---
// ---                                          ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients(std::vector< Table< 2, double > >& diffusion_coefficients) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       Table< 2, double > result(gases.size(), gases.size());
       result = get_ChapmanEnskog_diffusion_coefficients();

       for(unsigned int q = 0; q < diffusion_coefficients.size(); ++q)
       {
              diffusion_coefficients[q].reinit(gases.size(), gases.size());
              diffusion_coefficients[q] = result;
       }
}

// ---                                                               ---
// --- get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure ---
// ---                                                               ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure(const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       Table< 2, double > D(gases.size(), gases.size());
       D = get_ChapmanEnskog_isobaric_diffusion_coefficients(temp);

       Table< 2, double > result(gases.size(), gases.size());

       for(unsigned int i = 0; i < gases.size(); ++i)
              for(unsigned int j = 0; j < gases.size(); ++j)
                     if( i != j )
                            result(i,j) = D(i,j) / total_pressure;

       return result;
}

// ---                                                               ---
// --- get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure ---
// ---                                                               ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure(const std::vector<double>&         temp,
                                                                                std::vector< Table< 2, double > >& diffusion_coefficients) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       AssertThrow( temp.size() == diffusion_coefficients.size() , ExcDimensionMismatch(temp.size(), diffusion_coefficients.size()) );

       for(unsigned int q = 0; q < temp.size(); ++q)
       {
              diffusion_coefficients[q].reinit(gases.size(), gases.size());
              diffusion_coefficients[q] = get_ChapmanEnskog_diffusion_coefficients_at_constant_pressure(temp[q]);
       }
}

// ---                                                                  ---
// --- get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature ---
// ---                                                                  ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature(const double& total_pres) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       Table< 2, double > D(gases.size(), gases.size());
       D = get_ChapmanEnskog_isobaric_diffusion_coefficients();

       Table< 2, double > result(gases.size(), gases.size());

       for(unsigned int i = 0; i < gases.size(); ++i)
              for(unsigned int j = 0; j < gases.size(); ++j)
                     if( i != j )
                            result(i,j) = D(i,j) / total_pres;

       return result;
}

// ---                                                                  ---
// --- get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature ---
// ---                                                                  ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature(const std::vector<double>&         total_pres,
                                                                                   std::vector< Table< 2, double > >& diffusion_coefficients) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       AssertThrow( total_pres.size() == diffusion_coefficients.size() , ExcDimensionMismatch(total_pres.size(), diffusion_coefficients.size()) );

       for(unsigned int q = 0; q < total_pres.size(); ++q)
       {
              diffusion_coefficients[q].reinit(gases.size(), gases.size());
              diffusion_coefficients[q] = get_ChapmanEnskog_diffusion_coefficients_at_constant_temperature(total_pres[q]);
       }
}

// ---                                          ---
// --- get_ChapmanEnskog_diffusion_coefficients ---
// ---                                          ---

const Table< 2, double >
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients(const double& total_pres,
                                                           const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       Table< 2, double > D(gases.size(), gases.size());
       D = get_ChapmanEnskog_isobaric_diffusion_coefficients(temp);

       Table< 2, double > result(gases.size(), gases.size());

       for(unsigned int i = 0; i < gases.size(); ++i)
              for(unsigned int j = 0; j < gases.size(); ++j)
                     if( i != j )
                            result(i,j) = D(i,j) / total_pres;

       return result;
}

// ---                                          ---
// --- get_ChapmanEnskog_diffusion_coefficients ---
// ---                                          ---

void
NAME::GasMixture::get_ChapmanEnskog_diffusion_coefficients(const std::vector<double>&         total_pres,
                                                           const std::vector<double>&         temp,
                                                           std::vector< Table< 2, double > >& diffusion_coefficients) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( total_pres.size() == diffusion_coefficients.size() , ExcDimensionMismatch(total_pres.size(), diffusion_coefficients.size()) );
       AssertThrow( temp.size()       == diffusion_coefficients.size() , ExcDimensionMismatch(temp.size(),       diffusion_coefficients.size()) );

       for(unsigned int q = 0; q < temp.size(); ++q)
       {
              diffusion_coefficients[q].reinit(gases.size(), gases.size());
              diffusion_coefficients[q] = get_ChapmanEnskog_diffusion_coefficients(total_pres[q],
                                                                                   temp[q]);
       }
}

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Derivatives of Chapman Enskog diffusion coefficients - Ternary and more complicated gas mixtures //
  //////////////////////////////////////////////////////////////////////////////////////////////////////

// ---                                                     ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dpressure ---
// ---                                                     ---

const Table< 2, double >
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dpressure(const double& total_pres) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       Table< 2, double > D(gases.size(), gases.size());
       D = get_ChapmanEnskog_isobaric_diffusion_coefficients();

       Table< 2, double > result(gases.size(), gases.size());

       for(unsigned int i = 0; i < gases.size(); ++i)
              for(unsigned int j = 0; j < gases.size(); ++j)
                     if( i != j )
                            result(i,j) = - D(i,j) / (total_pres*total_pres);

       return result;
}

// ---                                                     ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dpressure ---
// ---                                                     ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dpressure(const std::vector<double>&         total_pres,
                                                                      std::vector< Table< 2, double > >& dst) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );

       for(unsigned int q = 0; q < total_pres.size(); ++q)
       {
              dst[q].reinit(gases.size(), gases.size());
              dst[q] = get_DChapmanEnskog_diffusion_coefficients_Dpressure(total_pres[q]);
       }
}

// ---                                                     ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dpressure ---
// ---                                                     ---

const Table< 2, double >
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dpressure(const double& total_pres,
                                                                      const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       Table< 2, double > D(gases.size(), gases.size());
       D = get_ChapmanEnskog_isobaric_diffusion_coefficients(temp);

       Table< 2, double > result(gases.size(), gases.size());

       for(unsigned int i = 0; i < gases.size(); ++i)
              for(unsigned int j = 0; j < gases.size(); ++j)
                     if( i != j )
                            result(i,j) = - D(i,j) / (total_pres*total_pres);

       return result;
}

// ---                                                     ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dpressure ---
// ---                                                     ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dpressure(const std::vector<double>&         total_pres,
                                                                      const std::vector<double>&         temp,
                                                                      std::vector< Table< 2, double > >& dst) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );
       AssertThrow( temp.size()       == dst.size() , ExcDimensionMismatch(temp.size(),       dst.size()) );

       for(unsigned int q = 0; q < temp.size(); ++q)
       {
              dst[q].reinit(gases.size(), gases.size());
              dst[q] = get_DChapmanEnskog_diffusion_coefficients_Dpressure(total_pres[q],
                                                                           temp[q]);
       }
}

// ---                                                        ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dtemperature ---
// ---                                                        ---

const Table< 2, double >
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       Table< 2, double > DD_DT(gases.size(), gases.size());
       DD_DT = get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(temp);

       Table< 2, double > result(gases.size(), gases.size());

       for(unsigned int i = 0; i < gases.size(); ++i)
              for(unsigned int j = 0; j < gases.size(); ++j)
                     if( i != j )
                            result(i,j) = DD_DT(i,j) / total_pressure;

       return result;
}

// ---                                                        ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dtemperature ---
// ---                                                        ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const std::vector<double>&         temp,
                                                                         std::vector< Table< 2, double > >& dst) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( total_pressure != _DUMMY_,
                    ExcMessage("The mixture is NOT isobaric. Use another function.") );

       AssertThrow( temp.size() == dst.size() , ExcDimensionMismatch(temp.size(), dst.size()) );

       for(unsigned int q = 0; q < temp.size(); ++q)
       {
              dst[q].reinit(gases.size(), gases.size());
              dst[q] = get_DChapmanEnskog_diffusion_coefficients_Dtemperature(temp[q]);
       }
}

// ---                                                        ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dtemperature ---
// ---                                                        ---

const Table< 2, double >
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const double& total_pres,
                                                                         const double& temp) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       Table< 2, double > DD_DT(gases.size(), gases.size());
       DD_DT = get_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature(temp);

       Table< 2, double > result(gases.size(), gases.size());

       for(unsigned int i = 0; i < gases.size(); ++i)
              for(unsigned int j = 0; j < gases.size(); ++j)
                     if( i != j )
                            result(i,j) = DD_DT(i,j) / total_pres;

       return result;
}

// ---                                                        ---
// --- get_DChapmanEnskog_diffusion_coefficients_Dtemperature ---
// ---                                                        ---

void
NAME::GasMixture::get_DChapmanEnskog_diffusion_coefficients_Dtemperature(const std::vector<double>&         total_pres,
                                                                         const std::vector<double>&         temp,
                                                                         std::vector< Table< 2, double > >& dst) const
{
       AssertThrow( !gases.empty(),
                    ExcMessage("std::vector gases does not contain data. Fill it up !") );

       AssertThrow( total_pres.size() == dst.size() , ExcDimensionMismatch(total_pres.size(), dst.size()) );
       AssertThrow( temp.size()       == dst.size() , ExcDimensionMismatch(temp.size(),       dst.size()) );

       for(unsigned int q = 0; q < temp.size(); ++q)
       {
              dst[q].reinit(gases.size(), gases.size());
              dst[q] = get_DChapmanEnskog_diffusion_coefficients_Dtemperature(total_pres[q],
                                                                              temp[q]);
       }
}

  ///////////////////////////////
  // Binary collision integral //
  ///////////////////////////////

// ---                               ---
// --- get_binary_collision_integral ---
// ---                               ---

const double
NAME::GasMixture::get_binary_collision_integral(const unsigned int& N1,
                                                const unsigned int& N2) const
{
       AssertThrow( gases.size() > 1,
                    ExcMessage("The mixture MUST contain at least 2 gases") );

       AssertThrow( std::max(N1, N2) < gases.size(),
                    ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       const double epsN1N2_BY_k = std::sqrt( gases[N1]->get_eps_BY_k() * gases[N2]->get_eps_BY_k() );

       const double tmp = temperature/epsN1N2_BY_k;

       return Constants::A_diff() / std::pow(tmp, Constants::B_diff())
              +
              Constants::C_diff() / std::exp(Constants::D_diff()*tmp)
              +
              Constants::E_diff() / std::exp(Constants::F_diff()*tmp)
              +
              Constants::G_diff() / std::exp(Constants::H_diff()*tmp);
}

// ---                               ---
// --- get_binary_collision_integral ---
// ---                               ---

void
NAME::GasMixture::get_binary_collision_integral(std::vector<double>& binary_collision_integral,
                                                const unsigned int&  N1,
                                                const unsigned int&  N2) const
{
       AssertThrow( gases.size() > 1,
                    ExcMessage("The mixture MUST contain at least 2 gases") );

       AssertThrow( std::max(N1, N2) < gases.size(),
                    ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );

       AssertThrow( temperature != _DUMMY_,
                    ExcMessage("The mixture is NOT isothermal. Use another function.") );

       const double result = get_binary_collision_integral(N1,N2);

       for(unsigned int q = 0; q < binary_collision_integral.size(); ++q)
              binary_collision_integral[q] = result;
}

// ---                               ---
// --- get_binary_collision_integral ---
// ---                               ---

const double
NAME::GasMixture::get_binary_collision_integral(const double&       temp,
                                                const unsigned int& N1,
                                                const unsigned int& N2) const
{
       AssertThrow( gases.size() > 1,
                    ExcMessage("The mixture MUST contain at least 2 gases") );

       AssertThrow( std::max(N1, N2) < gases.size(),
                    ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );

       const double epsN1N2_BY_k = std::sqrt( gases[N1]->get_eps_BY_k() * gases[N2]->get_eps_BY_k() );

       const double tmp = temp/epsN1N2_BY_k;

       return Constants::A_diff() / std::pow(tmp, Constants::B_diff())
              +
              Constants::C_diff() / std::exp(Constants::D_diff()*tmp)
              +
              Constants::E_diff() / std::exp(Constants::F_diff()*tmp)
              +
              Constants::G_diff() / std::exp(Constants::H_diff()*tmp);
}

// ---                               ---
// --- get_binary_collision_integral ---
// ---                               ---

void
NAME::GasMixture::get_binary_collision_integral(const std::vector<double>& temp,
                                                std::vector<double>&       binary_collision_integral,
                                                const unsigned int&        N1,
                                                const unsigned int&        N2) const
{
       AssertThrow( gases.size() > 1,
                    ExcMessage("The mixture MUST contain at least 2 gases") );

       AssertThrow( std::max(N1, N2) < gases.size(),
                    ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );

       AssertThrow( temp.size() == binary_collision_integral.size() , ExcDimensionMismatch(temp.size(), binary_collision_integral.size()) );

       for(unsigned int q = 0; q < temp.size(); ++q)
              binary_collision_integral[q] = get_binary_collision_integral(temp[q],
                                                                           N1,N2);
}

// ---                                             ---
// --- get_Dbinary_collision_integral_Dtemperature ---
// ---                                             ---

const double
NAME::GasMixture::get_Dbinary_collision_integral_Dtemperature(const double&       temp,
                                                              const unsigned int& N1,
                                                              const unsigned int& N2) const
{
       AssertThrow( gases.size() > 1,
                    ExcMessage("The mixture MUST contain at least 2 gases") );

       AssertThrow( std::max(N1, N2) < gases.size(),
                    ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );

       const double epsN1N2_BY_k = std::sqrt( gases[N1]->get_eps_BY_k() * gases[N2]->get_eps_BY_k() );

       const double tmp = temp/epsN1N2_BY_k;

       return - (1.0/epsN1N2_BY_k)*(   Constants::A_diff()*Constants::B_diff()*std::pow( tmp, (-Constants::B_diff()-1.0) )
                                       +
                                       Constants::C_diff()*Constants::D_diff()*std::exp( -Constants::D_diff()*tmp )
                                       +
                                       Constants::E_diff()*Constants::F_diff()*std::exp( -Constants::F_diff()*tmp )
                                       +
                                       Constants::G_diff()*Constants::H_diff()*std::exp( -Constants::H_diff()*tmp )   );
}

// ---                                             ---
// --- get_Dbinary_collision_integral_Dtemperature ---
// ---                                             ---

void
NAME::GasMixture::get_Dbinary_collision_integral_Dtemperature(const std::vector<double>& temp,
                                                              std::vector<double>&       dst,
                                                              const unsigned int&        N1,
                                                              const unsigned int&        N2) const
{
       AssertThrow( gases.size() > 1,
                    ExcMessage("The mixture MUST contain at least 2 gases") );

       AssertThrow( std::max(N1, N2) < gases.size(),
                    ExcMessage("The maximum of two last arguments MUST be smaller than the total number of gases") );

       AssertThrow( temp.size() == dst.size() , ExcDimensionMismatch(temp.size(), dst.size()) );

       for(unsigned int q = 0; q < temp.size(); ++q)
              dst[q] = get_Dbinary_collision_integral_Dtemperature(temp[q],
                                                                   N1,N2);
}

       ////////////////////////
       ////////////////////////
       // ACCESSORS AND INFO //
       ////////////////////////
       ////////////////////////

// ---                           ---
// --- print_material_properties ---
// ---                           ---

void
NAME::GasMixture::print_material_properties() const
{
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "------------------------------";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;

  FcstUtilities::log << "Parameters for pure gases :";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;

  for(unsigned int g = 0; g < gases.size(); ++g)
  {
         FcstUtilities::log << gases[g]->name_material() << " :";
         FcstUtilities::log << std::endl;

         FcstUtilities::log << "ID = "                                                                      << gases[g]->get_ID();
         FcstUtilities::log << std::endl;

         FcstUtilities::log << "Chemical formula = "                                                        << gases[g]->get_chemical_formula();
         FcstUtilities::log << std::endl;

         FcstUtilities::log << "Molar mass, [kg/mol] = "                                                    << gases[g]->get_molar_mass();
         FcstUtilities::log << std::endl;

         FcstUtilities::log << "Collision diameter, [A] = "                                                 << gases[g]->get_collision_diameter();
         FcstUtilities::log << std::endl;

         FcstUtilities::log << "The maximum energy of attraction divided by the Boltzmann constant, [K] = " << gases[g]->get_eps_BY_k();
         FcstUtilities::log << std::endl;

         FcstUtilities::log << "Prandtl number = "                                                          << gases[g]->get_Prandtl();
         FcstUtilities::log << std::endl;

         FcstUtilities::log << "Dynamic viscosity mode = "                                                  << gases[g]->get_dynamic_viscosity_mode();
         FcstUtilities::log << std::endl;

         FcstUtilities::log << "Bulk viscosity mode = "                                                     << gases[g]->get_bulk_viscosity_mode();
         FcstUtilities::log << std::endl;

         FcstUtilities::log << "Thermal conductivity mode = "                                               << gases[g]->get_thermal_conductivity_mode();
         FcstUtilities::log << std::endl;

         if( temperature != _DUMMY_ )
         {
                FcstUtilities::log << "Dynamic viscosity, [Pa sec] = "                             << gases[g]->get_dynamic_viscosity(temperature);
                FcstUtilities::log << std::endl;

                FcstUtilities::log << "Bulk viscosity, [Pa sec] = "                                << gases[g]->get_bulk_viscosity(   gases[g]->get_dynamic_viscosity(temperature)   );
                FcstUtilities::log << std::endl;

                FcstUtilities::log << "Thermal conductivity, [W/(m K)] = "                         << gases[g]->get_thermal_conductivity(temperature);
                FcstUtilities::log << std::endl;

                FcstUtilities::log << "Specific heat capacity at constant pressure, [J/(kg K)] = " << gases[g]->get_specific_heat_capacity(temperature);
                FcstUtilities::log << std::endl;

                FcstUtilities::log << "Molar enthalpy, [J/mol] = "                                 << gases[g]->get_molar_enthalpy(temperature);
                FcstUtilities::log << std::endl;

                if( gases[g]->get_chemical_formula() == "H2O Vapor" )
                {
                       FcstUtilities::log << "Water vapor saturation pressure, [Pa] = "            << gases[g]->get_water_vapor_saturation_pressure(temperature);
                       FcstUtilities::log << std::endl;
                }
         }

         FcstUtilities::log << std::endl;
  }

  FcstUtilities::log << "Parameters for the gas mixture :";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;

  FcstUtilities::log << this->name << " :";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;

  if( temperature != _DUMMY_ && total_pressure != _DUMMY_ )
         FcstUtilities::log << "The gas mixture is both ISOTHERMAL and ISOBARIC"         << std::endl;

  if( temperature == _DUMMY_ && total_pressure != _DUMMY_ )
         FcstUtilities::log << "The gas mixture is both NON-ISOTHERMAL and ISOBARIC"     << std::endl;

  if( temperature != _DUMMY_ && total_pressure == _DUMMY_ )
         FcstUtilities::log << "The gas mixture is both ISOTHERMAL and NON-ISOBARIC"     << std::endl;

  if( temperature == _DUMMY_ && total_pressure == _DUMMY_ )
         FcstUtilities::log << "The gas mixture is both NON-ISOTHERMAL and NON-ISOBARIC" << std::endl;

  if( temperature != _DUMMY_ )
  {
         FcstUtilities::log << "Chapman Enskog isobaric diffusion coefficients, [Pa m^2/sec] :";
         FcstUtilities::log << std::endl;

         Table< 2, double > D(gases.size(), gases.size());
         D = get_ChapmanEnskog_isobaric_diffusion_coefficients();

         for(unsigned int i = 0; i < gases.size(); ++i)
         {
                for(unsigned int j = 0; j < gases.size(); ++j)
                       FcstUtilities::log << D(i,j) << "   ";
                FcstUtilities::log << std::endl;
         }

         if( total_pressure != _DUMMY_ )
         {
                FcstUtilities::log << "Chapman Enskog diffusion coefficients, [m^2/sec] :";
                FcstUtilities::log << std::endl;

                D = get_ChapmanEnskog_diffusion_coefficients();

                for(unsigned int i = 0; i < gases.size(); ++i)
                {
                       for(unsigned int j = 0; j < gases.size(); ++j)
                              FcstUtilities::log << D(i,j) << "   ";
                       FcstUtilities::log << std::endl;
                }
         }
  }
}