// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: experimental_multi_fluid.cc
// - Description: This class describes a multi fluid
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: experimental_multi_fluid.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include "experimental_multi_fluid.h"

namespace NAME = FuelCellShop::Material;

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

NAME::ExperimentalMultiFluid::ExperimentalMultiFluid(const std::string& name)
:
NAME::BaseMaterial(name)
{ }

// ---            ---
// --- Destructor ---
// ---            ---

NAME::ExperimentalMultiFluid::~ExperimentalMultiFluid()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

void
NAME::ExperimentalMultiFluid::declare_parameters(ParameterHandler& param) const
{
  const unsigned int n_species_max = 5;

  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection("Materials");
    {
      param.enter_subsection(this->name);
      {
             param.declare_entry("Number of species",
                                 "1",
                                  Patterns::Integer(1,5),
                                 " ");

             param.declare_entry("Temperature of mixture",
                                 "300.0",
                                  Patterns::Double(),
                                 " ");

             for(unsigned int index = 1; index <= n_species_max; ++index)
             {
                    std::ostringstream streamOut;
                    streamOut << index;
                    const std::string name = "Molar mass " + streamOut.str();

                    param.declare_entry(name.c_str(),
                                       "1.0",
                                        Patterns::Double(),
                                       " ");
             }

             for(unsigned int index = 1; index <= n_species_max; ++index)
             {
                    std::ostringstream streamOut;
                    streamOut << index;
                    const std::string name = "Dynamic viscosity " + streamOut.str();

                    param.declare_entry(name.c_str(),
                                       "1.0e-2",
                                        Patterns::Double(),
                                       " ");
             }

             for(unsigned int index = 1; index <= n_species_max; ++index)
             {
                    std::ostringstream streamOut;
                    streamOut << index;
                    const std::string name = "Bulk viscosity " + streamOut.str();

                    param.declare_entry(name.c_str(),
                                       "0.0",
                                        Patterns::Double(),
                                       " ");
             }

             for(unsigned int line = 1; line <= n_species_max; ++line)
             {
                    for(unsigned int column = 1; column <= n_species_max; ++column)
                    {
                           std::ostringstream streamOut;
                           streamOut << line << column;
                           const std::string name = "Maxwell-Stefan isobaric diffusion coefficient " + streamOut.str();

                           param.declare_entry(name.c_str(),
                                              "0.0",
                                               Patterns::Double(),
                                              " ");
                    }
             }
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
NAME::ExperimentalMultiFluid::initialize(ParameterHandler& param)
{
  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection("Materials");
    {
      param.enter_subsection(this->name);
      {
             n_species = param.get_integer("Number of species");
             T_mixture = param.get_double("Temperature of mixture");

             for(unsigned int index = 1; index <= n_species; ++index)
             {
                    std::ostringstream streamOut;
                    streamOut << index;
                    const std::string name = "Molar mass " + streamOut.str();

                    molar_mass.push_back( param.get_double(name.c_str()) );
             }

             for(unsigned int index = 1; index <= n_species; ++index)
             {
                    std::ostringstream streamOut;
                    streamOut << index;
                    const std::string name = "Dynamic viscosity " + streamOut.str();

                    dynamic_viscosity.push_back( param.get_double(name.c_str()) );
             }

             for(unsigned int index = 1; index <= n_species; ++index)
             {
                    std::ostringstream streamOut;
                    streamOut << index;
                    const std::string name = "Bulk viscosity " + streamOut.str();

                    bulk_viscosity.push_back( param.get_double(name.c_str()) );
             }

             maxwell_stefan_isobaric_diffusion_coefficient.reinit(n_species, n_species);

             for(unsigned int line = 1; line <= n_species; ++line)
             {
                    for(unsigned int column = 1; column <= n_species; ++column)
                    {
                           std::ostringstream streamOut;
                           streamOut << line << column;
                           const std::string name = "Maxwell-Stefan isobaric diffusion coefficient " + streamOut.str();

                           maxwell_stefan_isobaric_diffusion_coefficient(line-1, column-1) = param.get_double(name.c_str());
                    }
             }
      }
      param.leave_subsection();
    }
    param.leave_subsection();
  }
  param.leave_subsection();
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
NAME::ExperimentalMultiFluid::print_material_properties() const
{
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "------------------------------";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "Parameters for " << this->name << ":";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "Number of species:          " << n_species;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "Temperature of mixture [K]: " << T_mixture;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "Molar mass [kg/mol]:";
  FcstUtilities::log << std::endl;
  for(unsigned int index = 0; index < n_species; ++index)
    FcstUtilities::log << molar_mass[index] << std::endl;
  FcstUtilities::log << "Dynamic viscosity [Pa sec]:";
  FcstUtilities::log << std::endl;
  for(unsigned int index = 0; index < n_species; ++index)
    FcstUtilities::log << dynamic_viscosity[index] << std::endl;
  FcstUtilities::log << "Bulk viscosity [Pa sec]:";
  FcstUtilities::log << std::endl;
  for(unsigned int index = 0; index < n_species; ++index)
    FcstUtilities::log << bulk_viscosity[index] << std::endl;
  FcstUtilities::log << "Maxwell-Stefan isobaric diffusion coefficient [Pa m^2/sec]:";
  FcstUtilities::log << std::endl;
  for(unsigned int line = 0; line < n_species; ++line)
  {
    for(unsigned int column = 0; column < n_species; ++column)
      FcstUtilities::log << maxwell_stefan_isobaric_diffusion_coefficient(line, column) << "   ";
    FcstUtilities::log << std::endl;
  }
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "------------------------------";
  FcstUtilities::log << std::endl;
}