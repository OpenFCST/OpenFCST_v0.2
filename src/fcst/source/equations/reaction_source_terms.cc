// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: reaction_source_terms.cc
// - Description: This class is used to assemble both cell matrix and cell residual
//                for reaction source terms in the catalyst layers for various equation classes
// - Developers: Madhur Bhaiya, Marc Secanell, Valentin N. Zingan
// - Id: $Id$
//
// ----------------------------------------------------------------------------

#include "equations/reaction_source_terms.h"

namespace NAME = FuelCellShop::Equation;

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::ReactionSourceTerms<dim>::ReactionSourceTerms(FuelCell::SystemManagement& system_management)
:
NAME::EquationBase<dim>(system_management),
electronic_electrical_potential_extractor(0),
protonic_electrical_potential_extractor(0)
{
    FcstUtilities::log << "->FuelCellShop::Equation::ReactionSourceTerms" << std::endl;

    anode_kinetics = NULL;
    cathode_kinetics = NULL;

    //-- Initializing ReactionHeat objects to NULL pointers. This avoids checking while "delete" whether they are allocated using "new" or not --------
    anode_reactionheat = NULL;
    cathode_reactionheat = NULL;

    //----Initializing VariableInfo Structs------------------------------------------
    //----Setting indices_exist to false --------------------------------------------
    x_oxygen.indices_exist = false;
    x_water.indices_exist = false;
    phi_s.indices_exist = false;
    phi_m.indices_exist = false;
    lambda.indices_exist = false;
    t_rev.indices_exist = false;
    s_liquid_water.indices_exist = false;

    this->counter.resize(4, true);

    indexO2  = -1; // does not exist in the problem
    indexH2  = -1; // does not exist in the problem
    indexH2O = -1; // does not exist in the problem
}

// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::ReactionSourceTerms<dim>::~ReactionSourceTerms()
{
    delete cathode_reactionheat;
    delete anode_reactionheat;
}

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::declare_parameters(ParameterHandler& param) const
{
    param.enter_subsection("Reaction Source Terms");
    {
        param.declare_entry("Irreversible heat source due to ORR",
                            "true",
                            Patterns::Bool(),
                            "Flag to include irreversible heating due to ORR inside the cathode catalyst layer.");
        param.declare_entry("Irreversible heat source due to HOR",
                            "true",
                            Patterns::Bool(),
                            "Flag to include irreversible heating due to HOR inside the anode catalyst layer.");
        param.declare_entry("Reversible heat source due to net reaction",
                            "true",
                            Patterns::Bool(),
                            "Flag to include reversible heating due to net reaction forming liquid water product.");
        param.declare_entry("Reversible heat fraction in ORR",
                            "1.0",
                            Patterns::Double(),
                            "Fraction of reversible heat released due to half-cell reaction of ORR inside the cathode catalyst layer.");
        param.declare_entry("Water produced during ORR in vapour phase",
                            "true",
                            Patterns::Bool(),
                            "Flag considering whether the water produced during ORR is in vapour phase or not.");
        param.declare_entry("Water vaporization heat sink in CCL",
                            "true",
                            Patterns::Bool(),
                            "If complete evaporation of water produced during ORR is considered, this flag is set to account for heat sink due to vaporization.");
    }
    param.leave_subsection();
}

// ---            ---
// --- initialize ---
// ---            ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::initialize(ParameterHandler& param)
{
    param.enter_subsection("Reaction Source Terms");
    {
        irrev_heat_ccl = param.get_bool("Irreversible heat source due to ORR");
        irrev_heat_acl = param.get_bool("Irreversible heat source due to HOR");
        rev_heat = param.get_bool("Reversible heat source due to net reaction");
        factor_rev_heat_ccl = param.get_double("Reversible heat fraction in ORR");
        water_vapour_phase = param.get_bool("Water produced during ORR in vapour phase");
        water_vap_heat_ccl = param.get_bool("Water vaporization heat sink in CCL");
    }
    param.leave_subsection();

    Assert( cathode_kinetics!=NULL || anode_kinetics!=NULL, ExcMessage("Atleast one of either cathode/anode kinetics should be set before ReactionSourceTerms::initialize.") );

    //----- Assertion checks whether water is produced in liquid phase or vapour phase (complete evaporated) in the CCL ----------------
    if (water_vapour_phase == false)
    {
        AssertThrow( this->system_management->solution_in_userlist("liquid_water_saturation"),
                    ExcMessage("Liquid water saturation transport equation should be considered, when \"Water produced during ORR in vapour phase\" is set to FALSE.") );
        AssertThrow( water_vap_heat_ccl == false, ExcMessage("\"Water vaporization heat sink in CCL\" should be set to FALSE, when \"Water produced during ORR in vapour phase\" is set to FALSE.") );
        AssertThrow( cathode_kinetics != NULL, ExcMessage("Cathode kinetics not initialized in the ReactionSourceTerms class for considering liquid water product in the CCL.") );
    }

    //-------Assertion checks for "Required" solution variables for Cathode/Anode kinetics --------------------
    AssertThrow( this->system_management->solution_in_userlist("electronic_electrical_potential"), VariableNotFoundForKinetics("Cathode/Anode", "electronic_electrical_potential") );
    AssertThrow( this->system_management->solution_in_userlist("protonic_electrical_potential"), VariableNotFoundForKinetics("Cathode/Anode", "protonic_electrical_potential") );

    //if ( cathode_kinetics != NULL )
        //AssertThrow( this->system_management->solution_in_userlist("oxygen_molar_fraction"), VariableNotFoundForKinetics("Cathode", "oxygen_molar_fraction") );

    //if ( anode_kinetics != NULL )
        //AssertThrow( this->system_management->solution_in_userlist("water_molar_fraction"), VariableNotFoundForKinetics("Anode", "water_molar_fraction") );

    //-------------- Initializing reaction heat objects if temperature is being solved for ---------------------------------------------
    if ( this->system_management->solution_in_userlist("temperature_of_REV") )
    {
        if ( cathode_kinetics != NULL )
        {
            cathode_reactionheat = new FuelCellShop::Equation::ReactionHeat;
            cathode_reactionheat->set_kinetics(cathode_kinetics);
            cathode_reactionheat->initialize_factors(irrev_heat_ccl, irrev_heat_acl, rev_heat, factor_rev_heat_ccl, water_vap_heat_ccl);
        }

        if ( anode_kinetics != NULL )
        {
            anode_reactionheat = new FuelCellShop::Equation::ReactionHeat;
            anode_reactionheat->set_kinetics(anode_kinetics);
            anode_reactionheat->initialize_factors(irrev_heat_ccl, irrev_heat_acl, rev_heat, factor_rev_heat_ccl, water_vap_heat_ccl);
        }
    }

    //----Setting cell_matrix and cell_residual counters to FALSE -------------------
    cell_matrix_counter   = false;
    cell_residual_counter = false;
}

// ---                                ---
// --- adjust_internal_cell_couplings ---
// ---                                ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::adjust_internal_cell_couplings(std::vector< couplings_map >& equation_map,
                                                               const std::vector< FuelCellShop::Material::PureGas* >& gases)
{
  Assert( equation_map.size() != 0, ExcMessage("Vector size should be greater than zero in ReactionSourceTerms::adjust_internal_cell_couplings") );

  if( gases.size() == 0 )
  {
    for (unsigned int i=0; i<equation_map.size(); ++i)
    {
        for ( couplings_map::iterator iter = equation_map[i].begin(); iter != equation_map[i].end(); ++iter )
        {
            if ( iter->first == "Thermal Transport Equation" )
            {
                (iter->second)["temperature_of_REV"] = DoFTools::always;
                (iter->second)["electronic_electrical_potential"] = DoFTools::always;
                (iter->second)["protonic_electrical_potential"] = DoFTools::always;

                if (cathode_kinetics != NULL)
                    (iter->second)["oxygen_molar_fraction"] = DoFTools::always;

                if (anode_kinetics != NULL)
                    (iter->second)["water_molar_fraction"] = DoFTools::always;
            }
            else if ( (iter->first == "Proton Transport Equation") || (iter->first == "Electron Transport Equation") )
            {
                (iter->second)["electronic_electrical_potential"] = DoFTools::always;
                (iter->second)["protonic_electrical_potential"] = DoFTools::always;

                if ( this->system_management->solution_in_userlist("temperature_of_REV") )
                    (iter->second)["temperature_of_REV"] = DoFTools::always;

                if (cathode_kinetics != NULL)
                    (iter->second)["oxygen_molar_fraction"] = DoFTools::always;

                if (anode_kinetics != NULL)
                    (iter->second)["water_molar_fraction"] = DoFTools::always;
            }
            else if ( iter->first == "Ficks Transport Equation - oxygen" )
            {
                (iter->second)["oxygen_molar_fraction"] = DoFTools::always;
                (iter->second)["electronic_electrical_potential"] = DoFTools::always;
                (iter->second)["protonic_electrical_potential"] = DoFTools::always;

                if (anode_kinetics != NULL)
                    (iter->second)["water_molar_fraction"] = DoFTools::always;

                if ( this->system_management->solution_in_userlist("temperature_of_REV") )
                    (iter->second)["temperature_of_REV"] = DoFTools::always;
            }
            else if ( ( iter->first == "Ficks Transport Equation - water" ) && water_vapour_phase )
            {
                (iter->second)["water_molar_fraction"] = DoFTools::always;
                (iter->second)["electronic_electrical_potential"] = DoFTools::always;
                (iter->second)["protonic_electrical_potential"] = DoFTools::always;

                if (cathode_kinetics != NULL)
                    (iter->second)["oxygen_molar_fraction"] = DoFTools::always;

                if ( this->system_management->solution_in_userlist("temperature_of_REV") )
                    (iter->second)["temperature_of_REV"] = DoFTools::always;
            }
            else if ( ( iter->first == "Liquid Water Saturation Transport Equation" ) && !water_vapour_phase )
            {
                (iter->second)["oxygen_molar_fraction"] = DoFTools::always;
                (iter->second)["water_molar_fraction"] = DoFTools::always;
                (iter->second)["electronic_electrical_potential"] = DoFTools::always;
                (iter->second)["protonic_electrical_potential"] = DoFTools::always;
                (iter->second)["temperature_of_REV"] = DoFTools::always;
            }
        }
    }
  }

  ////////////////////////////////////////////////////////////
  // In this case, the core of a problem is formed based on //
  // Kerkhof equations rather than Ficks equations          //
  ////////////////////////////////////////////////////////////

  else
  {
         AssertThrow( gases.size() <= 5,
                      ExcMessage("The total number of gases used in this problem can not exceed 5") );

         /////////////////////////////////////
         // LOTS OF COMPUTATIONAL CONSTANTS //
         // ARE ASSEMBLED HERE              //
         /////////////////////////////////////

         ///////
         // 1 //
         ///////

         n_species = gases.size();

         ///////
         // 2 //
         ///////

         for(unsigned int g = 1; g <= n_species; ++g)
         {
                if( gases[g-1]->name_material() == "oxygen"   )
                {
                       indexO2  = g;
                }

                if( gases[g-1]->name_material() == "hydrogen" )
                {
                       indexH2  = g;
                }

                if( gases[g-1]->name_material() == "water"    )
                {
                       indexH2O = g;
                }
         }

         ///////
         // 3 //
         ///////

         molar_mass.resize(n_species);

         for(unsigned int g = 0; g < n_species; ++g)
                molar_mass[g] = gases[g]->get_molar_mass() * 1.0e3;

         ///////
         // 4 //
         ///////

         if( indexO2  != -1 )
                multiplierO2  = - molar_mass[indexO2-1]  / ( 4.0 * Constants::F() );
         if( indexH2  != -1 )
                multiplierH2  = - molar_mass[indexH2-1]  / ( 2.0 * Constants::F() );
         if( indexH2O != -1 )
                multiplierH2O =   molar_mass[indexH2O-1] / ( 2.0 * Constants::F() );

         ///////
         // 5 //
         ///////

         std::vector<unsigned int> density_indices(n_species);

         for(unsigned int g = 0; g < n_species; ++g)
                density_indices[g] = g*(dim+1) + this->system_management->solution_name_to_index("density_species_1");

         for(unsigned int g = 0; g < n_species; ++g)
                density_extractors.push_back( FEValuesExtractors::Scalar(density_indices[g]) );

         electronic_electrical_potential_extractor.component = this->system_management->solution_name_to_index("electronic_electrical_potential");
         protonic_electrical_potential_extractor.component   = this->system_management->solution_name_to_index("protonic_electrical_potential");

         ///////
         // 7 //
         ///////

         eq_generic_prefix = "Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - ";

         ///////
         // 8 //
         ///////

         for(unsigned int g = 1; g <= n_species; ++g)
         {
                std::ostringstream streamOut;
                streamOut << g;
                const std::string name = "species " + streamOut.str();

                eq_postfixes.push_back( name.c_str() );
         }

         ///////
         // 9 //
         ///////

         for(unsigned int g = 1; g <= n_species; ++g)
         {
                std::ostringstream streamOut;
                streamOut << g;
                const std::string name = "species_" + streamOut.str();

                var_postfixes.push_back( name.c_str() );
         }

         //////////////////////////////////////////////////////////
         // WE WANT TO FIND WHICH ENTRY OF equation_map CONTAINS //
         // THE INTERNAL CELL COUPLINGS RELATED TO KERKHOF       //
         // EQUATIONS                                            //
         //////////////////////////////////////////////////////////

         unsigned int index; // keeps Kerkhof equations in equation_map
         bool         flag = false;

         for(unsigned int i = 0; i < equation_map.size(); ++i)
         {
                const couplings_map tmp = equation_map[i];

                couplings_map::const_iterator iter = tmp.find("Kerkhof-Geboers Fluid Transport Equations - steady-state - compressible - isothermal - single-phase - multi-component - mass conservation - species 1");

                if( iter != tmp.end() )
                {
                       index = i;
                       flag  = true;
                       break;
                }
         }

         AssertThrow( flag,
                      ExcMessage("Kerkhof equations are not found in this problem") );

         /////////////
         // KERKHOF //
         /////////////

         for(unsigned int g = 1; g <= n_species; ++g)
         {
                if( g == indexO2 )
                {
                       if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
                           ||
                           cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
                       {
                              //////////////////////////////
                              // SET # g DESCRIBES OXYGEN //
                              //////////////////////////////

                              ///////////////////////
                              // MASS CONSERVATION //
                              ///////////////////////

                              eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];

                              ////////////////////////////////////////////
                              // LET'S INCLUDE COUPLINGS WITH Fs AND Fm //
                              ////////////////////////////////////////////

                              couplings_map::iterator iter = equation_map[index].find(eq_name);

                              if( iter != equation_map[index].end() )
                              {
                                     iter->second["electronic_electrical_potential"] = DoFTools::always;
                                     iter->second["protonic_electrical_potential"]   = DoFTools::always;
                              }
                              else
                              {
                                     AssertThrow( false, ExcInternalError() );
                              }
                       }
                }

                if( g == indexH2 )
                {
                       if( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR
                           ||
                           cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )
                       {
                              ////////////////////////////////
                              // SET # g DESCRIBES HYDROGEN //
                              ////////////////////////////////

                              ///////////////////////
                              // MASS CONSERVATION //
                              ///////////////////////

                              eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];

                              ////////////////////////////////////////////
                              // LET'S INCLUDE COUPLINGS WITH Fs AND Fm //
                              ////////////////////////////////////////////

                              couplings_map::iterator iter = equation_map[index].find(eq_name);

                              if( iter != equation_map[index].end() )
                              {
                                     iter->second["electronic_electrical_potential"] = DoFTools::always;
                                     iter->second["protonic_electrical_potential"]   = DoFTools::always;
                              }
                              else
                              {
                                     AssertThrow( false, ExcInternalError() );
                              }
                       }
                }

                if( g == indexH2O )
                {
                       if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
                           ||
                           cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
                       {
                              ///////////////////////////////////
                              // SET # g DESCRIBES WATER VAPOR //
                              ///////////////////////////////////

                              ///////////////////////
                              // MASS CONSERVATION //
                              ///////////////////////

                              eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];

                              ////////////////////////////////////////////
                              // LET'S INCLUDE COUPLINGS WITH Fs AND Fm //
                              ////////////////////////////////////////////

                              couplings_map::iterator iter = equation_map[index].find(eq_name);

                              if( iter != equation_map[index].end() )
                              {
                                     iter->second["electronic_electrical_potential"] = DoFTools::always;
                                     iter->second["protonic_electrical_potential"]   = DoFTools::always;
                              }
                              else
                              {
                                     AssertThrow( false, ExcInternalError() );
                              }

                              ////////////////////////////////////////////////////
                              // LET'S INCLUDE COUPLINGS WITH DENSITY OF OXYGEN //
                              ////////////////////////////////////////////////////

                              if( indexO2 != -1 )
                              {
                                     var_name = "density_" + var_postfixes[indexO2-1];

                                     iter->second[var_name] = DoFTools::always;
                              }
                              else
                              {
                                     AssertThrow( false, ExcInternalError() );
                              }
                       }
                }
         }

         /////////////////////
         // OTHER EQUATIONS //
         /////////////////////

         for(unsigned int i = 0; i < equation_map.size(); ++i)
                if( i != index )
                       for( couplings_map::iterator iter  = equation_map[i].begin();
                                                    iter != equation_map[i].end();
                                                  ++iter )
                       {
                              if( iter->first == "Electron Transport Equation" )
                              {
                                     /////////////////////////////////////
                                     // LET'S INCLUDE COUPLINGS WITH Fm //
                                     /////////////////////////////////////

                                     iter->second["protonic_electrical_potential"] = DoFTools::always;

                                     ////////////////////////////////////////////////////
                                     // LET'S INCLUDE COUPLINGS WITH DENSITY OF OXYGEN //
                                     ////////////////////////////////////////////////////

                                     if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
                                         ||
                                         cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
                                     {
                                            if( indexO2 != -1 )
                                            {
                                                   var_name = "density_" + var_postfixes[indexO2-1];

                                                   iter->second[var_name] = DoFTools::always;
                                            }
                                            else
                                            {
                                                   AssertThrow( false, ExcNotImplemented() );
                                            }
                                     }

                                     //////////////////////////////////////////////////////
                                     // LET'S INCLUDE COUPLINGS WITH DENSITY OF HYDROGEN //
                                     //////////////////////////////////////////////////////

                                     if( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR
                                         ||
                                         cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )
                                     {
                                            if( indexH2 != -1 )
                                            {
                                                   var_name = "density_" + var_postfixes[indexH2-1];

                                                   iter->second[var_name] = DoFTools::always;
                                            }
                                            else
                                            {
                                                   AssertThrow( false, ExcNotImplemented() );
                                            }
                                     }
                              }

                              if( iter->first == "Proton Transport Equation" )
                              {
                                     /////////////////////////////////////
                                     // LET'S INCLUDE COUPLINGS WITH Fs //
                                     /////////////////////////////////////

                                     iter->second["electronic_electrical_potential"] = DoFTools::always;

                                     ////////////////////////////////////////////////////
                                     // LET'S INCLUDE COUPLINGS WITH DENSITY OF OXYGEN //
                                     ////////////////////////////////////////////////////

                                     if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
                                         ||
                                         cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
                                     {
                                            if( indexO2 != -1 )
                                            {
                                                   var_name = "density_" + var_postfixes[indexO2-1];

                                                   iter->second[var_name] = DoFTools::always;
                                            }
                                            else
                                            {
                                                   AssertThrow( false, ExcNotImplemented() );
                                            }
                                     }

                                     //////////////////////////////////////////////////////
                                     // LET'S INCLUDE COUPLINGS WITH DENSITY OF HYDROGEN //
                                     //////////////////////////////////////////////////////

                                     if( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR
                                         ||
                                         cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )
                                     {
                                            if( indexH2 != -1 )
                                            {
                                                   var_name = "density_" + var_postfixes[indexH2-1];

                                                   iter->second[var_name] = DoFTools::always;
                                            }
                                            else
                                            {
                                                   AssertThrow( false, ExcNotImplemented() );
                                            }
                                     }
                              }
                       }
  }
}

// ---                      ---
// --- assemble_cell_matrix ---
// ---                      ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                                     const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                     FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
       ////////////////
       // TYPES INFO //
       ////////////////

       ///////////////////////////////////
       // ONLY ANODE OR CATHODE         //
       // CATALYST LAYERS CAN BE        //
       // ACCEPTED AS THE LAST ARGUMENT //
       // OF THIS FUNCTION              //
       ///////////////////////////////////

       const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
       const std::type_info& info          = layer->get_base_type();

       ////////////////////
       // CATALYST LAYER //
       ////////////////////

       if( info == CatalystLayer )
       {
              try
              {
                     FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);

                     ////////////////////////////////////////////////////////
                     //                                                    //
                     // WE FIRST CHECK IF gas_mixture IS A NULL PTR OR NOT //
                     //                                                    //
                     // - IF IT IS A NULL PTR THEN THE CORE OF             //
                     //   MASS/MOMENTUM/ENERGY TRANSPORT IS FORMED AROUND  //
                     //   FICKS EQUATIONS                                  //
                     //                                                    //
                     // - IF IT IS NOT A NULL PTR THEN THE CORE OF         //
                     //   MASS/MOMENTUM/ENERGY TRANSPORT IS FORMED AROUND  //
                     //   KERKHOF EQUATIONS                                //
                     //                                                    //
                     ////////////////////////////////////////////////////////

                     const FuelCellShop::Material::GasMixture* gas_mixture = ptr->get_gas_mixture();

                     /////////////////////////////
                     // FICKS + OTHER EQUATIONS //
                     /////////////////////////////

                     if( !gas_mixture )
                     {
                            if ( this->counter[0] )
                            {
                                   this->make_assemblers_generic_constant_data();
                                   this->counter[0] = false;
                            }

                            if ( this->counter[1] )
                            {
                                   this->make_assemblers_cell_constant_data(cell_info);
                                   this->counter[1] = false;
                            }

                            cell_matrix_counter   = true;
                            cell_residual_counter = false;
                            this->make_assemblers_cell_variable_data(cell_info, layer);

                            assemble_matrix_for_equation(cell_matrices, cell_info, "Proton Transport Equation", cell_info.fe(phi_m.fetype_index), phi_phiM_cell, factor_protontranseq_cell);
                            assemble_matrix_for_equation(cell_matrices, cell_info, "Electron Transport Equation", cell_info.fe(phi_s.fetype_index), phi_phiS_cell, factor_electrontranseq_cell);
                            if ( x_oxygen.indices_exist )
                                   assemble_matrix_for_equation(cell_matrices, cell_info, "Ficks Transport Equation - oxygen", cell_info.fe(x_oxygen.fetype_index), phi_xOxygen_cell, factor_oxygentranseq_cell);
                            if ( water_vapour_phase && x_water.indices_exist )
                                   assemble_matrix_for_equation(cell_matrices, cell_info, "Ficks Transport Equation - water", cell_info.fe(x_water.fetype_index), phi_xWater_cell, factor_watertranseq_cell);
                            if ( !water_vapour_phase && s_liquid_water.indices_exist )
                                   assemble_matrix_for_equation(cell_matrices, cell_info, "Liquid Water Saturation Transport Equation", cell_info.fe(s_liquid_water.fetype_index), phi_s_cell, factor_saturationtranseq_cell);
                            if ( t_rev.indices_exist )
                                   assemble_matrix_for_equation(cell_matrices, cell_info, "Thermal Transport Equation", cell_info.fe(t_rev.fetype_index), phi_T_cell, -1.0);
                     }

                     ///////////////////////////////
                     // KERKHOF + OTHER EQUATIONS //
                     ///////////////////////////////

                     else
                     {
                            ///////////////////////////
                            // GENERIC CONSTANT DATA //
                            ///////////////////////////

                            if( this->counter[2] )
                            {
                                   this->make_assemblers_generic_constant_data(cell_info,
                                                                               layer);
                                   this->counter[2] = false;
                            }

                            ////////////////////////
                            // CELL CONSTANT DATA //
                            ////////////////////////

                            if( this->counter[3] )
                            {
                                   this->make_assemblers_cell_constant_data2(cell_info);
                                   this->counter[3] = false;
                            }

                            ////////////////////////
                            // CELL VARIABLE DATA //
                            ////////////////////////
                            this->make_assemblers_cell_variable_data2(cell_info,
                                                                             layer);

                            ///////////////////////
                            // LOCAL CELL MATRIX //
                            ///////////////////////

                            FullMatrix<double> local_matrix(this->dofs_per_cell, this->dofs_per_cell);

                            //////////////////////////////////////////////////
                            // SOURCE TERMS FOR MASS CONSERVATION EQUATIONS //
                            // OF OXYGEN, HYDROGEN, AND WATER VAPOR         //
                            // INSIDE ANODE OR/AND CATHODE CATALYST LAYERS  //
                            //////////////////////////////////////////////////

                            /////////// ///////////// /////////
                            // ANODE // // CATHODE // // HOR //
                            /////////// ///////////// /////////

                            if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == HOR
                                ||
                                ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == HOR )
                            {
                                   ///////////////////////
                                   // LOOP OVER SPECIES //
                                   ///////////////////////

                                   for(unsigned int s = 0; s < n_species; ++s)
                                   {
                                          if( s+1 == indexH2 )
                                          {
                                                 ///////////////////////////////// //////////////////////////////
                                                 // LOOP OVER QUADRATURE POINTS // // MASS CONSERVATION for H2 //
                                                 ///////////////////////////////// //////////////////////////////

                                                 for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                                 {
                                                        /////////////////
                                                        // LOOP OVER i //
                                                        /////////////////

                                                        for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                                        {
                                                               /////////////////
                                                               // LOOP OVER j //
                                                               /////////////////

                                                               for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                                                               {
                                                                      local_matrix(i,j) += - phi_density[s][q][i]*multiplierH2*( DHOR_current_density_Dhydrogen_concentration[q]*(1.0/molar_mass[s])*phi_density[s][q][j]
                                                                                                                                 +
                                                                                                                                 DHOR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                                                                                                                                 +
                                                                                                                                 DHOR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] )
                                                                                                                                 *
                                                                                                                                 this->JxW_cell[q];
                                                               }
                                                        }
                                                 }
                                          }
                                   }

                                   ///////////////////////////////// /////////////////////////////////
                                   // LOOP OVER QUADRATURE POINTS // // ELECTRON TRANSPORT EQUATION //
                                   ///////////////////////////////// /////////////////////////////////

                                   for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                   {
                                          /////////////////
                                          // LOOP OVER i //
                                          /////////////////

                                          for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                          {
                                                 /////////////////
                                                 // LOOP OVER j //
                                                 /////////////////

                                                 for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                                                 {
                                                        local_matrix(i,j) += - phi_electronic_electrical_potential[q][i]*(-1.0)*( DHOR_current_density_Dhydrogen_concentration[q]*(1.0/molar_mass[indexH2-1])*phi_density[indexH2-1][q][j]
                                                                                                                                  +
                                                                                                                                  DHOR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                                                                                                                                  +
                                                                                                                                  DHOR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] )
                                                                                                                                  *
                                                                                                                                  this->JxW_cell[q];
                                                 }
                                          }
                                   }

                                   ///////////////////////////////// ///////////////////////////////
                                   // LOOP OVER QUADRATURE POINTS // // PROTON TRANSPORT EQUATION //
                                   ///////////////////////////////// ///////////////////////////////

                                   for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                   {
                                          /////////////////
                                          // LOOP OVER i //
                                          /////////////////

                                          for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                          {
                                                 /////////////////
                                                 // LOOP OVER j //
                                                 /////////////////

                                                 for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                                                 {
                                                        local_matrix(i,j) += - phi_protonic_electrical_potential[q][i]*(1.0)*( DHOR_current_density_Dhydrogen_concentration[q]*(1.0/molar_mass[indexH2-1])*phi_density[indexH2-1][q][j]
                                                                                                                               +
                                                                                                                               DHOR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                                                                                                                               +
                                                                                                                               DHOR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] )
                                                                                                                               *
                                                                                                                               this->JxW_cell[q];
                                                 }
                                          }
                                   }
                            }

                            /////////// ///////////// /////////
                            // ANODE // // CATHODE // // ORR //
                            /////////// ///////////// /////////

                            if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == ORR
                                ||
                                ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == ORR )
                            {
                                   ///////////////////////
                                   // LOOP OVER SPECIES //
                                   ///////////////////////

                                   for(unsigned int s = 0; s < n_species; ++s)
                                   {
                                          if( s+1 == indexO2 )
                                          {
                                                 ///////////////////////////////// //////////////////////////////
                                                 // LOOP OVER QUADRATURE POINTS // // MASS CONSERVATION for O2 //
                                                 ///////////////////////////////// //////////////////////////////

                                                 for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                                 {
                                                        /////////////////
                                                        // LOOP OVER i //
                                                        /////////////////

                                                        for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                                        {
                                                               /////////////////
                                                               // LOOP OVER j //
                                                               /////////////////

                                                               for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                                                               {
                                                                      local_matrix(i,j) += - phi_density[s][q][i]*multiplierO2*( DORR_current_density_Doxygen_concentration[q]*(1.0/molar_mass[s])*phi_density[s][q][j]
                                                                                                                                 +
                                                                                                                                 DORR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                                                                                                                                 +
                                                                                                                                 DORR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] )
                                                                                                                                 *
                                                                                                                                 this->JxW_cell[q];
                                                               }
                                                        }
                                                 }
                                          }

                                          if( s+1 == indexH2O )
                                          {
                                                 ///////////////////////////////// ///////////////////////////////
                                                 // LOOP OVER QUADRATURE POINTS // // MASS CONSERVATION for H2O //
                                                 ///////////////////////////////// ///////////////////////////////

                                                 for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                                 {
                                                        /////////////////
                                                        // LOOP OVER i //
                                                        /////////////////

                                                        for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                                        {
                                                               /////////////////
                                                               // LOOP OVER j //
                                                               /////////////////

                                                               for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                                                               {
                                                                      local_matrix(i,j) += - phi_density[s][q][i]*multiplierH2O*( DORR_current_density_Doxygen_concentration[q]*(1.0/molar_mass[indexO2-1])*phi_density[indexO2-1][q][j]
                                                                                                                                  +
                                                                                                                                  DORR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                                                                                                                                  +
                                                                                                                                  DORR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] )
                                                                                                                                  *
                                                                                                                                  this->JxW_cell[q];
                                                               }
                                                        }
                                                 }
                                          }
                                   }

                                   ///////////////////////////////// /////////////////////////////////
                                   // LOOP OVER QUADRATURE POINTS // // ELECTRON TRANSPORT EQUATION //
                                   ///////////////////////////////// /////////////////////////////////

                                   for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                   {
                                          /////////////////
                                          // LOOP OVER i //
                                          /////////////////

                                          for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                          {
                                                 /////////////////
                                                 // LOOP OVER j //
                                                 /////////////////

                                                 for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                                                 {
                                                        local_matrix(i,j) += - phi_electronic_electrical_potential[q][i]*(-1.0)*( DORR_current_density_Doxygen_concentration[q]*(1.0/molar_mass[indexO2-1])*phi_density[indexO2-1][q][j]
                                                                                                                                  +
                                                                                                                                  DORR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                                                                                                                                  +
                                                                                                                                  DORR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] )
                                                                                                                                  *
                                                                                                                                  this->JxW_cell[q];
                                                 }
                                          }
                                   }

                                   ///////////////////////////////// ///////////////////////////////
                                   // LOOP OVER QUADRATURE POINTS // // PROTON TRANSPORT EQUATION //
                                   ///////////////////////////////// ///////////////////////////////

                                   for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                   {
                                          /////////////////
                                          // LOOP OVER i //
                                          /////////////////

                                          for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                          {
                                                 /////////////////
                                                 // LOOP OVER j //
                                                 /////////////////

                                                 for(unsigned int j = 0; j < this->dofs_per_cell; ++j)
                                                 {
                                                        local_matrix(i,j) += - phi_protonic_electrical_potential[q][i]*(1.0)*( DORR_current_density_Doxygen_concentration[q]*(1.0/molar_mass[indexO2-1])*phi_density[indexO2-1][q][j]
                                                                                                                               +
                                                                                                                               DORR_current_density_Delectronic_electrical_potential[q]*phi_electronic_electrical_potential[q][j]
                                                                                                                               +
                                                                                                                               DORR_current_density_Dprotonic_electrical_potential[q]*phi_protonic_electrical_potential[q][j] )
                                                                                                                               *
                                                                                                                               this->JxW_cell[q];
                                                 }
                                          }
                                   }
                            }

                            local_matrix *= -1.0;

                            ////////////////////////
                            // DEALII -> APPFRAME //
                            ////////////////////////

                            this->dealII_to_appframe(cell_matrices,
                                                     local_matrix,
                                                     this->matrix_block_indices);
                     }
              }

              catch(const std::bad_cast& e)
              {
                     FcstUtilities::log << "Object is not of type FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
                     FcstUtilities::log << e.what() << std::endl;
              }
       }

       //////////////////
       // OTHER LAYERS //
       //////////////////

       else
       {
              FcstUtilities::log << "Layer you specified is not FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
              AssertThrow( false , ExcInternalError() );
       }
}

// ---                        ---
// --- assemble_cell_residual ---
// ---                        ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_residual,
                                                       const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                       FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
       ////////////////
       // TYPES INFO //
       ////////////////

       ///////////////////////////////////
       // ONLY ANODE OR CATHODE         //
       // CATALYST LAYERS CAN BE        //
       // ACCEPTED AS THE LAST ARGUMENT //
       // OF THIS FUNCTION              //
       ///////////////////////////////////

       const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
       const std::type_info& info          = layer->get_base_type();

       ////////////////////
       // CATALYST LAYER //
       ////////////////////

       if( info == CatalystLayer )
       {
              try
              {
                     FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);

                     ////////////////////////////////////////////////////////
                     //                                                    //
                     // WE FIRST CHECK IF gas_mixture IS A NULL PTR OR NOT //
                     //                                                    //
                     // - IF IT IS A NULL PTR THEN THE CORE OF             //
                     //   MASS/MOMENTUM/ENERGY TRANSPORT IS FORMED AROUND  //
                     //   FICKS EQUATIONS                                  //
                     //                                                    //
                     // - IF IT IS NOT A NULL PTR THEN THE CORE OF         //
                     //   MASS/MOMENTUM/ENERGY TRANSPORT IS FORMED AROUND  //
                     //   KERKHOF EQUATIONS                                //
                     //                                                    //
                     ////////////////////////////////////////////////////////

                     const FuelCellShop::Material::GasMixture* gas_mixture = ptr->get_gas_mixture();

                     /////////////////////////////
                     // FICKS + OTHER EQUATIONS //
                     /////////////////////////////

                     if( !gas_mixture )
                     {
                            if ( this->counter[0] )
                            {
                                this->make_assemblers_generic_constant_data();
                                this->counter[0] = false;
                            }

                            if ( this->counter[1] )
                            {
                                this->make_assemblers_cell_constant_data(cell_info);
                                this->counter[1] = false;
                            }

                            cell_residual_counter = true;
                            cell_matrix_counter = false;
                            this->make_assemblers_cell_variable_data(cell_info, layer);

                            for (unsigned int q=0; q < this->n_q_points_cell; ++q)
                            {
                                // ---- Proton Transport Equation ------------------------------
                                for (unsigned int i=0; i < (cell_info.fe(phi_m.fetype_index)).dofs_per_cell; ++i)
                                    cell_residual.block(phi_m.solution_index)(i) += ( this->JxW_cell[q] * factor_protontranseq_cell * current_cell[q] * phi_phiM_cell[q][i] );

                                // ---- Electron Transport Equation ------------------------------
                                for (unsigned int i=0; i < (cell_info.fe(phi_s.fetype_index)).dofs_per_cell; ++i)
                                    cell_residual.block(phi_s.solution_index)(i) += ( this->JxW_cell[q] * factor_electrontranseq_cell * current_cell[q] * phi_phiS_cell[q][i] );

                                // ---- Ficks Transport Equation - oxygen ------------------------
                                if ( x_oxygen.indices_exist )
                                    for (unsigned int i=0; i < (cell_info.fe(x_oxygen.fetype_index)).dofs_per_cell; ++i)
                                        cell_residual.block(x_oxygen.solution_index)(i) += ( this->JxW_cell[q] * factor_oxygentranseq_cell * current_cell[q] * phi_xOxygen_cell[q][i] );

                                // ---- Ficks Transport Equation - water ------------------------
                                if ( water_vapour_phase && x_water.indices_exist )
                                    for (unsigned int i=0; i < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++i)
                                        cell_residual.block(x_water.solution_index)(i) += ( this->JxW_cell[q] * factor_watertranseq_cell * current_cell[q] * phi_xWater_cell[q][i] );

                                // ---- Liquid Water Saturation Transport Equation---------------
                                if ( !water_vapour_phase && s_liquid_water.indices_exist )
                                    for (unsigned int i=0; i < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++i)
                                        cell_residual.block(s_liquid_water.solution_index)(i) += ( this->JxW_cell[q] * factor_saturationtranseq_cell * current_cell[q] * phi_s_cell[q][i] );

                                // ---- Thermal Transport Equation ------------------------------
                                if ( t_rev.indices_exist )
                                    for (unsigned int i=0; i < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++i)
                                        cell_residual.block(t_rev.solution_index)(i) += ( this->JxW_cell[q] * (-1.0) * heat_cell[q] * phi_T_cell[q][i] );
                            }
                     }

                     ///////////////////////////////
                     // KERKHOF + OTHER EQUATIONS //
                     ///////////////////////////////

                     else
                     {
                            ///////////////////////////
                            // GENERIC CONSTANT DATA //
                            ///////////////////////////

                            if( this->counter[2] )
                            {
                                   this->make_assemblers_generic_constant_data(cell_info,
                                                                               layer);
                                   this->counter[2] = false;
                            }

                            ////////////////////////
                            // CELL CONSTANT DATA //
                            ////////////////////////

                            if( this->counter[3] )
                            {
                                   this->make_assemblers_cell_constant_data2(cell_info);
                                   this->counter[3] = false;
                            }

                            ////////////////////////
                            // CELL VARIABLE DATA //
                            ////////////////////////
                            this->make_assemblers_cell_variable_data2(cell_info,
                                                                             layer);

                            /////////////////////////
                            // LOCAL CELL RESIDUAL //
                            /////////////////////////

                            Vector<double> local_residual(this->dofs_per_cell);

                            //////////////////////////////////////////////////
                            // SOURCE TERMS FOR MASS CONSERVATION EQUATIONS //
                            // OF OXYGEN, HYDROGEN, AND WATER VAPOR         //
                            // INSIDE ANODE OR/AND CATHODE CATALYST LAYERS  //
                            //////////////////////////////////////////////////

                            /////////// ///////////// /////////
                            // ANODE // // CATHODE // // HOR //
                            /////////// ///////////// /////////

                            if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == HOR
                                ||
                                ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == HOR )
                            {
                                   ///////////////////////
                                   // LOOP OVER SPECIES //
                                   ///////////////////////

                                   for(unsigned int s = 0; s < n_species; ++s)
                                   {
                                          if( s+1 == indexH2 )
                                          {
                                                 ///////////////////////////////// //////////////////////////////
                                                 // LOOP OVER QUADRATURE POINTS // // MASS CONSERVATION for H2 //
                                                 ///////////////////////////////// //////////////////////////////

                                                 for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                                 {
                                                        /////////////////
                                                        // LOOP OVER i //
                                                        /////////////////

                                                        for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                                        {
                                                               local_residual(i) += phi_density[s][q][i]*multiplierH2*HOR_current_density[q]*this->JxW_cell[q];
                                                        }
                                                 }
                                          }
                                   }

                                   ///////////////////////////////// /////////////////////////////////
                                   // LOOP OVER QUADRATURE POINTS // // ELECTRON TRANSPORT EQUATION //
                                   ///////////////////////////////// /////////////////////////////////

                                   for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                   {
                                          /////////////////
                                          // LOOP OVER i //
                                          /////////////////

                                          for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                          {
                                                 local_residual(i) += phi_electronic_electrical_potential[q][i]*(-1.0)*HOR_current_density[q]*this->JxW_cell[q];
                                          }
                                   }

                                   ///////////////////////////////// ///////////////////////////////
                                   // LOOP OVER QUADRATURE POINTS // // PROTON TRANSPORT EQUATION //
                                   ///////////////////////////////// ///////////////////////////////

                                   for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                   {
                                          /////////////////
                                          // LOOP OVER i //
                                          /////////////////

                                          for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                          {
                                                 local_residual(i) += phi_protonic_electrical_potential[q][i]*(1.0)*HOR_current_density[q]*this->JxW_cell[q];
                                          }
                                   }
                            }

                            /////////// ///////////// /////////
                            // ANODE // // CATHODE // // ORR //
                            /////////// ///////////// /////////

                            if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == ORR
                                ||
                                ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == ORR )
                            {
                                   ///////////////////////
                                   // LOOP OVER SPECIES //
                                   ///////////////////////

                                   for(unsigned int s = 0; s < n_species; ++s)
                                   {
                                          if( s+1 == indexO2 )
                                          {
                                                 ///////////////////////////////// //////////////////////////////
                                                 // LOOP OVER QUADRATURE POINTS // // MASS CONSERVATION for O2 //
                                                 ///////////////////////////////// //////////////////////////////

                                                 for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                                 {
                                                        /////////////////
                                                        // LOOP OVER i //
                                                        /////////////////

                                                        for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                                        {
                                                               local_residual(i) += phi_density[s][q][i]*multiplierO2*ORR_current_density[q]*this->JxW_cell[q];
                                                        }
                                                 }
                                          }

                                          if( s+1 == indexH2O )
                                          {
                                                 ///////////////////////////////// ///////////////////////////////
                                                 // LOOP OVER QUADRATURE POINTS // // MASS CONSERVATION for H2O //
                                                 ///////////////////////////////// ///////////////////////////////

                                                 for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                                 {
                                                        /////////////////
                                                        // LOOP OVER i //
                                                        /////////////////

                                                        for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                                        {
                                                               local_residual(i) += phi_density[s][q][i]*multiplierH2O*ORR_current_density[q]*this->JxW_cell[q];
                                                        }
                                                 }
                                          }
                                   }

                                   ///////////////////////////////// /////////////////////////////////
                                   // LOOP OVER QUADRATURE POINTS // // ELECTRON TRANSPORT EQUATION //
                                   ///////////////////////////////// /////////////////////////////////

                                   for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                   {
                                          /////////////////
                                          // LOOP OVER i //
                                          /////////////////

                                          for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                          {
                                                 local_residual(i) += phi_electronic_electrical_potential[q][i]*(-1.0)*ORR_current_density[q]*this->JxW_cell[q];
                                          }
                                   }

                                   ///////////////////////////////// ///////////////////////////////
                                   // LOOP OVER QUADRATURE POINTS // // PROTON TRANSPORT EQUATION //
                                   ///////////////////////////////// ///////////////////////////////

                                   for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                   {
                                          /////////////////
                                          // LOOP OVER i //
                                          /////////////////

                                          for(unsigned int i = 0; i < this->dofs_per_cell; ++i)
                                          {
                                                 local_residual(i) += phi_protonic_electrical_potential[q][i]*(1.0)*ORR_current_density[q]*this->JxW_cell[q];
                                          }
                                   }
                            }

                            ////////////////////////
                            // DEALII -> APPFRAME //
                            ////////////////////////

                            this->dealII_to_appframe(cell_residual,
                                                     local_residual,
                                                     this->residual_indices);
                     }
              }

              catch(const std::bad_cast& e)
              {
                     FcstUtilities::log << "Object is not of type FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
                     FcstUtilities::log << e.what() << std::endl;
              }
       }

       //////////////////
       // OTHER LAYERS //
       //////////////////

       else
       {
              FcstUtilities::log << "Layer you specified is not FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
              AssertThrow( false , ExcInternalError() );
       }
}

       /////////////////////////////////////////////////////
       /////////////////////////////////////////////////////
       // LOCAL CG FEM BASED ASSEMBLERS - make_ FUNCTIONS //
       /////////////////////////////////////////////////////
       /////////////////////////////////////////////////////

// ---                                       ---
// --- make_assemblers_generic_constant_data ---
// ---                                       ---

template<int dim>
template<typename INFO>
void
NAME::ReactionSourceTerms<dim>::make_assemblers_generic_constant_data(const INFO&                                InFo,
                                                                      FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
       this->dofs_per_cell   = InFo.get_fe_val_unsplit().dofs_per_cell;
       this->n_q_points_cell = InFo.get_fe_val_unsplit().n_quadrature_points;

       this->make_matrix_block_indices();
       this->make_residual_indices();

       ////////////////
       // TYPES INFO //
       ////////////////

       ///////////////////////////////////
       // ONLY ANODE OR CATHODE         //
       // CATALYST LAYERS CAN BE        //
       // ACCEPTED AS THE LAST ARGUMENT //
       // OF THIS FUNCTION              //
       ///////////////////////////////////

       const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
       const std::type_info& info          = layer->get_base_type();

       ////////////////////
       // CATALYST LAYER //
       ////////////////////

       if( info == CatalystLayer )
       {
              try
              {
                     FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);

                     T_mixture.resize(this->n_q_points_cell);
                     for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                            T_mixture[q] = ptr->get_gas_mixture()->get_temperature();
              }

              catch(const std::bad_cast& e)
              {
                     FcstUtilities::log << "Object is not of type FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
                     FcstUtilities::log << e.what() << std::endl;
              }
       }

       //////////////////
       // OTHER LAYERS //
       //////////////////

       else
       {
              FcstUtilities::log << "Layer you specified is not FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
              AssertThrow( false , ExcInternalError() );
       }
}

// ---                           ---
// --- make_matrix_block_indices ---
// ---                           ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::make_matrix_block_indices()
{
       ///////
       // 1 //
       ///////

       unsigned int index;

       //////////////////////////
       // MATRIX BLOCK INDICES //
       //////////////////////////

       ///////////////////////
       // KERKHOF EQUATIONS //
       ///////////////////////

       for(unsigned int g = 1; g <= n_species; ++g)
       {
              if( g == indexO2 )
              {
                     if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
                         ||
                         cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
                     {
                            ///////////////////////
                            // MASS CONSERVATION //
                            ///////////////////////

                            eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];

                            ///////////////////////
                            // DENSITY OF OXYGEN //
                            ///////////////////////

                            var_name = "density_" + var_postfixes[g-1];

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->matrix_block_index(eq_name,
                                                                                var_name);
                            this->matrix_block_indices.push_back(index);

                            ////////
                            // Fs //
                            ////////

                            var_name = "electronic_electrical_potential";

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->matrix_block_index(eq_name,
                                                                                var_name);
                            this->matrix_block_indices.push_back(index);

                            ////////
                            // Fm //
                            ////////

                            var_name = "protonic_electrical_potential";

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->matrix_block_index(eq_name,
                                                                                var_name);
                            this->matrix_block_indices.push_back(index);
                     }
              }

              if( g == indexH2 )
              {
                     if( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR
                         ||
                         cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )
                     {
                            ///////////////////////
                            // MASS CONSERVATION //
                            ///////////////////////

                            eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];

                            /////////////////////////
                            // DENSITY OF HYDROGEN //
                            /////////////////////////

                            var_name = "density_" + var_postfixes[g-1];

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->matrix_block_index(eq_name,
                                                                                var_name);
                            this->matrix_block_indices.push_back(index);

                            ////////
                            // Fs //
                            ////////

                            var_name = "electronic_electrical_potential";

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->matrix_block_index(eq_name,
                                                                                var_name);
                            this->matrix_block_indices.push_back(index);

                            ////////
                            // Fm //
                            ////////

                            var_name = "protonic_electrical_potential";

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->matrix_block_index(eq_name,
                                                                                var_name);
                            this->matrix_block_indices.push_back(index);
                     }
              }

              if( g == indexH2O )
              {
                     if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
                         ||
                         cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
                     {
                            ///////////////////////
                            // MASS CONSERVATION //
                            ///////////////////////

                            eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];

                            ///////////////////////
                            // DENSITY OF OXYGEN //
                            ///////////////////////

                            var_name = "density_" + var_postfixes[indexO2-1];

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->matrix_block_index(eq_name,
                                                                                var_name);
                            this->matrix_block_indices.push_back(index);

                            ////////////////////////////
                            // DENSITY OF WATER VAPOR //
                            ////////////////////////////

                            var_name = "density_" + var_postfixes[g-1];

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->matrix_block_index(eq_name,
                                                                                var_name);
                            this->matrix_block_indices.push_back(index);

                            ////////
                            // Fs //
                            ////////

                            var_name = "electronic_electrical_potential";

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->matrix_block_index(eq_name,
                                                                                var_name);
                            this->matrix_block_indices.push_back(index);

                            ////////
                            // Fm //
                            ////////

                            var_name = "protonic_electrical_potential";

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->matrix_block_index(eq_name,
                                                                                var_name);
                            this->matrix_block_indices.push_back(index);
                     }
              }
       }

       /////////////////////////////////
       // ELECTRON TRANSPORT EQUATION //
       /////////////////////////////////

       eq_name = "Electron Transport Equation";

       if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
           ||
           cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
       {
              var_name = "density_" + var_postfixes[indexO2-1];

              index = this->system_management->matrix_block_index(eq_name,
                                                                  var_name);
              this->matrix_block_indices.push_back(index);
       }

       if( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR
           ||
           cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )
       {
              var_name = "density_" + var_postfixes[indexH2-1];

              index = this->system_management->matrix_block_index(eq_name,
                                                                  var_name);
              this->matrix_block_indices.push_back(index);
       }

       var_name = "electronic_electrical_potential";

       index = this->system_management->matrix_block_index(eq_name,
                                                           var_name);
       this->matrix_block_indices.push_back(index);

       var_name = "protonic_electrical_potential";

       index = this->system_management->matrix_block_index(eq_name,
                                                           var_name);
       this->matrix_block_indices.push_back(index);

       ///////////////////////////////
       // PROTON TRANSPORT EQUATION //
       ///////////////////////////////

       eq_name = "Proton Transport Equation";

       if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
           ||
           cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
       {
              var_name = "density_" + var_postfixes[indexO2-1];

              index = this->system_management->matrix_block_index(eq_name,
                                                                  var_name);
              this->matrix_block_indices.push_back(index);
       }

       if( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR
           ||
           cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )
       {
              var_name = "density_" + var_postfixes[indexH2-1];

              index = this->system_management->matrix_block_index(eq_name,
                                                                  var_name);
              this->matrix_block_indices.push_back(index);
       }

       var_name = "electronic_electrical_potential";

       index = this->system_management->matrix_block_index(eq_name,
                                                           var_name);
       this->matrix_block_indices.push_back(index);

       var_name = "protonic_electrical_potential";

       index = this->system_management->matrix_block_index(eq_name,
                                                           var_name);
       this->matrix_block_indices.push_back(index);
}

// ---                       ---
// --- make_residual_indices ---
// ---                       ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::make_residual_indices()
{
       ///////
       // 1 //
       ///////

       unsigned int index;

       //////////////////////
       // RESIDUAL INDICES //
       //////////////////////

       ///////////////////////
       // KERKHOF EQUATIONS //
       ///////////////////////

       for(unsigned int g = 1; g <= n_species; ++g)
       {
              if( g == indexO2 )
              {
                     if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
                         ||
                         cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
                     {
                            ///////////////////////
                            // MASS CONSERVATION //
                            ///////////////////////

                            eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->equation_name_to_index(eq_name);
                            this->residual_indices.push_back(index);
                     }
              }

              if( g == indexH2 )
              {
                     if( anode_kinetics   && anode_kinetics->get_reaction_name()   == HOR
                         ||
                         cathode_kinetics && cathode_kinetics->get_reaction_name() == HOR )
                     {
                            ///////////////////////
                            // MASS CONSERVATION //
                            ///////////////////////

                            eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->equation_name_to_index(eq_name);
                            this->residual_indices.push_back(index);
                     }
              }

              if( g == indexH2O )
              {
                     if( anode_kinetics   && anode_kinetics->get_reaction_name()   == ORR
                         ||
                         cathode_kinetics && cathode_kinetics->get_reaction_name() == ORR )
                     {
                            ///////////////////////
                            // MASS CONSERVATION //
                            ///////////////////////

                            eq_name = eq_generic_prefix + "mass conservation - " + eq_postfixes[g-1];

                            ///////////
                            // index //
                            ///////////

                            index = this->system_management->equation_name_to_index(eq_name);
                            this->residual_indices.push_back(index);
                     }
              }
       }

       /////////////////////////////////
       // ELECTRON TRANSPORT EQUATION //
       /////////////////////////////////

       eq_name = "Electron Transport Equation";

       index = this->system_management->equation_name_to_index(eq_name);
       this->residual_indices.push_back(index);

       ///////////////////////////////
       // PROTON TRANSPORT EQUATION //
       ///////////////////////////////

       eq_name = "Proton Transport Equation";

       index = this->system_management->equation_name_to_index(eq_name);
       this->residual_indices.push_back(index);
}

// ---                                     ---
// --- make_assemblers_cell_constant_data2 ---
// ---                                     ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::make_assemblers_cell_constant_data2(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
       this->JxW_cell.resize(this->n_q_points_cell);

       ORR_current_density.resize(this->n_q_points_cell);
       DORR_current_density_Doxygen_concentration.resize(this->n_q_points_cell);
       DORR_current_density_Delectronic_electrical_potential.resize(this->n_q_points_cell);
       DORR_current_density_Dprotonic_electrical_potential.resize(this->n_q_points_cell);

       HOR_current_density.resize(this->n_q_points_cell);
       DHOR_current_density_Dhydrogen_concentration.resize(this->n_q_points_cell);
       DHOR_current_density_Delectronic_electrical_potential.resize(this->n_q_points_cell);
       DHOR_current_density_Dprotonic_electrical_potential.resize(this->n_q_points_cell);

       density_old.resize(this->n_q_points_cell);
       electronic_electrical_potential_old.resize(this->n_q_points_cell);
       protonic_electrical_potential_old.resize(this->n_q_points_cell);

       phi_density.resize( n_species,
                           std::vector< std::vector<double> >( this->n_q_points_cell,
                                                               std::vector<double>( this->dofs_per_cell )
                                                             )
                         );

       phi_electronic_electrical_potential.resize( this->n_q_points_cell,
                                                   std::vector<double>( this->dofs_per_cell ) );

       phi_protonic_electrical_potential.resize( this->n_q_points_cell,
                                                 std::vector<double>( this->dofs_per_cell ) );

       last_iter_cell = cell_info.global_data->find_vector("Newton iterate");
}

// ---                                     ---
// --- make_assemblers_cell_variable_data2 ---
// ---                                     ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::make_assemblers_cell_variable_data2(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                    FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
       ///////
       // 1 //
       ///////

       for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
       {
              this->JxW_cell[q] = cell_info.get_fe_val_unsplit().JxW(q);
       }

       ///////
       // 2 //
       ///////

       for(unsigned int g = 1; g <= n_species; ++g)
       {
              var_name = "density_" + var_postfixes[g-1];

              density_old[g-1] = cell_info.values[last_iter_cell][this->system_management->solution_name_to_index(var_name)];
       }

       ///////
       // 3 //
       ///////

       electronic_electrical_potential_old = cell_info.values[last_iter_cell][this->system_management->solution_name_to_index("electronic_electrical_potential")];
       protonic_electrical_potential_old   = cell_info.values[last_iter_cell][this->system_management->solution_name_to_index("protonic_electrical_potential")];

       ///////
       // 4 //
       ///////

       for(unsigned int g = 0; g < n_species; ++g)
       {
              for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
              {
                     for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
                     {
                            phi_density[g][q][k] = cell_info.get_fe_val_unsplit()[ density_extractors[g] ].value(k,q);
                     }
              }
       }

       ///////
       // 5 //
       ///////

       for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
       {
              for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
              {
                     phi_electronic_electrical_potential[q][k] = cell_info.get_fe_val_unsplit()[ electronic_electrical_potential_extractor ].value(k,q);
              }
       }

       ///////
       // 7 //
       ///////

       for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
       {
              for(unsigned int k = 0; k < this->dofs_per_cell; ++k)
              {
                     phi_protonic_electrical_potential[q][k] = cell_info.get_fe_val_unsplit()[ protonic_electrical_potential_extractor ].value(k,q);
              }
       }

       ////////////////
       // TYPES INFO //
       ////////////////

       ///////////////////////////////////
       // ONLY ANODE OR CATHODE         //
       // CATALYST LAYERS CAN BE        //
       // ACCEPTED AS THE LAST ARGUMENT //
       // OF THIS FUNCTION              //
       ///////////////////////////////////

       const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);
       const std::type_info& info          = layer->get_base_type();

       ////////////////////
       // CATALYST LAYER //
       ////////////////////

       if( info == CatalystLayer )
       {
              try
              {
                     FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);

                     /////////// ///////////// /////////
                     // ANODE // // CATHODE // // HOR //
                     /////////// ///////////// /////////

                     if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == HOR
                         ||
                         ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == HOR )
                     {
                            ///////
                            // 1 //
                            ///////

                            std::vector<double> hydrogen_concentration_old(this->n_q_points_cell);

                            if( indexH2 != -1 )
                            {
                                   for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                          hydrogen_concentration_old[q] = density_old[indexH2-1][q]/molar_mass[indexH2-1];
                            }
                            else
                            {
                                   AssertThrow( false, ExcInternalError() );
                            }

                            ///////
                            // 2 //
                            ///////

                            std::vector<FuelCellShop::SolutionVariable> solution_variables;

                            solution_variables.push_back( FuelCellShop::SolutionVariable( hydrogen_concentration_old,
                                                                                          VariableNames::hydrogen_concentration)
                                                        );

                            solution_variables.push_back( FuelCellShop::SolutionVariable( electronic_electrical_potential_old,
                                                                                          VariableNames::electronic_electrical_potential)
                                                        );

                            solution_variables.push_back( FuelCellShop::SolutionVariable( protonic_electrical_potential_old,
                                                                                          VariableNames::protonic_electrical_potential)
                                                        );

                            ///////
                            // 3 //
                            ///////

                            std::vector<VariableNames> derivative_flags;
                            derivative_flags.push_back(VariableNames::hydrogen_concentration);
                            derivative_flags.push_back(VariableNames::electronic_electrical_potential);
                            derivative_flags.push_back(VariableNames::protonic_electrical_potential);

                            ///////
                            // 4 //
                            ///////

                            ptr->set_solution(solution_variables);
                            ptr->set_cell_id(cell_info.dof_active_cell->index());
                            ptr->current_density( HOR_current_density );

                            ///////
                            // 5 //
                            ///////

                            std::map< VariableNames, std::vector<double> > Dcurrent;
                            ptr->set_derivative_flags(derivative_flags);
                            ptr->derivative_current_density(Dcurrent);

                            ///////
                            // 7 //
                            ///////

                            DHOR_current_density_Dhydrogen_concentration = Dcurrent.at(VariableNames::hydrogen_concentration);

                            for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                   DHOR_current_density_Dhydrogen_concentration[q] *= (Constants::R()*T_mixture[q])/(ptr->get_electrolyte()->get_H_H2()*1.0e-6);

                            DHOR_current_density_Delectronic_electrical_potential = Dcurrent.at(VariableNames::electronic_electrical_potential);
                            DHOR_current_density_Dprotonic_electrical_potential   = Dcurrent.at(VariableNames::protonic_electrical_potential);
                     }

                     /////////// ///////////// /////////
                     // ANODE // // CATHODE // // ORR //
                     /////////// ///////////// /////////

                     if( ptr->get_kinetics() == anode_kinetics   && ptr->get_kinetics()->get_reaction_name() == ORR
                         ||
                         ptr->get_kinetics() == cathode_kinetics && ptr->get_kinetics()->get_reaction_name() == ORR )
                     {
                            ///////
                            // 1 //
                            ///////

                            std::vector<double> oxygen_concentration_old(this->n_q_points_cell);

                            if( indexO2 != -1 )
                            {
                                   for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                          oxygen_concentration_old[q] = density_old[indexO2-1][q]/molar_mass[indexO2-1];
                            }
                            else
                            {
                                   AssertThrow( false, ExcInternalError() );
                            }

                            ///////
                            // 2 //
                            ///////

                            std::vector<FuelCellShop::SolutionVariable> solution_variables;

                            solution_variables.push_back( FuelCellShop::SolutionVariable( oxygen_concentration_old,
                                                                                          VariableNames::oxygen_concentration)
                                                        );

                            solution_variables.push_back( FuelCellShop::SolutionVariable( electronic_electrical_potential_old,
                                                                                          VariableNames::electronic_electrical_potential)
                                                        );

                            solution_variables.push_back( FuelCellShop::SolutionVariable( protonic_electrical_potential_old,
                                                                                          VariableNames::protonic_electrical_potential)
                                                        );

                            ///////
                            // 3 //
                            ///////

                            std::vector<VariableNames> derivative_flags;
                            derivative_flags.push_back(VariableNames::oxygen_concentration);
                            derivative_flags.push_back(VariableNames::electronic_electrical_potential);
                            derivative_flags.push_back(VariableNames::protonic_electrical_potential);

                            ///////
                            // 4 //
                            ///////

                            ptr->set_solution(solution_variables);
                            ptr->set_cell_id(cell_info.dof_active_cell->index());
                            ptr->current_density( ORR_current_density );

                            ///////
                            // 5 //
                            ///////

                            std::map< VariableNames, std::vector<double> > Dcurrent;
                            ptr->set_derivative_flags(derivative_flags);
                            ptr->derivative_current_density(Dcurrent);

                            ///////
                            // 7 //
                            ///////

                            DORR_current_density_Doxygen_concentration = Dcurrent.at(VariableNames::oxygen_concentration);

                            for(unsigned int q = 0; q < this->n_q_points_cell; ++q)
                                   DORR_current_density_Doxygen_concentration[q] *= (Constants::R()*T_mixture[q])/(ptr->get_electrolyte()->get_H_O2()*1.0e-6);

                            DORR_current_density_Delectronic_electrical_potential = Dcurrent.at(VariableNames::electronic_electrical_potential);
                            DORR_current_density_Dprotonic_electrical_potential   = Dcurrent.at(VariableNames::protonic_electrical_potential);
                     }
              }

              catch(const std::bad_cast& e)
              {
                     FcstUtilities::log << "Object is not of type FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
                     FcstUtilities::log << e.what() << std::endl;
              }
       }

       //////////////////
       // OTHER LAYERS //
       //////////////////

       else
       {
              FcstUtilities::log << "Layer you specified is not FuelCellShop::Layer::CatalystLayer<dim>" << std::endl;
              AssertThrow( false , ExcInternalError() );
       }
}

// ---                                       ---
// --- make_assemblers_generic_constant_data ---
// ---                                       ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::make_assemblers_generic_constant_data()
{
    FuelCellShop::Material::WaterVapor water;

    F = Constants::F();
    M_water = water.get_molar_mass() * 1.0e3;        // gm/mol

    //--------Block indices can't be filled here, as they depend on what equation we are in -------------------------------------------
    //--------While doing cell_matrix assembly, developers need to be wary of the fact that block_indices are still not filled yet-----

    //-----------Filling VariableInfo structures----------------------------------------------------------
    //-----------electronic_electrical_potential----------------------------------------------------------
    phi_s.solution_index = this->system_management->solution_name_to_index("electronic_electrical_potential");
    phi_s.fetype_index = this->system_management->block_info->base_element[phi_s.solution_index];
    phi_s.indices_exist = true;

    //----------protonic_electrical_potential-------------------------------------------------------------
    phi_m.solution_index = this->system_management->solution_name_to_index("protonic_electrical_potential");
    phi_m.fetype_index = this->system_management->block_info->base_element[phi_m.solution_index];
    phi_m.indices_exist = true;

    //-----------membrane water content----------------------------------------------------------
    if(this->system_management->solution_in_userlist("membrane_water_content")){
        lambda.solution_index = this->system_management->solution_name_to_index("membrane_water_content");
        lambda.fetype_index = this->system_management->block_info->base_element[lambda.solution_index];
        lambda.indices_exist = true;
    }


    if ( (cathode_reactionheat != NULL) || (anode_reactionheat != NULL) ) // It indirectly checks whether temperature_of_solid_phase is in user-defined list or not.
    {
        //-----------temperature_of_solid_phase-------------------------------------------------------
        t_rev.solution_index = this->system_management->solution_name_to_index("temperature_of_REV");
        t_rev.fetype_index = this->system_management->block_info->base_element[t_rev.solution_index];
        t_rev.indices_exist = true;
    }

    if ( cathode_kinetics != NULL ) // It indirectly checks whether oxygen_molar_fraction is in user-defined list or not.
    {
        //----------oxygen_molar_fraction-------------------------------------------------------------
        x_oxygen.solution_index = this->system_management->solution_name_to_index("oxygen_molar_fraction");
        x_oxygen.fetype_index = this->system_management->block_info->base_element[x_oxygen.solution_index];
        x_oxygen.indices_exist = true;
    }

    if ( anode_kinetics != NULL ) // It indirectly checks whether water_molar_fraction is in user-defined list or not.
    {
        //----------water_molar_fraction--------------------------------------------------------------
        x_water.solution_index = this->system_management->solution_name_to_index("water_molar_fraction");
        x_water.fetype_index = this->system_management->block_info->base_element[x_water.solution_index];
        x_water.indices_exist = true;
    }

    if ( water_vapour_phase == false ) // Liquid water saturation is needed only when water is produced in liquid phase during ORR
    {
        //----------liquid_water_saturation-------------------------------------------------------------
        s_liquid_water.solution_index = this->system_management->solution_name_to_index("liquid_water_saturation");
        s_liquid_water.fetype_index = this->system_management->block_info->base_element[s_liquid_water.solution_index];
        s_liquid_water.indices_exist = true;
    }
}

// ---                                    ---
// --- make_assemblers_cell_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    Assert( (phi_s.indices_exist && phi_m.indices_exist), ExcMessage("make_assemblers_generic_constant_data function not called before.") );

    this->n_q_points_cell = cell_info.get_fe_val_unsplit().n_quadrature_points;
    last_iter_cell = cell_info.global_data->find_vector("Newton iterate");

    //-------Allocation------------------------------------------------------------------------
    // ----- All containers intialized to zero by default -------------------------------------
    phi_phiS_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(phi_s.fetype_index)).dofs_per_cell ) );
    phi_phiM_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(phi_m.fetype_index)).dofs_per_cell ) );

    current_cell.resize( this->n_q_points_cell );

    if (t_rev.indices_exist)
    {
        phi_T_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(t_rev.fetype_index)).dofs_per_cell ) );
        heat_cell.resize( this->n_q_points_cell );
    }

    if (x_oxygen.indices_exist)
        phi_xOxygen_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(x_oxygen.fetype_index)).dofs_per_cell ) );

    if (x_water.indices_exist)
        phi_xWater_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(x_water.fetype_index)).dofs_per_cell ) );

    if (s_liquid_water.indices_exist)
        phi_s_cell.resize( this->n_q_points_cell, std::vector<double>( (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell ) );

    //-----------------------------------------------------------------
    this->JxW_cell.resize(this->n_q_points_cell);
}

// ---                                    ---
// --- make_assemblers_cell_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                   FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
    Assert( this->n_q_points_cell != 0, ExcMessage("make_assemblers_cell_constant_data function not called before.") );
    Assert( cell_matrix_counter != cell_residual_counter, ExcMessage("cell_matrix_counter should not match with cell_residual_counter, code has some bugs.") );

    // ----- type infos -------------
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);

    const std::type_info& base_layer = layer->get_base_type();

    // ----- dynamic cast and filling the containers -----------------
    try
    {
        if ( base_layer == CatalystLayer )
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);
            active_area_cell = ptr->get_active_area_Pt();       //--- Storing active area (cm^2/cm^3) for the catalyst layer ---------------

            std::vector<FuelCellShop::SolutionVariable> solution_variables;
            std::vector<VariableNames> deriv_flags;

            solution_variables.push_back( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][phi_s.solution_index], electronic_electrical_potential) );
            deriv_flags.push_back(electronic_electrical_potential);

            solution_variables.push_back( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][phi_m.solution_index], protonic_electrical_potential) );
            deriv_flags.push_back(protonic_electrical_potential);

            if(this->system_management->solution_in_userlist("membrane_water_content")){
                solution_variables.push_back( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][lambda.solution_index], membrane_water_content) );

                //No derivatives implemented, derivatives will be practically zero for all but extreme cases (i.e. conductivity scaling in agglomerateCL)
                //deriv_flags.push_back(membrane_water_content);

            }

            if (ptr->get_kinetics() == cathode_kinetics)
            {
                solution_variables.push_back( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][x_oxygen.solution_index], oxygen_molar_fraction) );
                deriv_flags.push_back(oxygen_molar_fraction);

                //----- Filling in factors for the source terms--------------------------
                factor_protontranseq_cell         = 1.0;
                factor_electrontranseq_cell       = (-1.0);
                factor_oxygentranseq_cell         = 1.0/(4.0*F);
                factor_watertranseq_cell          = ((-1.0)/(2.0*F)) * double(water_vapour_phase);
                factor_saturationtranseq_cell     = ((-1.0)/(2.0*F)) * M_water * double(!water_vapour_phase);

                if (t_rev.indices_exist)
                {
                    solution_variables.push_back( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );
                    deriv_flags.push_back(temperature_of_REV);

                    cathode_reactionheat->set_electrolyte_potential( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][phi_m.solution_index],
                                                                                                    protonic_electrical_potential) );
                    cathode_reactionheat->set_solid_potential( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][phi_s.solution_index],
                                                                                              electronic_electrical_potential) );
                    cathode_reactionheat->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index],
                                                                                          temperature_of_REV) );
                }
            }
            else if (ptr->get_kinetics() == anode_kinetics)
            {
                //---- Converting x_water into x_hydrogen ---------------------------
                std::vector<double> x_hydrogen( (cell_info.values[last_iter_cell][x_water.solution_index]).size() );
                for (unsigned int i=0; i < x_hydrogen.size(); ++i)
                    x_hydrogen[i] = 1.0 - (cell_info.values[last_iter_cell][x_water.solution_index][i]);

                solution_variables.push_back( FuelCellShop::SolutionVariable(x_hydrogen, hydrogen_molar_fraction) );
                deriv_flags.push_back(hydrogen_molar_fraction);

                //----- Filling in factors for the source terms--------------------------
                factor_protontranseq_cell        = (-1.0);
                factor_electrontranseq_cell      = 1.0;
                factor_oxygentranseq_cell        = 0.0;
                factor_watertranseq_cell         = 0.0;
                factor_saturationtranseq_cell    = 0.0;

                if (t_rev.indices_exist)
                {
                    solution_variables.push_back( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );
                    deriv_flags.push_back(temperature_of_REV);

                    anode_reactionheat->set_electrolyte_potential( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][phi_m.solution_index],
                                                                                                  protonic_electrical_potential) );
                    anode_reactionheat->set_solid_potential( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][phi_s.solution_index],
                                                                                            electronic_electrical_potential) );
                    anode_reactionheat->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index],
                                                                                        temperature_of_REV) );
                }
            }
            else
                AssertThrow(false, ExcMessage("Anode/Cathode kinetics are not set in the catalyst layer in the initialize method of the application.") );

            //------ Setting solution in the catalyst layer -----------------------------------------------
            ptr->set_solution(solution_variables);
            ptr->set_cell_id(cell_info.dof_active_cell->index());
            ptr->current_density( current_cell );

            //-------- Filling derivatives of source terms for cell_matrix assembly -----------------------
            if (cell_matrix_counter)
            {
                std::map< VariableNames, std::vector<double> > Dcurrent;
                ptr->set_derivative_flags(deriv_flags);
                ptr->derivative_current_density(Dcurrent);

                if (t_rev.indices_exist)
                {
                    dcurrent_dT_cell = Dcurrent.at(temperature_of_REV);

                    std::map< VariableNames, std::vector<double> > Dreactionheat;
                    if (ptr->get_kinetics() == cathode_kinetics)
                    {
                        cathode_reactionheat->set_derivative_flags(deriv_flags);
                        cathode_reactionheat->derivative_heat_source(Dreactionheat, Dcurrent, current_cell);

                        dheat_dxOxygen_cell = Dreactionheat.at(oxygen_molar_fraction);
                        dheat_dxWater_cell.assign(this->n_q_points_cell, 0.);
                    }
                    else if (ptr->get_kinetics() == anode_kinetics)
                    {
                        anode_reactionheat->set_derivative_flags(deriv_flags);
                        anode_reactionheat->derivative_heat_source(Dreactionheat, Dcurrent, current_cell);

                        dheat_dxOxygen_cell.assign(this->n_q_points_cell, 0.);
                        dheat_dxWater_cell.resize(this->n_q_points_cell, 0.);
                        for (unsigned int k=0; k < Dreactionheat.at(hydrogen_molar_fraction).size(); ++k)
                            dheat_dxWater_cell[k] = Dreactionheat.at(hydrogen_molar_fraction)[k] * (-1.0);      // Accounting for the fact that Dq/DxH2 = -Dq/DxH2O
                    }
                    dheat_dphiS_cell = Dreactionheat.at(electronic_electrical_potential);
                    dheat_dphiM_cell = Dreactionheat.at(protonic_electrical_potential);
                    dheat_dT_cell = Dreactionheat.at(temperature_of_REV);
                }

                dcurrent_dphiS_cell = Dcurrent.at(electronic_electrical_potential);
                dcurrent_dphiM_cell = Dcurrent.at(protonic_electrical_potential);

                if (ptr->get_kinetics() == cathode_kinetics)
                {
                    dcurrent_dxOxygen_cell = Dcurrent.at(oxygen_molar_fraction);
                    dcurrent_dxWater_cell.assign(this->n_q_points_cell, 0.);
                }
                else if (ptr->get_kinetics() == anode_kinetics)
                {
                    dcurrent_dxOxygen_cell.assign(this->n_q_points_cell, 0.);
                    dcurrent_dxWater_cell.resize(this->n_q_points_cell, 0.);
                    for (unsigned int k=0; k < Dcurrent.at(hydrogen_molar_fraction).size(); ++k)
                        dcurrent_dxWater_cell[k] = Dcurrent.at(hydrogen_molar_fraction)[k] * (-1.0);            // Accounting for the fact that Dj/DxH2 = -Dj/DxH2O
                }
            }

            else if (cell_residual_counter)
            {
                if ( t_rev.indices_exist )
                {
                    if ( ptr->get_kinetics() == cathode_kinetics )
                        cathode_reactionheat->heat_source(heat_cell, current_cell);

                    else if ( ptr->get_kinetics() == anode_kinetics )
                        anode_reactionheat->heat_source(heat_cell, current_cell);
                }
            }
        }

        else
            AssertThrow( false, ExcNotImplemented() );
    }

    catch(const std::bad_cast& e)
    {
        const std::type_info& info = typeid(*layer);
        FcstUtilities::log << "Object of type " << info.name() << " not implemented" << std::endl;
        FcstUtilities::log << e.what() << std::endl;
    }
        // --- continue ---

    //---------------------------------------------------------------------------------------------------------------
    //------------Looping over quadrature points in the cell --------------------------------------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        //-------JxW----------
        this->JxW_cell[q] = (cell_info.fe(phi_s.fetype_index)).JxW(q);

        //------ Filling shape functions etc ----------------------------------------------------------------------
        //------ This avoids recalculating shape functions etc for efficiency -------------------------------------

        for (unsigned int k=0; k < (cell_info.fe(phi_s.fetype_index)).dofs_per_cell; ++k)
            phi_phiS_cell[q][k] = (cell_info.fe(phi_s.fetype_index)).shape_value(k,q);

        for (unsigned int k=0; k < (cell_info.fe(phi_m.fetype_index)).dofs_per_cell; ++k)
            phi_phiM_cell[q][k] = (cell_info.fe(phi_m.fetype_index)).shape_value(k,q);

        //------- Checking based on boolean flags for other fe elements--------------------------------------------
        if ( x_oxygen.indices_exist )
            for (unsigned int k=0; k < (cell_info.fe(x_oxygen.fetype_index)).dofs_per_cell; ++k)
                phi_xOxygen_cell[q][k] = (cell_info.fe(x_oxygen.fetype_index)).shape_value(k,q);

        if ( x_water.indices_exist )
            for (unsigned int k=0; k < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++k)
                phi_xWater_cell[q][k] = (cell_info.fe(x_water.fetype_index)).shape_value(k,q);

        if ( t_rev.indices_exist )
            for (unsigned int k=0; k < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++k)
                phi_T_cell[q][k] = (cell_info.fe(t_rev.fetype_index)).shape_value(k,q);

        if ( s_liquid_water.indices_exist )
            for (unsigned int k=0; k < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++k)
                phi_s_cell[q][k] = (cell_info.fe(s_liquid_water.fetype_index)).shape_value(k,q);
    }
}

// ---                              ---
// --- assemble_matrix_for_equation ---
// ---                              ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::assemble_matrix_for_equation(FuelCell::ApplicationCore::MatrixVector& cell_matrices,
                                                            const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                            const std::string& eq_name,
                                                            const FEValuesBase<dim>& test_fe,
                                                            const std::vector< std::vector<double> >& test_shape_functions,
                                                            const double& sourceterm_factor)
{
    Assert( ((cell_matrix_counter == true) && (cell_residual_counter == false)), ExcInternalError() );
    Assert( this->n_q_points_cell != 0, ExcMessage("make_assemblers_cell_constant_data function not called before.") );
    Assert( ((eq_name == "Proton Transport Equation") || (eq_name == "Electron Transport Equation") || (eq_name == "Thermal Transport Equation") ||
            (eq_name == "Ficks Transport Equation - oxygen") || (eq_name == "Ficks Transport Equation - water") ||
            (eq_name == "Liquid Water Saturation Transport Equation")), ExcNotImplemented() );

    // --- Filling block indices --------------------------------------------
    phi_s.block_index = this->system_management->matrix_block_index(eq_name,"electronic_electrical_potential");
    phi_m.block_index = this->system_management->matrix_block_index(eq_name,"protonic_electrical_potential");
    if (x_oxygen.indices_exist)
        x_oxygen.block_index = this->system_management->matrix_block_index(eq_name,"oxygen_molar_fraction");
    if (x_water.indices_exist)
        x_water.block_index = this->system_management->matrix_block_index(eq_name,"water_molar_fraction");
    if (t_rev.indices_exist)
        t_rev.block_index = this->system_management->matrix_block_index(eq_name,"temperature_of_REV");

    std::vector<double> *dsource_dphiS, *dsource_dphiM, *dsource_dxOxygen, *dsource_dxWater, *dsource_dT;
    if ( eq_name == "Thermal Transport Equation" )
    {
        dsource_dphiS = &dheat_dphiS_cell;
        dsource_dphiM = &dheat_dphiM_cell;
        dsource_dxOxygen = &dheat_dxOxygen_cell;
        dsource_dxWater = &dheat_dxWater_cell;
        dsource_dT = &dheat_dT_cell;
    }
    else
    {
        dsource_dphiS = &dcurrent_dphiS_cell;
        dsource_dphiM = &dcurrent_dphiM_cell;
        dsource_dxOxygen = &dcurrent_dxOxygen_cell;
        dsource_dxWater = &dcurrent_dxWater_cell;
        dsource_dT = &dcurrent_dT_cell;
    }

    //-------- Looping over Quadrature points ----------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        //---------------LOOP over i -----------------------------------------------------------------
        for (unsigned int i=0; i < test_fe.dofs_per_cell; ++i)
        {
            //--------------LOOP(s) over j-------------------------------------------------------------
            for (unsigned int j=0; j < (cell_info.fe(phi_s.fetype_index)).dofs_per_cell; ++j)
                cell_matrices[phi_s.block_index].matrix(i,j) += ( this->JxW_cell[q] * test_shape_functions[q][i] * phi_phiS_cell[q][j] *
                sourceterm_factor * dsource_dphiS->at(q) );

            for (unsigned int j=0; j < (cell_info.fe(phi_m.fetype_index)).dofs_per_cell; ++j)
                cell_matrices[phi_m.block_index].matrix(i,j) += ( this->JxW_cell[q] * test_shape_functions[q][i] * phi_phiM_cell[q][j] *
                sourceterm_factor * dsource_dphiM->at(q) );

            if ( x_oxygen.indices_exist )
                for (unsigned int j=0; j < (cell_info.fe(x_oxygen.fetype_index)).dofs_per_cell; ++j)
                    cell_matrices[x_oxygen.block_index].matrix(i,j) += ( this->JxW_cell[q] * test_shape_functions[q][i] * phi_xOxygen_cell[q][j] *
                    sourceterm_factor * dsource_dxOxygen->at(q) );

            if ( x_water.indices_exist )
                for (unsigned int j=0; j < (cell_info.fe(x_water.fetype_index)).dofs_per_cell; ++j)
                    cell_matrices[x_water.block_index].matrix(i,j) += ( this->JxW_cell[q] * test_shape_functions[q][i] * phi_xWater_cell[q][j] *
                    sourceterm_factor * dsource_dxWater->at(q) );

            if ( t_rev.indices_exist )
                for (unsigned int j=0; j < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
                    cell_matrices[t_rev.block_index].matrix(i,j) += ( this->JxW_cell[q] * test_shape_functions[q][i] * phi_T_cell[q][j] *
                    sourceterm_factor * dsource_dT->at(q) );
        } // End Loop over "i"
    } // End Loop over "q"
}

// ---                     ---
// --- print_equation_info ---
// ---                     ---

template<int dim>
void
NAME::ReactionSourceTerms<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;

    FcstUtilities::log << "\"Reaction Source Terms \" PARAMETERS FOR \" Thermal Transport Equation \":" << std::endl;
    FcstUtilities::log << std::endl;

    FcstUtilities::log << "Irreversible heat source due to ORR:         " << irrev_heat_ccl << std::endl;
    FcstUtilities::log << "Irreversible heat source due to HOR:         " << irrev_heat_acl << std::endl;
    FcstUtilities::log << "Reversible heat source due to net reaction:  " << rev_heat << std::endl;
    FcstUtilities::log << "Reversible heat fraction in ORR:             " << factor_rev_heat_ccl << std::endl;
    FcstUtilities::log << "Water produced during ORR in vapour phase:   " << water_vapour_phase << std::endl;
    FcstUtilities::log << "Water vaporization heat sink in CCL:         " << water_vap_heat_ccl << std::endl;

    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
}

// ---                           ---
// ---  EXPLICIT INSTANTIATIONS  ---
// ---                           ---

template class NAME::ReactionSourceTerms<deal_II_dimension>;