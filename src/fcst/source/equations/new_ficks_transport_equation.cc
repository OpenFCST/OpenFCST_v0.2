//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: new_ficks_transport_equation.cc
//    - Description: Implementation of Ficks diffusion model.
//    - Developers: M. Secanell, Madhur Bhaiya, Valentin N. Zingan
//    - $Id: new_ficks_transport_equation.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include "new_ficks_transport_equation.h"

namespace NAME = FuelCellShop::Equation;

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::NewFicksTransportEquation<dim>::NewFicksTransportEquation(FuelCell::SystemManagement& system_management,
                                                                FuelCellShop::Material::PureGas* solute,
                                                                FuelCellShop::Material::PureGas* solvent)
:
EquationBase<dim>(system_management),
gas(solute),
solvent(solvent)
{
    std::stringstream s;
    s <<"Ficks Transport Equation - "<<gas->name_material();
    this->equation_name = s.str();

    std::stringstream ss;
    ss <<gas->name_material()<<"_molar_fraction";
    this->name_base_variable = ss.str();

    FcstUtilities::log << "->FuelCellShop::Equation::NewFicksTransportEquation" << std::endl;

    xi.indices_exist = false;
    t_rev.indices_exist = false;
    s_liquid_water.indices_exist = false;

    this->counter.resize(3, true);
}

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::NewFicksTransportEquation<dim>::NewFicksTransportEquation(FuelCell::SystemManagement& system_management, std::string& name_section)
:
EquationBase<dim>(system_management)
{
    
    this->equation_name = "Ficks Transport Equation - " + name_section;
    
    
    this->name_base_variable = name_section + "_molar_fraction";
    
    FcstUtilities::log << "->FuelCellShop::Equation::NewFicksTransportEquation" << std::endl;
    
    xi.indices_exist = false;
    t_rev.indices_exist = false;
    s_liquid_water.indices_exist = false;
    
    this->counter.resize(3, true);
}
// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::NewFicksTransportEquation<dim>::~NewFicksTransportEquation()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::declare_parameters(ParameterHandler& param)
{
    NAME::EquationBase<dim>::declare_parameters(param);
    
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {

            param.enter_subsection("Boundary data");
            {
                   param.declare_entry("species_flux",
                                       """",
                                        Patterns::Map(   Patterns::Integer(0) , Patterns::Double()   ),
                                       "Enter the molar flux in units of moles/(s*cm^2) that you would like to use as the Neumann boundary condition for each"
                                       "boundary_id. The format should be as follows: boundary_id1:value1, boundary_id2:value2. For example, using"
                                       "3:0.1, 41:0.25 will set the molar flux in boundary 3 to 0.1 and in boundary 41 to 0.25");
            }
            param.leave_subsection();

            param.enter_subsection("Boundary conditions");
            {
                param.declare_entry("Dirichlet Boundary Indicators",
                                    "",
                                    Patterns::List( Patterns::Integer(0) ),
                                    "Boundary_id(s) corresponding to Dirichlet boundaries (where molar fraction is to be specified), "
                                    "provided by comma-separated list of unsigned integers.");
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

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::initialize(ParameterHandler& param)
{
    NAME::EquationBase<dim>::initialize(param);
    
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {
            param.enter_subsection("Boundary data");
            {
                if( !param.get("species_flux").empty() )
                {
                    species_flux = FcstUtilities::string_to_map<unsigned int, double>( param.get("species_flux") );
                }
            }
            param.leave_subsection();
            
            param.enter_subsection("Boundary conditions");
            {
                dirichlet_bdry_ids = FcstUtilities::string_to_number<unsigned int>( Utilities::split_string_list( param.get("Dirichlet Boundary Indicators") ) );
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    this->make_internal_cell_couplings();
    this->make_boundary_types();
}

// ---                      ---
// --- assemble_cell_matrix ---
// ---                      ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::assemble_cell_matrix(FuelCell::ApplicationCore::MatrixVector& cell_matrices,
                                                           const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                           FuelCellShop::Layer::BaseLayer<dim>* const layer)
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

    this->make_assemblers_cell_variable_data(cell_info, layer);

    //-------- Looping over Quadrature points ----------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        //---------------LOOP over i -----------------------------------------------------------------
        for (unsigned int i=0; i < (cell_info.fe(xi.fetype_index)).dofs_per_cell; ++i)
        {
            //--------------LOOP(s) over j-------------------------------------------------------------

            //-----------Assembling Matrix for terms corresponding to "xi" BLOCK------------------------
            for (unsigned int j=0; j < (cell_info.fe(xi.fetype_index)).dofs_per_cell; ++j)
            {
                cell_matrices[xi.block_index].matrix(i,j) += grad_phi_xi_cell[q][i] * conc_Deff_cell[q] * grad_phi_xi_cell[q][j] * this->JxW_cell[q];
            }
            
            //----------TERM CORRESPONDING TO "T" BLOCK-----------------------------------------
            if (t_rev.indices_exist)
            {
                //-----Assembling Matrix------------------------------------------------------------
                for (unsigned int j=0; j < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++j)
                {
                    cell_matrices[t_rev.block_index].matrix(i,j) += grad_phi_xi_cell[q][i] * dconc_Deff_dT_cell[q] * cell_info.gradients[last_iter_cell][xi.solution_index][q] * phi_T_cell[q][j] *
                                                                    this->JxW_cell[q];
                }
            }

            //----------TERM CORRESPONDING TO "s" BLOCK-----------------------------------------
            if (s_liquid_water.indices_exist)
            {
                //-----Assembling Matrix------------------------------------------------------------
                for (unsigned int j=0; j < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++j)
                {
                    cell_matrices[s_liquid_water.block_index].matrix(i,j) += grad_phi_xi_cell[q][i] * dconc_Deff_ds_cell[q] * cell_info.gradients[last_iter_cell][xi.solution_index][q] * phi_s_cell[q][j] *
                                                                             this->JxW_cell[q];
                }
            }
        }
    }
}

// ---                        ---
// --- assemble_cell_residual ---
// ---                        ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::assemble_cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_rhs,
                                                             const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                             FuelCellShop::Layer::BaseLayer<dim>* const              layer)
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

    this->make_assemblers_cell_variable_data(cell_info, layer);

    for (unsigned int q=0; q < this->n_q_points_cell; ++q)
    {
        for (unsigned int i=0; i < (cell_info.fe(xi.fetype_index)).dofs_per_cell; ++i)
        {
            cell_rhs.block(xi.solution_index)(i) += grad_phi_xi_cell[q][i] * conc_Deff_cell[q] * cell_info.gradients[last_iter_cell][xi.solution_index][q] * this->JxW_cell[q];
        }
    }
}

// ---                      ---
// --- assemble_bdry_matrix ---
// ---                      ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::assemble_bdry_matrix(FuelCell::ApplicationCore::MatrixVector&                                 bdry_matrices,
                                                           const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                           FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
       // The boundary integral " q * [ ... ] " is always zero because of the boundary conditions we use:
       // - if xi = known then q = 0
       // - if [-C_tot * Di * grad_xi * n] = known then [ ... ] = 0
}

// ---                        ---
// --- assemble_bdry_residual ---
// ---                        ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::assemble_bdry_residual(FuelCell::ApplicationCore::FEVector&                                     bdry_residual,
                                                             const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                             FuelCellShop::Layer::BaseLayer<dim>* const              layer)
{
    if ( this->counter[0] )
    {
        this->make_assemblers_generic_constant_data();
        this->counter[0] = false;
    }

    if ( this->counter[2] )
    {
        this->make_assemblers_bdry_constant_data(bdry_info);
        this->counter[2] = false;
    }

    std::map<unsigned int, double>::const_iterator iter = species_flux.find( bdry_info.dof_face->boundary_indicator() );
    if (iter != species_flux.end() )
    {
        this->make_assemblers_bdry_variable_data(bdry_info, layer);

        //-------- Looping over Quadrature points ----------------------------
        for (unsigned int q = 0; q < this->n_q_points_bdry; ++q)
        {
            for (unsigned int i = 0; i < (bdry_info.fe(xi.fetype_index)).dofs_per_cell; ++i)
            {
                bdry_residual.block(xi.solution_index)(i) += phi_xi_bdry[q][i] * iter->second * this->JxW_bdry[q];
            }
        }
    }
}

// ---                     ---
// --- print_equation_info ---
// ---                     ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::print_equation_info() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "INTERNAL CELL COUPLINGS FOR \"" << this->equation_name << "\":" << std::endl;
    FcstUtilities::log << std::endl;

    couplings_map::const_iterator iter;

    for( iter = this->internal_cell_couplings.begin(); iter != this->internal_cell_couplings.end(); ++iter )
    {
        FcstUtilities::log << "\"" << iter->first << "\"" << ":";

        std::map<std::string, DoFTools::Coupling> int_map = iter->second;
        std::map<std::string, DoFTools::Coupling>::const_iterator int_iter;

        for( int_iter = int_map.begin(); int_iter != int_map.end(); ++int_iter )
            FcstUtilities::log << "\"" << int_iter->first << "\"" << " is coupled as " << int_iter->second << std::endl;

        FcstUtilities::log << std::endl;
    }

  FcstUtilities::log << "Initial data:";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;

  for(component_materialID_value_map::const_iterator iter  = this->component_materialID_value.begin();
                                                     iter != this->component_materialID_value.end();
                                                   ++iter)
  {
         std::map<types::material_id, double> tmp = iter->second;

         for(std::map<types::material_id, double>::const_iterator iter2  = tmp.begin();
                                                                  iter2 != tmp.end();
                                                                ++iter2)
         {
                FcstUtilities::log << "Name of the solution component  = " << iter->first   << std::endl;
                FcstUtilities::log << "Material id                     = " << iter2->first  << std::endl;
                FcstUtilities::log << "Value of the solution component = " << iter2->second << std::endl;
                FcstUtilities::log << std::endl;
         }
  }

  FcstUtilities::log << "Boundary data:";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;

  for(component_boundaryID_value_map::const_iterator iter  = this->component_boundaryID_value.begin();
                                                     iter != this->component_boundaryID_value.end();
                                                   ++iter)
  {
         std::map<types::boundary_id, double> tmp = iter->second;

         for(std::map<types::boundary_id, double>::const_iterator iter2  = tmp.begin();
                                                                  iter2 != tmp.end();
                                                                ++iter2)
         {
                FcstUtilities::log << "Name of the solution component  = " << iter->first   << std::endl;
                FcstUtilities::log << "Boundary id                     = " << iter2->first  << std::endl;
                FcstUtilities::log << "Value of the solution component = " << iter2->second << std::endl;
                FcstUtilities::log << std::endl;
         }
  }

    FcstUtilities::log << std::endl;
    FcstUtilities::log << "-------------------------------------------------------------------------------" << std::endl;
    FcstUtilities::log << std::endl;
}

// ---                              ---
// --- make_internal_cell_couplings ---
// ---                              ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::make_internal_cell_couplings()
{
    AssertThrow(this->system_management->solution_in_userlist(this->name_base_variable), VariableShouldExistForEquation(this->name_base_variable, this->equation_name) );
    AssertThrow(this->system_management->equation_name_to_index(this->equation_name) == this->system_management->solution_name_to_index(this->name_base_variable),
                IndexDoNotMatch(this->name_base_variable, this->equation_name) );

    std::map< std::string, DoFTools::Coupling > tmp;

    std::vector< std::string> sol_names = this->system_management->get_solution_names();

    for (unsigned int i = 0; i < sol_names.size(); ++i)
    {
        if (sol_names[i] == this->name_base_variable)
            tmp[this->name_base_variable] = DoFTools::always;

        else if (sol_names[i] == "temperature_of_REV")
            tmp["temperature_of_REV"] = DoFTools::always;

        else if (sol_names[i] == "liquid_water_saturation")
            tmp["liquid_water_saturation"] = DoFTools::always;

        else
            tmp[sol_names[i]] = DoFTools::none;
    }

    this->internal_cell_couplings[this->equation_name] = tmp;
}

// ---                     ---
// --- make_boundary_types ---
// ---                     ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::make_boundary_types()
{
    for (unsigned int index = 1; index <= dirichlet_bdry_ids.size(); ++index)
    {
        BoundaryType temp_dirich;

        std::ostringstream streamOut;
        streamOut << index;
        std::string temp_name = "Dirichlet_" + streamOut.str();

        temp_dirich.boundary_name = temp_name;
        temp_dirich.boundary_id = dirichlet_bdry_ids[index-1];
        temp_dirich.boundary_condition = "Dirichlet";

        this->boundary_types.push_back(temp_dirich);
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
void
NAME::NewFicksTransportEquation<dim>::make_assemblers_generic_constant_data()
{
    //-----------Filling VariableInfo structures------------------------------------------
    xi.solution_index = this->system_management->solution_name_to_index(this->name_base_variable);
    xi.block_index = this->system_management->matrix_block_index(this->equation_name, this->name_base_variable);
    xi.fetype_index = this->system_management->block_info->base_element[xi.solution_index];
    xi.indices_exist = true;

    //----------temperature_of_solid_phase----------------------------------------------------
    if ( this->system_management->solution_in_userlist("temperature_of_REV") )
    {
        t_rev.solution_index = this->system_management->solution_name_to_index("temperature_of_REV");
        t_rev.block_index = this->system_management->matrix_block_index(this->equation_name, "temperature_of_REV");
        t_rev.fetype_index = this->system_management->block_info->base_element[t_rev.solution_index];
        t_rev.indices_exist = true;
    }

    //----------liquid_water_saturation----------------------------------------------------
    if ( this->system_management->solution_in_userlist("liquid_water_saturation") )
    {
        s_liquid_water.solution_index = this->system_management->solution_name_to_index("liquid_water_saturation");
        s_liquid_water.block_index = this->system_management->matrix_block_index(this->equation_name, "liquid_water_saturation");
        s_liquid_water.fetype_index = this->system_management->block_info->base_element[s_liquid_water.solution_index];
        s_liquid_water.indices_exist = true;
    }
}

// ---                                    ---
// --- make_assemblers_cell_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::make_assemblers_cell_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    Assert( xi.indices_exist, ExcMessage("make_assemblers_generic_constant_data function not called before.") );

    this->n_q_points_cell = (cell_info.fe(xi.fetype_index)).n_quadrature_points;
    this->last_iter_cell = cell_info.global_data->find_vector("Newton iterate");

    conc_Deff_cell.resize(this->n_q_points_cell);
    dconc_Deff_dT_cell.resize(this->n_q_points_cell);
    dconc_Deff_ds_cell.resize(this->n_q_points_cell);

    //-------------Allocation------------------------------------------
    grad_phi_xi_cell.resize( this->n_q_points_cell, std::vector< Tensor<1,dim> >( (cell_info.fe(xi.fetype_index)).dofs_per_cell ) );

    if (t_rev.indices_exist)
        phi_T_cell.resize( this->n_q_points_cell, std::vector< double >( (cell_info.fe(t_rev.fetype_index)).dofs_per_cell ) );

    if (s_liquid_water.indices_exist)
        phi_s_cell.resize( this->n_q_points_cell, std::vector< double >( (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell ) );

    this->JxW_cell.resize(this->n_q_points_cell);

    //-----------------------------------------------------------------
}

// ---                                    ---
// --- make_assemblers_bdry_constant_data ---
// ---                                    ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::make_assemblers_bdry_constant_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info)
{
    Assert( xi.indices_exist, ExcMessage("make_assemblers_generic_constant_data function not called before.") );

    this->n_q_points_bdry = (bdry_info.fe(xi.fetype_index)).n_quadrature_points;
    last_iter_bdry = bdry_info.global_data->find_vector("Newton iterate");

    //-------------Allocation------------------------------------------
    phi_xi_bdry.resize( this->n_q_points_bdry, std::vector<double>( (bdry_info.fe(xi.fetype_index)).dofs_per_cell ) );

    //-----------------------------------------------------------------
    this->JxW_bdry.resize(this->n_q_points_bdry);
}

// ---                                    ---
// --- make_assemblers_cell_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::make_assemblers_cell_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info,
                                                                         FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    Assert( this->n_q_points_cell != 0, ExcMessage("make_assemblers_cell_constant_data function not called before.") );

    //---------------Effective Transport Properties---------------------------------------------------------------
    // ----- type infos -------------
    const std::type_info& GasDiffusionLayer = typeid(FuelCellShop::Layer::GasDiffusionLayer<dim>);
    const std::type_info& MicroPorousLayer  = typeid(FuelCellShop::Layer::MicroPorousLayer<dim>);
    const std::type_info& CatalystLayer = typeid(FuelCellShop::Layer::CatalystLayer<dim>);

    const std::type_info& base_layer = layer->get_base_type();

    // Creating some internal variables -- to be used commonly by all the dynamically casted layers
    Table< 2, Tensor<2,dim> > Deff_iso;
    double p, T;
    int index_gas, index_solvent;
    //------
    std::vector< Tensor<2,dim> > Deff_noniso;
    std::map< VariableNames, std::vector< Tensor<2,dim> > > dDeff_du;
    std::vector<VariableNames> deriv_flags;

    // ----- dynamic cast and filling the containers -----------------
    try
    {
        if (base_layer == GasDiffusionLayer)
        {
            FuelCellShop::Layer::GasDiffusionLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::GasDiffusionLayer<dim>* >(layer);

            if (!t_rev.indices_exist)         // t_rev.indices_exist = false
            {
                ptr->get_T_and_p(T,p);
                double concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);       // mol/cm^3
                ptr->effective_gas_diffusivity(Deff_iso);               // m^2/s
                ptr->get_gas_index(this->gas, index_gas);
                ptr->get_gas_index(this->solvent, index_solvent);
                for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                    conc_Deff_cell[q] = concentration * Deff_iso(index_gas,index_solvent)*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                         // mol/(cm-s)
            }
            else                                // t_rev.indices_exist = true
            {
                ptr->get_p(p);
                ptr->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );
                deriv_flags.push_back(temperature_of_REV);
                if (s_liquid_water.indices_exist)
                {
                    ptr->set_saturation( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][s_liquid_water.solution_index], liquid_water_saturation) );
                    deriv_flags.push_back(liquid_water_saturation);
                }

                ptr->compute_gas_diffusion(gas, solvent);
                ptr->effective_gas_diffusivity(Deff_noniso);            // m^2/s

                ptr->set_derivative_flags(deriv_flags);
                ptr->derivative_effective_gas_diffusivity(dDeff_du);
                for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                {
                    T = cell_info.values[last_iter_cell][t_rev.solution_index][q];
                    double concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/cm^3
                    double dC_dT = (-1.)*((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/(cm^3-K)

                    conc_Deff_cell[q] = concentration * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                                            // mol/(cm-s)
                    dconc_Deff_dT_cell[q] = dC_dT * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2) +
                                            concentration * dDeff_du[temperature_of_REV][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                       // mol/(cm-s-K)

                    if (s_liquid_water.indices_exist)
                        dconc_Deff_ds_cell[q] = concentration * dDeff_du[liquid_water_saturation][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);              // mol/(cm-s)
                }
            }
        }

        else if (base_layer == MicroPorousLayer)
        {
            FuelCellShop::Layer::MicroPorousLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::MicroPorousLayer<dim>* >(layer);

            if (!t_rev.indices_exist)         // t_rev.indices_exist = false
            {
                ptr->get_T_and_p(T,p);
                double concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);       // mol/cm^3
                ptr->effective_gas_diffusivity(Deff_iso);               // m^2/s
                ptr->get_gas_index(this->gas, index_gas);
                ptr->get_gas_index(this->solvent, index_solvent);
                for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                    conc_Deff_cell[q] = concentration * Deff_iso(index_gas,index_solvent)*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                         // mol/(cm-s)
            }
            else                                // t_rev.indices_exist = true
            {
                ptr->get_p(p);
                ptr->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );
                deriv_flags.push_back(temperature_of_REV);
                if (s_liquid_water.indices_exist)
                {
                    ptr->set_saturation( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][s_liquid_water.solution_index], liquid_water_saturation) );
                    deriv_flags.push_back(liquid_water_saturation);
                }

                ptr->compute_gas_diffusion(gas, solvent);
                ptr->effective_gas_diffusivity(Deff_noniso);            // m^2/s

                ptr->set_derivative_flags(deriv_flags);
                ptr->derivative_effective_gas_diffusivity(dDeff_du);
                for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                {
                    T = cell_info.values[last_iter_cell][t_rev.solution_index][q];
                    double concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/cm^3
                    double dC_dT = (-1.)*((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/(cm^3-K)

                    conc_Deff_cell[q] = concentration * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                                            // mol/(cm-s)
                    dconc_Deff_dT_cell[q] = dC_dT * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2) +
                                            concentration * dDeff_du[temperature_of_REV][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                       // mol/(cm-s-K)

                    if (s_liquid_water.indices_exist)
                        dconc_Deff_ds_cell[q] = concentration * dDeff_du[liquid_water_saturation][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);              // mol/(cm-s)
                }
            }
        }

        else if (base_layer == CatalystLayer)
        {
            FuelCellShop::Layer::CatalystLayer<dim>* ptr = dynamic_cast< FuelCellShop::Layer::CatalystLayer<dim>* >(layer);

            if (!t_rev.indices_exist)         // t_rev.indices_exist = false
            {
                ptr->get_T_and_p(T,p);
                double concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T))*Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3);       // mol/cm^3
                ptr->effective_gas_diffusivity(Deff_iso);               // m^2/s
                ptr->get_gas_index(this->gas, index_gas);
                ptr->get_gas_index(this->solvent, index_solvent);
                for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                    conc_Deff_cell[q] = concentration * Deff_iso(index_gas,index_solvent)*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                         // mol/(cm-s)
            }
            else                                // t_rev.indices_exist = true
            {
                ptr->get_p(p);
                ptr->set_temperature( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][t_rev.solution_index], temperature_of_REV) );
                deriv_flags.push_back(temperature_of_REV);
                if (s_liquid_water.indices_exist)
                {
                    ptr->set_saturation( FuelCellShop::SolutionVariable(&cell_info.values[last_iter_cell][s_liquid_water.solution_index], liquid_water_saturation) );
                    deriv_flags.push_back(liquid_water_saturation);
                }

                ptr->compute_gas_diffusion(gas, solvent);
                ptr->effective_gas_diffusivity(Deff_noniso);            // m^2/s

                ptr->set_derivative_flags(deriv_flags);
                ptr->derivative_effective_gas_diffusivity(dDeff_du);
                for (unsigned int q=0; q<this->n_q_points_cell; ++q)
                {
                    T = cell_info.values[last_iter_cell][t_rev.solution_index][q];
                    double concentration = ((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/cm^3
                    double dC_dT = (-1.)*((p*Units::convert(1.,Units::ATM_to_PA))/(Constants::R()*T*T)) * Units::convert(1.,Units::PER_C_UNIT3, Units::PER_UNIT3); // mol/(cm^3-K)

                    conc_Deff_cell[q] = concentration * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                                            // mol/(cm-s)
                    dconc_Deff_dT_cell[q] = dC_dT * Deff_noniso[q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2) +
                                            concentration * dDeff_du[temperature_of_REV][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);                       // mol/(cm-s-K)

                    if (s_liquid_water.indices_exist)
                        dconc_Deff_ds_cell[q] = concentration * dDeff_du[liquid_water_saturation][q]*Units::convert(1.,Units::C_UNIT2, Units::UNIT2);              // mol/(cm-s)
                }
            }
        }
        else
            AssertThrow( false, ExcNotImplemented() );
    }
    catch(const std::bad_cast& e)
    {
        const std::type_info& info = typeid(*layer);
        FcstUtilities::log << "Object of type "<<info.name()<<" not implemented"<< std::endl;
        FcstUtilities::log << e.what() << std::endl;
    }

    //---------------------------------------------------------------------------------------------------------------
    //------------Looping over quadrature points in the cell --------------------------------------------------------
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
        //-------JxW----------
        this->JxW_cell[q] = (cell_info.fe(xi.fetype_index)).JxW(q);

        //------ Filling shape functions etc ----------------------------------------------------------------------
        //------ This avoids recalculating shape functions etc for efficiency -------------------------------------
        for (unsigned int k=0; k < (cell_info.fe(xi.fetype_index)).dofs_per_cell; ++k)
        {
            grad_phi_xi_cell[q][k] = (cell_info.fe(xi.fetype_index)).shape_grad(k,q);
        }

        if (t_rev.indices_exist)
            for (unsigned int k=0; k < (cell_info.fe(t_rev.fetype_index)).dofs_per_cell; ++k)
                phi_T_cell[q][k] = (cell_info.fe(t_rev.fetype_index)).shape_value(k,q);

        if (s_liquid_water.indices_exist)
            for (unsigned int k=0; k < (cell_info.fe(s_liquid_water.fetype_index)).dofs_per_cell; ++k)
                phi_s_cell[q][k] = (cell_info.fe(s_liquid_water.fetype_index)).shape_value(k,q);
    }
}

// ---                                    ---
// --- make_assemblers_bdry_variable_data ---
// ---                                    ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::make_assemblers_bdry_variable_data(const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                                                         FuelCellShop::Layer::BaseLayer<dim>* const layer)
{
    Assert( this->n_q_points_bdry != 0, ExcMessage("make_assemblers_bdry_constant_data function not called before.") );

    //---------------------------------------------------------------------------------------------------------------
    //------------Looping over quadrature points in the cell --------------------------------------------------------
    for (unsigned int q=0; q < this->n_q_points_bdry; ++q)
    {
        this->JxW_bdry[q] = (bdry_info.fe(xi.fetype_index)).JxW(q);

        for (unsigned int k=0; k < (bdry_info.fe(xi.fetype_index)).dofs_per_cell; ++k)
            phi_xi_bdry[q][k] = (bdry_info.fe(xi.fetype_index)).shape_value(k,q);
    }
}

// ---            ---
// --- class_test ---
// ---            ---

template<int dim>
void
NAME::NewFicksTransportEquation<dim>::class_test()
{
  FuelCellShop::Material::Hydrogen H2;
  FuelCellShop::Material::Oxygen O2;
  FuelCell::ApplicationCore::BlockInfo block_info;
  Table< 2, DoFTools::Coupling > cell_couplings;
  Table< 2, DoFTools::Coupling > flux_couplings;
  FuelCell::SystemManagement sys(block_info, cell_couplings, flux_couplings);
  FcstUtilities::log<<"Create object NewFicksTransportEquation:"<<std::endl;
  FuelCellShop::Equation::NewFicksTransportEquation<dim> test(sys,&O2,&H2);
  test.print_equation_info();
}

// ---                         ---
// --- explicit instantiations ---
// ---                         ---

template class NAME::NewFicksTransportEquation<deal_II_dimension>;
