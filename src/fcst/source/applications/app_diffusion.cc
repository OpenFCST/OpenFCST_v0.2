// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2014 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: app_diffusion.cc
// - Description: This class describes diffusion in fuel cell cathodes
//                Ficks, one gas
// - Developers: Mayank Sabharwal,University of Alberta
//               Marc Secanell, University of Alberta
// - Id: $Id: app_diffusion.cc 2616 2014-08-15 22:57:14Z secanell $
//
// ----------------------------------------------------------------------------

#include "app_diffusion.h"

namespace NAME  = FuelCell::Application;

// ---              ---
// ---              ---
// --- AppDiffusion ---
// ---              ---
// ---              ---

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- constructor ---
// ---             ---
std::string name_section = "oxygen";
template<int dim>
NAME::AppDiffusion<dim>::AppDiffusion( boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data )
:
FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>(data),
ficks_transport_equation(this->system_management,name_section)
{
  this->repair_diagonal = true;
  FcstUtilities::log <<  "->FuelCell::Application::AppDiffusion-" << dim << "D" << std::endl;

}

// ---            ---
// --- destructor ---
// ---            ---

template<int dim>
NAME::AppDiffusion<dim>::~AppDiffusion()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::AppDiffusion<dim>::declare_parameters(ParameterHandler& param)
{

    OptimizationBlockMatrixApplication<dim>::declare_parameters(param);

    // Declare parameters in system management:
    this->system_management.declare_parameters(param);

    // Declare parameters in operating conditions:
    OC.declare_parameters(param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Materials");
        {   param.declare_entry("solute",
                                "oxygen",
                                Patterns::Selection("oxygen | nitrogen | hydrogen | water | air | helium"),
                                "Gas species to be treated as the solute for solving the Fick's Law ");
            param.declare_entry("solvent",
                                "nitrogen",
                                Patterns::Selection("oxygen | nitrogen | hydrogen | water | air | helium"),
                                "Gas species to be treated as the solvent for solving the Fick's Law ");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
   
    // Declare layer class:
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Cathode gas diffusion layer", param);
    
    // Declare equation class:
    ficks_transport_equation.declare_parameters(param);
    
    
}

// ---            ---
// --- initialize ---
// ---            ---
template<int dim>
void
NAME::AppDiffusion<dim>::initialize(ParameterHandler& param)
{
    OptimizationBlockMatrixApplication<dim>::initialize(param);
    
    std::string solutename;
    std::string solventname;
    
    // Initialize parameters in system management:FuelCellShop::Material::PureGas
    this->system_management.initialize(param);
    
    // Initialize parameters in operating conditions:
    OC.initialize(param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Materials");
        {   
            solutename = param.get("solute");
            solventname = param.get("solvent");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    if (solutename=="oxygen")
        solute = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Oxygen);
    else if (solutename=="nitrogen")
        solute = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Nitrogen);
    else if (solutename=="helium")
        solute = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Helium);
    else if (solutename=="hydrogen")
        solute = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Hydrogen);
    else if (solutename=="water")
        solute = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::WaterVapor);
    else
    {
        FcstUtilities::log << "Solute specified does not match any of components defined in the database" << std::endl;
        AssertThrow(false, ExcInternalError());
    }

    if (solventname=="oxygen")
        solvent = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Oxygen);
    else if (solventname=="nitrogen")
        solvent = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Nitrogen);
    else if (solventname=="helium")
        solvent = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Helium);
    else if (solventname=="hydrogen")
        solvent = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::Hydrogen);
    else if (solventname=="water")
        solvent = boost::shared_ptr<FuelCellShop::Material::PureGas>(new FuelCellShop::Material::WaterVapor);
    else
    {
        FcstUtilities::log << "Solvent specified does not match any of components defined in the database" << std::endl;
        AssertThrow(false, ExcInternalError());
    }
    // Initialize materials and layers:
    solute->declare_parameters(param);
    solvent->declare_parameters(param);
    solute->initialize(param);
    solvent->initialize(param);
    
    // Initialize gases and material classes:  
    std::vector< FuelCellShop::Material::PureGas* > gases;
    gases.push_back(solute.get());
    gases.push_back(solvent.get());
        
    // Initialize layer classes:
    CGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Cathode gas diffusion layer", param);
    CGDL->set_gases_and_compute(gases, OC.get_pc_atm(), OC.get_T());
    
    // Initialize parameters for physics classes:
    ficks_transport_equation.set_solute_and_solvent (solute.get(), solvent.get(), param);
    ficks_transport_equation.initialize(param);
    
    // --- we make cell couplings for this problem ---
    std::vector<couplings_map> tmp;
    tmp.push_back( ficks_transport_equation.get_internal_cell_couplings()    );
    this->system_management.make_cell_couplings(tmp);
    
    // Now, initialize object that are used to setup initial solution and boundary conditions:    
    this->component_materialID_value_maps.push_back( ficks_transport_equation.get_component_materialID_value()    );
    OC.adjust_initial_solution(this->component_materialID_value_maps, this->mesh_generator);
    
    this->component_boundaryID_value_maps.push_back( ficks_transport_equation.get_component_boundaryID_value() );
    OC.adjust_boundary_conditions(this->component_boundaryID_value_maps, this->mesh_generator);
    
    // --- and then allocate memory for matrices ---
    this->remesh_matrices();
    
    // Output options:
    // - system info:
    //this->system_management.print_system_info();
    
    // - layers info:
    CGDL->print_layer_properties();
        
    // - equations info:
    ficks_transport_equation.print_equation_info();
    
}

// ---               ---
// --- init_solution ---
// ---               ---

template<int dim>
void
NAME::AppDiffusion<dim>::initialize_solution(FuelCell::ApplicationCore::FEVector& initial_guess,
                                             std::shared_ptr<Function<dim> > initial_function)
{
    DoFApplication<dim>::initialize_solution(initial_guess); 
}

       ///////////////////////////////////
       ///////////////////////////////////
       // LOCAL CG FEM BASED ASSEMBLERS //
       ///////////////////////////////////
       ///////////////////////////////////

// ---             ---
// --- cell_matrix ---
// ---             ---

template<int dim>
void
NAME::AppDiffusion<dim>::cell_matrix(FuelCell::ApplicationCore::MatrixVector&                                 cell_matrices,
                                     const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{

    if(      CGDL->belongs_to_material(cell_info.cell->material_id())   )
    {
        ficks_transport_equation.assemble_cell_matrix(cell_matrices, cell_info, CGDL.get());
        

    }
    
    else
    {
      
      std::stringstream ss;
      ss <<solute->name_material()<<"_molar_fraction";
      std::string speciesname = ss.str();
      const int index = this->system_management.matrix_block_index(ficks_transport_equation.get_equation_name(), speciesname);
      cell_matrices[index].matrix.all_zero();
      
      
    }

}

// ---               ---
// --- cell_residual ---
// ---               ---

template<int dim>
void
NAME::AppDiffusion<dim>::cell_residual(FuelCell::ApplicationCore::FEVector&                                     cell_res,
                                       const typename FuelCell::ApplicationCore::DoFApplication<dim>::CellInfo& cell_info)
{
    if(      CGDL->belongs_to_material(cell_info.cell->material_id())   )
    {
        cell_res = 0;//ficks_transport_equation.assemble_cell_residual(cell_res, cell_info, CGDL.get());
    
    }
    else
    {
        cell_res = 0;
        
    }
}


       /////////////////////
       /////////////////////
       // OTHER FUNCTIONS //
       /////////////////////
       /////////////////////

// ---              ---
// --- dirichlet_bc ---
// ---              ---

template<int dim>
void
NAME::AppDiffusion<dim>::dirichlet_bc(std::map<unsigned int, double>& boundary_values) const
{
  FuelCell::InitialAndBoundaryData::make_constant_DirichletBC_values( boundary_values,
                                                                      *this->mapping,
                                                                      *this->dof,
                                                                      this->system_management,
                                                                      this->component_boundaryID_value_maps );
} 

/////////////////////
/////////////////////
// POST-PROCESSING //
/////////////////////
/////////////////////


/////////////////////
/////////////////////
// POST-PROCESSING //
/////////////////////
/////////////////////

// ---                ---
// --- bdry_responses ---
// ---                ---


template<int dim>
void
NAME::AppDiffusion<dim>::bdry_responses(std::vector<double>&                                                     dst,
                                      const typename FuelCell::ApplicationCore::DoFApplication<dim>::FaceInfo& bdry_info,
                                      const FuelCell::ApplicationCore::FEVector& src)
{
    //User input bdry id which will be later added to the class parameters when flux computation class is created
    int user_input_bdry = 2;
    int index_gas, index_solvent;
    CGDL->get_gas_index(this->solute.get(), index_gas);
    CGDL->get_gas_index(this->solvent.get(), index_solvent);
    // -- Find out what material is the cell made of, i.e. MEA layer)
    const unsigned int bdry_id = bdry_info.dof_face->boundary_indicator();
    
    unsigned solution_index = bdry_info.global_data->find_vector("Solution");
    Table< 2, Tensor<2,dim> > Deff_iso;
    
    FuelCellShop::Equation::VariableInfo xi;
    
    if ( this->system_management.solution_in_userlist("oxygen_molar_fraction") )
    {
        xi.solution_index = this->system_management.solution_name_to_index("oxygen_molar_fraction"); 
        xi.fetype_index = this->system_management.block_info->base_element[xi.solution_index];
        xi.indices_exist = true;
    }
     
    CGDL->effective_gas_diffusivity(Deff_iso);
    
    int n_q_points_bdry = (bdry_info.fe(xi.fetype_index)).n_quadrature_points;
    
    //-------- Looping over Quadrature points ----------------------------
    std::vector<double> JxW_bdry;
    JxW_bdry.resize(n_q_points_bdry);
    
    if (bdry_id==user_input_bdry)
    {
        for (unsigned int q = 0; q < n_q_points_bdry; ++q)
        {
            JxW_bdry[q] = (bdry_info.fe(xi.fetype_index)).JxW(q);
            
            dst[0]+= -Deff_iso(index_gas,index_solvent)*bdry_info.gradients[solution_index][xi.solution_index][q] * bdry_info.fe(xi.fetype_index).normal_vector(q)* JxW_bdry[q];
	    
        }   
    }

    
}

// ---          ---
// --- evaluate ---
// ---          ---
// This function 
template <int dim>
double
NAME::AppDiffusion<dim>::evaluate (const FuelCell::ApplicationCore::FEVectors& src)
{
    return 0.0;
}

// ---          ---
// --- data_out ---
// ---          ---

template<int dim>
void
NAME::AppDiffusion<dim>::data_out(const std::string& filename,
                                  const FuelCell::ApplicationCore::FEVectors& src)
{
  //////////////
  // SOLUTION //
  //////////////

  // --- Find solution ---
  FuelCell::ApplicationCore::FEVector solution = src.vector( src.find_vector("Solution") );

  // --- Assign solution names ---
  std::vector<std::string> solution_names;

  solution_names.push_back("oxygen_molar_fraction");

  // --- Assign solution interpretations ---
  this->solution_interpretations.clear();
  this->solution_interpretations.resize(this->element->n_blocks(),
                                        DataComponentInterpretation::component_is_scalar);

  ///////////////////////////////////
  // Do further POST-PROCESSING    //
  ///////////////////////////////////

  // --- Create vector of PostProcessing objects ---
  std::vector< DataPostprocessor<dim>* > PostProcessing;

  // --- output ---
  DoFApplication<dim>::data_out( filename,
                                 solution,
                                 solution_names);
}

/////////////////////////////
/////////////////////////////
// EXPLICIT INSTANTIATIONS //
/////////////////////////////
/////////////////////////////

// ---              ---
// --- AppDiffusion ---
// ---              ---
template class NAME::AppDiffusion<deal_II_dimension>;
