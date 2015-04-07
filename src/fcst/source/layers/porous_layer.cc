//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: porous_layer.cc
//    - Description: Child of base layer that implements common functionality for handling gases.
//    - Developers: M. Secanell, Madhur Bhaiya and Valentin Zingan
//    - $Id: porous_layer.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include "layers/porous_layer.h"

namespace NAME = FuelCellShop::Layer;

template <int dim>
void
NAME::PorousLayer<dim>::declare_parameters (const std::string &name, ParameterHandler &param) const
{
  FuelCellShop::Layer::BaseLayer<dim>::declare_parameters(name,param);

  param.enter_subsection("Fuel cell data");
  {
      param.enter_subsection(name);
      {
          param.enter_subsection("Generic data");
          {
              param.declare_entry("Porosity",
                                  "1.0",
                                  Patterns::Double(),
                                  "Volume fraction of void space to total volume in the layer, [-]");           
              param.declare_entry("Knudsen pore radius [um]",
                                  "1.0e6",
                                  Patterns::Double(),
                                  "Pore radius used to compute Knudsen diffusivity");
              
              param.declare_entry("Permeability_XX",
                                  "1.0",
                                  Patterns::Double(),
                                  "XX component of the absolute permeability tensor of the layer, Units [cm^2]");
              param.declare_entry("Permeability_XY",
                                  "0.0",
                                  Patterns::Double(),
                                  "XY component of the absolute permeability tensor of the layer, Units [cm^2]");
              param.declare_entry("Permeability_XZ",
                                  "0.0",
                                  Patterns::Double(),
                                  "XZ component of the absolute permeability tensor of the layer, Units [cm^2]");
              
              param.declare_entry("Permeability_YX",
                                  "0.0",
                                  Patterns::Double(),
                                  "YX component of the absolute permeability tensor of the layer, Units [cm^2]");
              param.declare_entry("Permeability_YY",
                                  "1.0",
                                  Patterns::Double(),
                                  "YY component of the absolute permeability tensor of the layer, Units [cm^2]");
              param.declare_entry("Permeability_YZ",
                                  "0.0",
                                  Patterns::Double(),
                                  "YZ component of the absolute permeability tensor of the layer, Units [cm^2]");
              
              param.declare_entry("Permeability_ZX",
                                  "0.0",
                                  Patterns::Double(),
                                  "ZX component of the absolute permeability tensor of the layer, Units [cm^2]");
              param.declare_entry("Permeability_ZY",
                                  "0.0",
                                  Patterns::Double(),
                                  "ZY component of the absolute permeability tensor of the layer, Units [cm^2]");
              param.declare_entry("Permeability_ZZ",
                                  "1.0",
                                  Patterns::Double(),
                                  "ZZ component of the absolute permeability tensor of the layer, Units [cm^2]");
              
              param.declare_entry("Forchheimer_permeability_XX",
                                  "1.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Forchheimer_permeability_XY",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Forchheimer_permeability_XZ",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");                
              param.declare_entry("Forchheimer_permeability_YX",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Forchheimer_permeability_YY",
                                  "1.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Forchheimer_permeability_YZ",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");
              
              param.declare_entry("Forchheimer_permeability_ZX",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Forchheimer_permeability_ZY",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Forchheimer_permeability_ZZ",
                                  "1.0",
                                  Patterns::Double(),
                                  " ");
              
              param.declare_entry("Tortuosity_XX",
                                  "1.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Tortuosity_XY",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Tortuosity_XZ",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");
              
              param.declare_entry("Tortuosity_YX",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Tortuosity_YY",
                                  "1.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Tortuosity_YZ",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");
              
              param.declare_entry("Tortuosity_ZX",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Tortuosity_ZY",
                                  "0.0",
                                  Patterns::Double(),
                                  " ");
              param.declare_entry("Tortuosity_ZZ",
                                  "1.0",
                                  Patterns::Double(),
                                  " ");
          }
          param.leave_subsection();
      }
      param.leave_subsection();
  }
  param.leave_subsection();
  
  
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::PorousLayer<dim>::initialize (ParameterHandler &param)
{
  NAME::BaseLayer<dim>::initialize(param);

  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection(this->name);
    {        
        param.enter_subsection("Generic data");
        {
            porosity = param.get_double("Porosity");   
            
            Knudsen_radius = param.get_double("Knudsen pore radius [um]");
            
            SymmetricTensor<2,3> permeability_tmp;
            
            permeability_tmp[0][0] = param.get_double("Permeability_XX");
            permeability_tmp[0][1] = param.get_double("Permeability_XY");
            permeability_tmp[0][2] = param.get_double("Permeability_XZ");
            
            permeability_tmp[1][0] = param.get_double("Permeability_YX");
            permeability_tmp[1][1] = param.get_double("Permeability_YY");
            permeability_tmp[1][2] = param.get_double("Permeability_YZ");
            
            permeability_tmp[2][0] = param.get_double("Permeability_ZX");
            permeability_tmp[2][1] = param.get_double("Permeability_ZY");
            permeability_tmp[2][2] = param.get_double("Permeability_ZZ");
            
            for(unsigned int i = 0; i < dim; ++i)
                for(unsigned int j = 0; j < dim; ++j)
                    permeability[i][j] = permeability_tmp[i][j];
                
                for(unsigned int i = 0; i < dim; ++i)
                    for(unsigned int j = 0; j < dim; ++j)
                        SQRT_permeability[i][j] = std::sqrt( permeability[i][j] );
                    
                    FullMatrix<double> m(dim,dim);
                
                for(unsigned int i = 0; i < dim; ++i)
                    for(unsigned int j = 0; j < dim; ++j)
                        m(i,j) = permeability[i][j];
                    
                    m.gauss_jordan();
                
                for(unsigned int i = 0; i < dim; ++i)
                    for(unsigned int j = 0; j < dim; ++j)
                        permeability_INV[i][j] = m(i,j);
                    
                    for(unsigned int i = 0; i < dim; ++i)
                        for(unsigned int j = 0; j < dim; ++j)
                            m(i,j) = SQRT_permeability[i][j];
                        
                        m.gauss_jordan();
                    
                    for(unsigned int i = 0; i < dim; ++i)
                        for(unsigned int j = 0; j < dim; ++j)
                            SQRT_permeability_INV[i][j] = m(i,j);
                        
                        SymmetricTensor<2,3> Forchheimer_permeability_tmp;
                    
                    Forchheimer_permeability_tmp[0][0] = param.get_double("Forchheimer_permeability_XX");
                    Forchheimer_permeability_tmp[0][1] = param.get_double("Forchheimer_permeability_XY");
                    Forchheimer_permeability_tmp[0][2] = param.get_double("Forchheimer_permeability_XZ");
                    
                    Forchheimer_permeability_tmp[1][0] = param.get_double("Forchheimer_permeability_YX");
                    Forchheimer_permeability_tmp[1][1] = param.get_double("Forchheimer_permeability_YY");
                    Forchheimer_permeability_tmp[1][2] = param.get_double("Forchheimer_permeability_YZ");
                    
                    Forchheimer_permeability_tmp[2][0] = param.get_double("Forchheimer_permeability_ZX");
                    Forchheimer_permeability_tmp[2][1] = param.get_double("Forchheimer_permeability_ZY");
                    Forchheimer_permeability_tmp[2][2] = param.get_double("Forchheimer_permeability_ZZ");
                    
                    for(unsigned int i = 0; i < dim; ++i)
                        for(unsigned int j = 0; j < dim; ++j)
                            Forchheimer_permeability[i][j] = Forchheimer_permeability_tmp[i][j];
                        
                        SymmetricTensor<2,3> tortuosity_tmp;
                    
                    tortuosity_tmp[0][0] = param.get_double("Tortuosity_XX");
                    tortuosity_tmp[0][1] = param.get_double("Tortuosity_XY");
                    tortuosity_tmp[0][2] = param.get_double("Tortuosity_XZ");
                    
                    tortuosity_tmp[1][0] = param.get_double("Tortuosity_YX");
                    tortuosity_tmp[1][1] = param.get_double("Tortuosity_YY");
                    tortuosity_tmp[1][2] = param.get_double("Tortuosity_YZ");
                    
                    tortuosity_tmp[2][0] = param.get_double("Tortuosity_ZX");
                    tortuosity_tmp[2][1] = param.get_double("Tortuosity_ZY");
                    tortuosity_tmp[2][2] = param.get_double("Tortuosity_ZZ");
                    
                    for(unsigned int i = 0; i < dim; ++i)
                        for(unsigned int j = 0; j < dim; ++j)
                            tortuosity[i][j] = tortuosity_tmp[i][j];
        }
        param.leave_subsection();
    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

// ---              ---
// --- get_porosity ---
// ---              ---

template<int dim>
void
NAME::PorousLayer<dim>::get_porosity(std::vector<double>&             dst,
                                                 const std::vector< Point<dim> >& points) const
{
  AssertThrow( !porosity_is_constant , ExcInternalError() );
  AssertThrow( dst.size() == points.size() , ExcDimensionMismatch(dst.size(), points.size()) );

  print_caller_name(__FUNCTION__);
}

// ---                  ---
// --- get_permeability ---
// ---                  ---

template<int dim>
void
NAME::PorousLayer<dim>::get_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const
{
  AssertThrow( permeability_is_constant , ExcInternalError() );

  for(unsigned int q = 0; q < dst.size(); ++q)
    dst[q] = permeability;
}

// ---                  ---
// --- get_permeability ---
// ---                  ---

template<int dim>
void
NAME::PorousLayer<dim>::get_permeability(std::vector< SymmetricTensor<2,dim> >& dst,
                                                     const std::vector< Point<dim> >&       points) const
{
  AssertThrow( !permeability_is_constant , ExcInternalError() );
  AssertThrow( dst.size() == points.size() , ExcDimensionMismatch(dst.size(), points.size()) );

  print_caller_name(__FUNCTION__);
}

// ---                       ---
// --- get_SQRT_permeability ---
// ---                       ---

template<int dim>
void
NAME::PorousLayer<dim>::get_SQRT_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const
{
  AssertThrow( permeability_is_constant , ExcInternalError() );

  for(unsigned int q = 0; q < dst.size(); ++q)
    dst[q] = SQRT_permeability;
}

// ---                       ---
// --- get_SQRT_permeability ---
// ---                       ---

template<int dim>
void
NAME::PorousLayer<dim>::get_SQRT_permeability(std::vector< SymmetricTensor<2,dim> >& dst,
                                                          const std::vector< Point<dim> >&       points) const
{
  AssertThrow( !permeability_is_constant , ExcInternalError() );
  AssertThrow( dst.size() == points.size() , ExcDimensionMismatch(dst.size(), points.size()) );

  std::vector< SymmetricTensor<2,dim> > tmp( points.size() );
  this->get_permeability(tmp,
                         points);

  for(unsigned int q = 0; q < dst.size(); ++q)
    for(unsigned int i = 0; i < dim; ++i)
      for(unsigned int j = 0; j < dim; ++j)
        dst[q][i][j] = std::sqrt( tmp[q][i][j] );
}

// ---                      ---
// --- get_permeability_INV ---
// ---                      ---

template<int dim>
void
NAME::PorousLayer<dim>::get_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst) const
{
  AssertThrow( permeability_is_constant , ExcInternalError() );

  for(unsigned int q = 0; q < dst.size(); ++q)
    dst[q] = permeability_INV;
}

// ---                      ---
// --- get_permeability_INV ---
// ---                      ---

template<int dim>
void
NAME::PorousLayer<dim>::get_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst,
                                                         const std::vector< Point<dim> >&       points) const
{
  AssertThrow( !permeability_is_constant , ExcInternalError() );
  AssertThrow( dst.size() == points.size() , ExcDimensionMismatch(dst.size(), points.size()) );

  std::vector< SymmetricTensor<2,dim> > tmp( points.size() );
  this->get_permeability(tmp,
                         points);

  for(unsigned int q = 0; q < dst.size(); ++q)
  {
    FullMatrix<double> m(dim,dim);

    for(unsigned int i = 0; i < dim; ++i)
      for(unsigned int j = 0; j < dim; ++j)
        m(i,j) = tmp[q][i][j];

    m.gauss_jordan();

    for(unsigned int i = 0; i < dim; ++i)
      for(unsigned int j = 0; j < dim; ++j)
        dst[q][i][j] = m(i,j);
  }
}

// ---                           ---
// --- get_SQRT_permeability_INV ---
// ---                           ---

template<int dim>
void
NAME::PorousLayer<dim>::get_SQRT_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst) const
{
  AssertThrow( permeability_is_constant , ExcInternalError() );

  for(unsigned int q = 0; q < dst.size(); ++q)
    dst[q] = SQRT_permeability_INV;
}

// ---                           ---
// --- get_SQRT_permeability_INV ---
// ---                           ---

template<int dim>
void
NAME::PorousLayer<dim>::get_SQRT_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst,
                                                              const std::vector< Point<dim> >&       points) const
{
  AssertThrow( !permeability_is_constant , ExcInternalError() );
  AssertThrow( dst.size() == points.size() , ExcDimensionMismatch(dst.size(), points.size()) );

  std::vector< SymmetricTensor<2,dim> > tmp( points.size() );
  this->get_permeability(tmp,
                         points);

  for(unsigned int q = 0; q < dst.size(); ++q)
  {
    FullMatrix<double> m(dim,dim);

    for(unsigned int i = 0; i < dim; ++i)
      for(unsigned int j = 0; j < dim; ++j)
        m(i,j) = std::sqrt( tmp[q][i][j] );

    m.gauss_jordan();

    for(unsigned int i = 0; i < dim; ++i)
      for(unsigned int j = 0; j < dim; ++j)
        dst[q][i][j] = m(i,j);
  }
}

// ---                              ---
// --- get_Forchheimer_permeability ---
// ---                              ---

template<int dim>
void
NAME::PorousLayer<dim>::get_Forchheimer_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const
{
  AssertThrow( permeability_is_constant , ExcInternalError() );

  for(unsigned int q = 0; q < dst.size(); ++q)
    dst[q] = Forchheimer_permeability;
}

// ---                              ---
// --- get_Forchheimer_permeability ---
// ---                              ---

template<int dim>
void
NAME::PorousLayer<dim>::get_Forchheimer_permeability(std::vector< SymmetricTensor<2,dim> >& dst,
                                                     const std::vector< Point<dim> >&       points) const
{
  AssertThrow( !permeability_is_constant , ExcInternalError() );
  AssertThrow( dst.size() == points.size() , ExcDimensionMismatch(dst.size(), points.size()) );

  print_caller_name(__FUNCTION__);
}

// ---                ---
// --- get_tortuosity ---
// ---                ---

template<int dim>
void
NAME::PorousLayer<dim>::get_tortuosity(std::vector< SymmetricTensor<2,dim> >& dst) const
{
  AssertThrow( tortuosity_is_constant , ExcInternalError() );

  for(unsigned int q = 0; q < dst.size(); ++q)
    dst[q] = tortuosity;
}

// ---                ---
// --- get_tortuosity ---
// ---                ---

template<int dim>
void
NAME::PorousLayer<dim>::get_tortuosity(std::vector< SymmetricTensor<2,dim> >& dst,
                                       const std::vector< Point<dim> >&       points) const
{
  AssertThrow( !tortuosity_is_constant , ExcInternalError() );
  AssertThrow( dst.size() == points.size() , ExcDimensionMismatch(dst.size(), points.size()) );

  print_caller_name(__FUNCTION__);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::PorousLayer<dim>::set_gases_and_compute (std::vector<FuelCellShop::Material::PureGas*>& gases_in,
                                               const double& pressure_in,
                                               const double& temperature_in)
{
    Assert(gases_in.size() >= 2, ExcMessage("Number of gases should be more than or equal to two in PorousLayer::set_gases_and_compute method."));

    this->gases = gases_in;
    this->pressure = pressure_in;
    this->temperature = temperature_in;

    this->D_ECtheory.reinit(gases_in.size(), gases_in.size());

    this->dD_ECtheory_dx.resize(this->derivative_flags.size() );
    for (unsigned int d = 0; d< this->derivative_flags.size(); ++d)
        this->dD_ECtheory_dx[d].reinit(gases_in.size(), gases_in.size());

    FuelCellShop::Material::GasMixture mixture("noname");
    mixture.set_gases(gases_in);   

    this->D_ECtheory = mixture.get_ChapmanEnskog_diffusion_coefficients(101325.0*pressure_in, temperature_in);
    
    for(unsigned int d = 0; d < this->derivative_flags.size(); ++d)
    {
        if( this->derivative_flags[d] == total_pressure )
            this->dD_ECtheory_dx[d] = mixture.get_DChapmanEnskog_diffusion_coefficients_Dpressure(101325.0*pressure_in, temperature_in);
        else if( this->derivative_flags[d] == temperature_of_REV )
            this->dD_ECtheory_dx[d] = mixture.get_DChapmanEnskog_diffusion_coefficients_Dtemperature(101325.0*pressure_in, temperature_in);
        else
            AssertThrow(false, ExcNotImplemented());
    }
    
    this->gas_mixture = nullptr;
};

//---------------------------------------------------------------------------
template <int dim>
void
NAME::PorousLayer<dim>::compute_gas_diffusion(FuelCellShop::Material::PureGas* solute_gas,
                                              FuelCellShop::Material::PureGas* solvent_gas)
{
    Assert( this->pressure != 0., ExcMessage("Pressure is not set inside the porous layer. Use set_gases method in initialization of the application.") );
    Assert( this->T_vector.is_initialized(), ExcMessage("Temperature is not set inside the porous layer. Use set_temperature method before calling this function.") );
    Assert( std::find(this->gases.begin(), this->gases.end(), solute_gas)!=this->gases.end(), ExcMessage("Solute gas not set using set_gases method in initialization of the application.") );
    Assert( std::find(this->gases.begin(), this->gases.end(), solvent_gas)!=this->gases.end(), ExcMessage("Solvent gas not set using set_gases method in initialization of the application.") );

    std::vector< FuelCellShop::Material::PureGas* > gases_in;
    gases_in.push_back(solute_gas);
    gases_in.push_back(solvent_gas);

    FuelCellShop::Material::GasMixture mixture("noname");
    mixture.set_gases(gases_in);

    std::vector<double> press;
    press.resize(this->T_vector.size());

    for(unsigned int q = 0; q < press.size(); ++q)
           press[q] = 101325.0*this->pressure;

    std::vector<double> temp;
    temp.resize(this->T_vector.size());

    for(unsigned int q = 0; q < temp.size(); ++q)
           temp[q] = this->T_vector[q];

    this->D_bulk.clear();
    this->D_bulk.resize(this->T_vector.size());
    this->dD_bulk_dT.clear();
    this->dD_bulk_dT.resize(this->T_vector.size());

    mixture.get_ChapmanEnskog_diffusion_coefficient(press,
                                                    temp,
                                                    this->D_bulk);

    mixture.get_DChapmanEnskog_diffusion_coefficient_Dtemperature(press,
                                                                  temp,
                                                                  this->dD_bulk_dT);
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::PorousLayer<dim>::get_gas_index(FuelCellShop::Material::PureGas* gas_type,
                                      int& index) const
{
    Assert( std::find(this->gases.begin(), this->gases.end(), gas_type)!=this->gases.end(), ExcMessage("Gas not set using set_gases_and_compute method in initialization of the application.") );

    for (unsigned int i=0; i<this->gases.size(); i++)
    {
        if (this->gases[i]->name_material() == gas_type->name_material())
        {
            index = i;
            break;
        }
    }
}


//---------------------------------------------------------------------------
template <int dim>
void
NAME::PorousLayer<dim>::gas_diffusion_coefficients(Table< 2, double > &D_eff ) const
{
    D_eff = this->D_ECtheory;
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::PorousLayer<dim>::derivative_gas_diffusion_coefficients(std::vector<Table< 2, double > > &dD_dx) const
{
    dD_dx = this->dD_ECtheory_dx;
}


//---------------------------------------------------------------------------
template <int dim>
void
NAME::PorousLayer<dim>::Knudsen_diffusion(const FuelCellShop::Material::PureGas* solute_gas,
                                          const SolutionVariable& T_in,
                                          std::vector<double>& D_k) const
{
    
    D_k.resize(T_in.size());
    
    for (unsigned int q=0; q<T_in.size(); ++q)
        D_k[q] = (2.0*Knudsen_radius/3.0)*sqrt(8.0*Constants::R()*T_in[q]/(Constants::Pi()*solute_gas->get_molar_mass()));
    
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::PorousLayer<dim>::Knudsen_diffusion(const FuelCellShop::Material::PureGas* solute_gas,
                                          const SolutionVariable& T_in,
                                          std::vector<double>& D_k,
                                          std::vector<double>& dD_k_dT) const
{
    this->Knudsen_diffusion(solute_gas, T_in, D_k);
    
    dD_k_dT.resize(T_in.size());
    
    for (unsigned int q=0; q<T_in.size(); ++q)
    {
        dD_k_dT[q] = (Knudsen_radius/3)*(1/sqrt(8*Constants::R()*T_in[q]/(Constants::Pi()*solute_gas->get_molar_mass())))
        *(8*Constants::R()/(Constants::Pi()*solute_gas->get_molar_mass()));
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::PorousLayer<dim>::print_layer_properties() const
{
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "------------------------------";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    FcstUtilities::log << "Parameters for " << this->name << ":";
    FcstUtilities::log << std::endl;
    FcstUtilities::log << std::endl;
    for (unsigned int i = 0; i < this->material_ids.size(); i++)
        FcstUtilities::log << "Material ids: "   << this->material_ids.at(i)<< std::endl;
    FcstUtilities::log << std::endl;
    
    if( porosity_is_constant )
    {
        FcstUtilities::log << "Constant porosity: " << porosity;
        FcstUtilities::log << std::endl;
    }
    else
    {
        FcstUtilities::log << "Porosity is a function of space";
        FcstUtilities::log << std::endl;
    }
    
    if( permeability_is_constant )
    {
        FcstUtilities::log << "Constant permeability [m^2]:";
        FcstUtilities::log << std::endl;
        for(unsigned int i = 0; i < dim; ++i)
        {
            for(unsigned int j = 0; j < dim; ++j)
                FcstUtilities::log << permeability[i][j] << "   ";
            FcstUtilities::log << std::endl;
        }
        
        FcstUtilities::log << "Square root of constant permeability [m]:";
        FcstUtilities::log << std::endl;
        for(unsigned int i = 0; i < dim; ++i)
        {
            for(unsigned int j = 0; j < dim; ++j)
                FcstUtilities::log << SQRT_permeability[i][j] << "   ";
            FcstUtilities::log << std::endl;
        }
        
        FcstUtilities::log << "Inverse of constant permeability [1/m^2]:";
        FcstUtilities::log << std::endl;
        for(unsigned int i = 0; i < dim; ++i)
        {
            for(unsigned int j = 0; j < dim; ++j)
                FcstUtilities::log << permeability_INV[i][j] << "   ";
            FcstUtilities::log << std::endl;
        }
        
        FcstUtilities::log << "Inverse of square root of constant permeability [1/m]:";
        FcstUtilities::log << std::endl;
        for(unsigned int i = 0; i < dim; ++i)
        {
            for(unsigned int j = 0; j < dim; ++j)
                FcstUtilities::log << SQRT_permeability_INV[i][j] << "   ";
            FcstUtilities::log << std::endl;
        }
        
        FcstUtilities::log << "Constant Forchheimer permeability [1/m]:";
        FcstUtilities::log << std::endl;
        for(unsigned int i = 0; i < dim; ++i)
        {
            for(unsigned int j = 0; j < dim; ++j)
                FcstUtilities::log << Forchheimer_permeability[i][j] << "   ";
            FcstUtilities::log << std::endl;
        }
    }
    else
    {
        FcstUtilities::log << "Permeability is a function of space";
        FcstUtilities::log << std::endl;
    }
    
    if( tortuosity_is_constant )
    {
        FcstUtilities::log << "Constant tortuosity:";
        FcstUtilities::log << std::endl;
        for(unsigned int i = 0; i < dim; ++i)
        {
            for(unsigned int j = 0; j < dim; ++j)
                FcstUtilities::log << tortuosity[i][j] << "   ";
            FcstUtilities::log << std::endl;
        }
    }
    else
    {
        FcstUtilities::log << "Tortuosity is a function of space";
        FcstUtilities::log << std::endl;
    }
    
    if( gas_mixture )
        gas_mixture->print_material_properties();
}

// ---                   ---
// --- print_caller_name ---
// ---                   ---

template<int dim>
void
NAME::PorousLayer<dim>::print_caller_name(const std::string& caller_name) const
{
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << caller_name << " called in Class " << info.name() << std::endl;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::PorousLayer<deal_II_dimension>;