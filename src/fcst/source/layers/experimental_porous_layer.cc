// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2013 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT License
// For more information, see the README file in /doc/LICENSE
//
// - Class: experimental_porous_layer.cc
// - Description: This class describes a porous layer
// - Developers: Valentin N. Zingan, University of Alberta
// - Id: $Id: experimental_porous_layer.cc 2605 2014-08-15 03:36:44Z secanell $
//
// ----------------------------------------------------------------------------

#include "experimental_porous_layer.h"

namespace NAME = FuelCellShop::Layer;

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::ExperimentalPorousLayer<dim>::ExperimentalPorousLayer(const std::string& name)
:
NAME::BaseLayer<dim>(name),
fluid(nullptr),
multi_fluid(nullptr),
gas_mixture(nullptr),
solid(nullptr),
porosity_is_constant(true),
permeability_is_constant(true),
tortuosity_is_constant(true)
{ }

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::ExperimentalPorousLayer<dim>::ExperimentalPorousLayer(const std::string&                         name,
                                                            FuelCellShop::Material::ExperimentalFluid& fluid,
                                                            FuelCellShop::Material::ExperimentalSolid& solid)
:
NAME::BaseLayer<dim>(name),
fluid(&fluid),
multi_fluid(nullptr),
gas_mixture(nullptr),
solid(&solid),
porosity_is_constant(true),
permeability_is_constant(true),
tortuosity_is_constant(true)
{ }

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::ExperimentalPorousLayer<dim>::ExperimentalPorousLayer(const std::string&                              name,
                                                            FuelCellShop::Material::ExperimentalMultiFluid& multi_fluid,
                                                            FuelCellShop::Material::ExperimentalSolid&      solid)
:
NAME::BaseLayer<dim>(name),
fluid(nullptr),
multi_fluid(&multi_fluid),
gas_mixture(nullptr),
solid(&solid),
porosity_is_constant(true),
permeability_is_constant(true),
tortuosity_is_constant(true)
{ }

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::ExperimentalPorousLayer<dim>::ExperimentalPorousLayer(const std::string&                         name,
                                                            FuelCellShop::Material::GasMixture&        gas_mixture,
                                                            FuelCellShop::Material::ExperimentalSolid& solid)
:
NAME::BaseLayer<dim>(name),
fluid(nullptr),
multi_fluid(nullptr),
gas_mixture(&gas_mixture),
solid(&solid),
porosity_is_constant(true),
permeability_is_constant(true),
tortuosity_is_constant(true)
{ }

// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::ExperimentalPorousLayer<dim>::~ExperimentalPorousLayer()
{ }

// ---                    ---
// --- declare_parameters ---
// ---                    ---

template<int dim>
void
NAME::ExperimentalPorousLayer<dim>::declare_parameters(ParameterHandler& param) const
{
  NAME::BaseLayer<dim>::declare_parameters(this->name,
                                           param);

  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection(this->name);
    {
           param.declare_entry("Porosity",
                               "1.0",
                                Patterns::Double(),
                               " ");

           param.declare_entry("Permeability_XX",
                               "1.0",
                                Patterns::Double(),
                               " ");
           param.declare_entry("Permeability_XY",
                               "0.0",
                                Patterns::Double(),
                               " ");
           param.declare_entry("Permeability_XZ",
                               "0.0",
                                Patterns::Double(),
                               " ");

           param.declare_entry("Permeability_YX",
                               "0.0",
                                Patterns::Double(),
                               " ");
           param.declare_entry("Permeability_YY",
                               "1.0",
                                Patterns::Double(),
                               " ");
           param.declare_entry("Permeability_YZ",
                               "0.0",
                                Patterns::Double(),
                               " ");

           param.declare_entry("Permeability_ZX",
                               "0.0",
                                Patterns::Double(),
                               " ");
           param.declare_entry("Permeability_ZY",
                               "0.0",
                                Patterns::Double(),
                               " ");
           param.declare_entry("Permeability_ZZ",
                               "1.0",
                                Patterns::Double(),
                               " ");

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

// ---            ---
// --- initialize ---
// ---            ---

template<int dim>
void
NAME::ExperimentalPorousLayer<dim>::initialize(ParameterHandler& param)
{
  NAME::BaseLayer<dim>::initialize(param);

  param.enter_subsection("Fuel cell data");
  {
    param.enter_subsection(this->name);
    {
           porosity = param.get_double("Porosity");

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

       ////////////////////////
       ////////////////////////
       // ACCESSORS AND INFO //
       ////////////////////////
       ////////////////////////

// ---              ---
// --- get_porosity ---
// ---              ---

template<int dim>
void
NAME::ExperimentalPorousLayer<dim>::get_porosity(std::vector<double>& dst) const
{
  AssertThrow( porosity_is_constant , ExcInternalError() );

  for(unsigned int q = 0; q < dst.size(); ++q)
    dst[q] = porosity;
}

// ---              ---
// --- get_porosity ---
// ---              ---

template<int dim>
void
NAME::ExperimentalPorousLayer<dim>::get_porosity(std::vector<double>&             dst,
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
NAME::ExperimentalPorousLayer<dim>::get_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const
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
NAME::ExperimentalPorousLayer<dim>::get_permeability(std::vector< SymmetricTensor<2,dim> >& dst,
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
NAME::ExperimentalPorousLayer<dim>::get_SQRT_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const
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
NAME::ExperimentalPorousLayer<dim>::get_SQRT_permeability(std::vector< SymmetricTensor<2,dim> >& dst,
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
NAME::ExperimentalPorousLayer<dim>::get_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst) const
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
NAME::ExperimentalPorousLayer<dim>::get_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst,
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
NAME::ExperimentalPorousLayer<dim>::get_SQRT_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst) const
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
NAME::ExperimentalPorousLayer<dim>::get_SQRT_permeability_INV(std::vector< SymmetricTensor<2,dim> >& dst,
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
NAME::ExperimentalPorousLayer<dim>::get_Forchheimer_permeability(std::vector< SymmetricTensor<2,dim> >& dst) const
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
NAME::ExperimentalPorousLayer<dim>::get_Forchheimer_permeability(std::vector< SymmetricTensor<2,dim> >& dst,
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
NAME::ExperimentalPorousLayer<dim>::get_tortuosity(std::vector< SymmetricTensor<2,dim> >& dst) const
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
NAME::ExperimentalPorousLayer<dim>::get_tortuosity(std::vector< SymmetricTensor<2,dim> >& dst,
                                                   const std::vector< Point<dim> >&       points) const
{
  AssertThrow( !tortuosity_is_constant , ExcInternalError() );
  AssertThrow( dst.size() == points.size() , ExcDimensionMismatch(dst.size(), points.size()) );

  print_caller_name(__FUNCTION__);
}

// ---                        ---
// --- print_layer_properties ---
// ---                        ---

template<int dim>
void
NAME::ExperimentalPorousLayer<dim>::print_layer_properties() const
{
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "------------------------------";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;
  FcstUtilities::log << "Parameters for "  << this->name << ":";
  FcstUtilities::log << std::endl;
  FcstUtilities::log << std::endl;
  for(unsigned int i = 0; i < this->material_ids.size(); ++i)
    FcstUtilities::log << "Material ids: " << this->material_ids.at(i) << std::endl;
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

  if( fluid )
    fluid->print_material_properties();

  if( multi_fluid )
    multi_fluid->print_material_properties();

  if( gas_mixture )
    gas_mixture->print_material_properties();

  solid->print_material_properties();
}

       /////////////////////
       /////////////////////
       // MINOR FUNCTIONS //
       /////////////////////
       /////////////////////

// ---                   ---
// --- print_caller_name ---
// ---                   ---

template<int dim>
void
NAME::ExperimentalPorousLayer<dim>::print_caller_name(const std::string& caller_name) const
{
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << caller_name << " called in Class " << info.name() << std::endl;
}

       /////////////////////////////
       /////////////////////////////
       // EXPLICIT INSTANTIATIONS //
       /////////////////////////////
       /////////////////////////////

// ---                         ---
// --- ExperimentalPorousLayer ---
// ---                         ---

template class NAME::ExperimentalPorousLayer<deal_II_dimension>;