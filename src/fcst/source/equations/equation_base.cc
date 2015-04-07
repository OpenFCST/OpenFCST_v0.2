// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: equation_base.cc
// - Description: This is a base class for all available openFCST equations
// - Developers: Valentin N. Zingan,    University of Alberta
//               Marc Secanell Gallart, University of Alberta
// - Id: $Id$
//
// ----------------------------------------------------------------------------

#include "equation_base.h"

namespace NAME = FuelCellShop::Equation;

       //////////////////////////////////////////////////
       //////////////////////////////////////////////////
       // CONSTRUCTORS, DESTRUCTOR, AND INITIALIZATION //
       //////////////////////////////////////////////////
       //////////////////////////////////////////////////

// ---             ---
// --- Constructor ---
// ---             ---

template<int dim>
NAME::EquationBase<dim>::EquationBase(FuelCell::SystemManagement& sys_management)
:
system_management(&sys_management)
{ }

// ---            ---
// --- Destructor ---
// ---            ---

template<int dim>
NAME::EquationBase<dim>::~EquationBase()
{ }

       ////////////////////
       ////////////////////
       // GENERIC INITIALIZATION //
       ////////////////////


//====================================
template<int dim>
void
NAME::EquationBase<dim>::declare_parameters(ParameterHandler& param) const
{
param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {
            //--
            param.enter_subsection("Initial data");
            {
                   param.declare_entry("Variable initial data",
                                       "false",
                                        Patterns::Bool(),
                                       "True, if initial data is NOT constant at least in one sub-domain. False, if initial data is constant or piece-wise constant.");

                   param.declare_entry( name_base_variable,
                                       """",
                                        Patterns::Map(   Patterns::Integer(0) , Patterns::Double()   ),
                                       "This information is used to setup an initial solution. Note this can be overwritten by Operating Conditions. \n"
                                       "Enter the mole fraction that you would like to use as initial value for each material in your mesh. \n"
                                       "The format should be as follows: material_id1:value1, material_id2:value2. For example, using "
                                       "2:0.1, 4:0.5 will set the initial mole fraction in material 2 to 0.1 and in material 4 to 0.5");
            }
            param.leave_subsection();
            //--
            param.enter_subsection("Boundary data");
            {
                   param.declare_entry("Variable boundary data",
                                       "false",
                                        Patterns::Bool(),
                                       "True, if Dirichlet boundary data is NOT constant at least at one sub-boundary. False, if Dirichlet boundary data is constant or piece-wise constant.");

                   param.declare_entry( name_base_variable,
                                       """",
                                        Patterns::Map(   Patterns::Integer(0) , Patterns::Double()   ),
                                       "This information is used to setup an initial solution. Note this can be overwritten by Operating Conditions.\n"
                                       "Enter the value that you would like to use as the Dirichlet boundary condition for each"
                                       "boundary_id. The format should be as follows: boundary_id1:value1, boundary_id2:value2. \n"
                                       "For example, using 3:0.1, 41:0.25 will set the value in boundary 3 to 0.1 and in boundary 41 to 0.25");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//======================================
template<int dim>
void
NAME::EquationBase<dim>::initialize(ParameterHandler& param)
{
    param.enter_subsection("Equations");
    {
        param.enter_subsection(this->equation_name);
        {
            param.enter_subsection("Initial data");
            {
                   variable_initial_data = param.get_bool("Variable initial data");

                   if( !param.get(name_base_variable).empty() )
                   {
                          const std::map<types::material_id, double> tmp = FcstUtilities::string_to_map<types::material_id, double>( param.get(name_base_variable) );
                          this->component_materialID_value[name_base_variable] = tmp;
                   }
            }
            param.leave_subsection();

            param.enter_subsection("Boundary data");
            {
                   variable_boundary_data = param.get_bool("Variable boundary data");

                   if( !param.get(name_base_variable).empty() )
                   {
                          const std::map<types::boundary_id, double> tmp = FcstUtilities::string_to_map<types::boundary_id, double>( param.get(name_base_variable) );
                          this->component_boundaryID_value[name_base_variable] = tmp;
                   }
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}


       ////////////////
       ////////////////
       // CONVERTERS //
       ////////////////
       ////////////////

// ---                        ---
// --- standard_to_block_wise ---
// ---                        ---

template<int dim>
void
NAME::EquationBase<dim>::standard_to_block_wise(FullMatrix<double>& target) const
{
  const std::vector<unsigned int> local_renumbering( this->system_management->block_info->local_renumbering );

  FullMatrix<double> tmp(target.m(), target.n());

  for(unsigned int i = 0; i < target.m(); ++i)
    for(unsigned int j = 0; j < target.n(); ++j)
      tmp(local_renumbering[i], local_renumbering[j]) = target(i,j);

  target = tmp;
}

// ---                        ---
// --- standard_to_block_wise ---
// ---                        ---

template<int dim>
void
NAME::EquationBase<dim>::standard_to_block_wise(Vector<double>& target) const
{
  const std::vector<unsigned int> local_renumbering( this->system_management->block_info->local_renumbering );

  Vector<double> tmp(target.size());

  for(unsigned int i = 0; i < target.size(); ++i)
    tmp(local_renumbering[i]) = target(i);

  target = tmp;
}

// ---                    ---
// --- dealII_to_appframe ---
// ---                    ---

template<int dim>
void
NAME::EquationBase<dim>::dealII_to_appframe(FuelCell::ApplicationCore::MatrixVector&          dst,
                                            const FullMatrix<double>&        src,
                                            const std::vector<unsigned int>& matrix_block_indices) const
{
  FullMatrix<double> block_wise(src);
  this->standard_to_block_wise(block_wise);

  BlockIndices local = this->system_management->block_info->local;
  std::vector<unsigned int> sizes(local.size());
  for(unsigned int i = 0; i < sizes.size(); ++i)
    sizes[i] = local.block_size(i);

  unsigned int number = 0;

  unsigned int i_from = 0;
  unsigned int i_to   = sizes[0] - 1;

  for(unsigned int line = 0; line < sizes.size(); ++line)
  {
    unsigned int j_from = 0;
    unsigned int j_to   = sizes[0] - 1;

    for(unsigned int column = 0; column < sizes.size(); ++column)
    {
           // ---

           const std::vector<unsigned int>::const_iterator
             iter = std::find(matrix_block_indices.begin(), matrix_block_indices.end(), number);

           if( iter == matrix_block_indices.end() )
           {
                  if( (*(this->system_management->cell_couplings))(line,column) != 0 )
                  {
                         if( number != dst.size() - 1 )
                         {
                                number++;
                         }
                  }
           }
           else
           {
                  if( (*(this->system_management->cell_couplings))(line,column) != 0 )
                  {
                         const unsigned int n_lines   = i_to - i_from + 1;
                         const unsigned int n_columns = j_to - j_from + 1;

                         FullMatrix<double> tmp(n_lines, n_columns);

                         unsigned int line_index = i_from;

                         for(unsigned int i = 0; i < n_lines; ++i)
                         {
                                unsigned int column_index = j_from;

                                for(unsigned int j = 0; j < n_columns; ++j)
                                {
                                       tmp(i,j) = block_wise(line_index, column_index);
                                       column_index++;
                                }

                                line_index++;
                         }

                         // new procedure which takes into account the accumulation

                         for(unsigned int i = 0; i < tmp.m(); ++i)
                                for(unsigned int j = 0; j < tmp.n(); ++j)
                                       dst[number].matrix(i,j) += tmp(i,j);

                         // dst[number].matrix = tmp; // it does not work if accumulation is going to take place
                                                      // typically it happens if the equations are assembled separately from sources

                         if( number != dst.size() - 1 )
                         {
                                number++;
                         }
                  }
           }

           if( column != sizes.size() - 1 )
           {
                  j_from = j_to + 1;
                  j_to   = j_to + sizes[column+1];
           }

           // ---
    }

    if( line != sizes.size() - 1 )
    {
           i_from = i_to + 1;
           i_to   = i_to + sizes[line+1];
    }
  }
}

// ---                    ---
// --- dealII_to_appframe ---
// ---                    ---

template<int dim>
void
NAME::EquationBase<dim>::dealII_to_appframe(FuelCell::ApplicationCore::FEVector&              dst,
                                            const Vector<double>&            src,
                                            const std::vector<unsigned int>& residual_indices) const
{
  AssertThrow( dst.size() == src.size() , ExcDimensionMismatch( dst.size() , src.size() ) );

  Vector<double> block_wise(src);
  this->standard_to_block_wise(block_wise);

  BlockIndices local = this->system_management->block_info->local;
  std::vector<unsigned int> sizes(local.size());
  for(unsigned int i = 0; i < sizes.size(); ++i)
    sizes[i] = local.block_size(i);

  unsigned int number = 0;

  unsigned int i_from = 0;
  unsigned int i_to   = sizes[0] - 1;

  for(unsigned int line = 0; line < sizes.size(); ++line)
  {
         const std::vector<unsigned int>::const_iterator
           iter = std::find(residual_indices.begin(), residual_indices.end(), number);

         if( iter != residual_indices.end() )
         {
                const unsigned int n_lines = i_to - i_from + 1;

                Vector<double> tmp(n_lines);

                unsigned int line_index = i_from;

                for(unsigned int i = 0; i < n_lines; ++i)
                {
                       tmp(i) = block_wise(line_index);
                       line_index++;
                }

                // new procedure which takes into account the accumulation

                for(unsigned int i = 0; i < tmp.size(); ++i)
                       dst.block(line)(i) += tmp(i);

                // dst.block(line) = tmp; // it does not work if accumulation is going to take place
                                          // typically it happens if the equations are assembled separately from sources
         }

         if( number != sizes.size() - 1 )
         {
                number++;
         }

         if( line != sizes.size() - 1 )
         {
                i_from = i_to + 1;
                i_to   = i_to + sizes[line+1];
         }
  }
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
NAME::EquationBase<dim>::print_caller_name(const std::string& caller_name) const
{
  const std::type_info& info = typeid(*this);
  FcstUtilities::log << "Pure function " << caller_name << " called in Class " << info.name() << std::endl;
}

       /////////////////////////////
       /////////////////////////////
       // EXPLICIT INSTANTIATIONS //
       /////////////////////////////
       /////////////////////////////

// ---              ---
// --- EquationBase ---
// ---              ---

template class NAME::EquationBase<deal_II_dimension>;