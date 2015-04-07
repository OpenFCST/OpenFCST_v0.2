//---------------------------------------------------------------------------
//    $Id: micro_porous_layer.cc 2605 2014-08-15 03:36:44Z secanell $
//
//    Copyright (C) 2009 by Marc Secanell, University of Alberta
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <micro_porous_layer.h>

namespace NAME = FuelCellShop::Layer;

template <int dim>
const std::string NAME::MicroPorousLayer<dim>::concrete_name ("MicroPorousLayer");

//---------------------------------------------------------------------------
template <int dim>
NAME::MicroPorousLayer<dim>::MicroPorousLayer()
  : NAME::PorousLayer<dim>()
{}

//---------------------------------------------------------------------------
template <int dim>
NAME::MicroPorousLayer<dim>::MicroPorousLayer(const std::string& name)
  : NAME::PorousLayer<dim>(name)
{}

//---------------------------------------------------------------------------
template <int dim>
NAME::MicroPorousLayer<dim>::~MicroPorousLayer()
{}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::MicroPorousLayer<dim>::declare_parameters (const std::string& name, 
                                                 ParameterHandler &param) const
{
    FuelCellShop::Layer::PorousLayer<dim>::declare_parameters(name,param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(name);
        {
            param.declare_entry("Micro porous layer type",
                                "SGL24BC",
                                Patterns::Selection("SGL24BC | DesignMPL"),
                                "The type of the Micro Porous Layer ");
            
            FuelCellShop::MicroScale::BasePSD<dim>::declare_PSD_parameters(param);
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}   

//---------------------------------------------------------------------------
template <int dim>
void
NAME::MicroPorousLayer<dim>::initialize (ParameterHandler &param)
{
  NAME::PorousLayer<dim>::initialize(param);
  
  param.enter_subsection("Fuel cell data"); 
  {
    param.enter_subsection(this->name); 
    { 
        param.enter_subsection("PSD parameters");
        {
            param.enter_subsection("BasePSD");
            {   
                PSD_type = param.get("psd type");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
        
        PSD = FuelCellShop::MicroScale::BasePSD<dim>::create_PSD(PSD_type,param);

    }
    param.leave_subsection();
  }
  param.leave_subsection();
}

//---------------------------------------------------------------------------
//-------------------- OPTIMIZATION ROUTINES --------------------------------
//---------------------------------------------------------------------------
/*
template <int dim>
Tensor<2,dim>
NAME::MicroPorousLayer<dim>::Deffective_transport_property_pores_Dporosity(const double prop) const
{

}
*/

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::MicroPorousLayer<deal_II_dimension>;
