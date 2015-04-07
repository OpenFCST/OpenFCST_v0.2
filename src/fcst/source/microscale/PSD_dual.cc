//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PSD_base.h
//    - Description: Base class for pore size distribution model.
//    - Developers: 2009-13 by Marc Secanell, University of Alberta
//                  2013-14 by Jie Zhou, University of Alberta
//    - $ $
//
//---------------------------------------------------------------------------

#include "PSD_dual.h"

namespace NAME = FuelCellShop::MicroScale;

template<int dim>
const std::string NAME::DualPSD<dim>::concrete_name("DualPSD");

template<int dim>
NAME::DualPSD<dim> const* NAME::DualPSD<dim>::PROTOTYPE = new NAME::DualPSD<dim>();

//---------------------------------------------------------------------------
template<int dim>
NAME::DualPSD<dim>::DualPSD() :
NAME::BasePSD<dim>() 
{
    this->get_mapFactory()->insert(
            std::pair<std::string, FuelCellShop::MicroScale::BasePSD<dim>*>(
                    concrete_name, this));
}

//---------------------------------------------------------------------------
template<int dim>
NAME::DualPSD<dim>::DualPSD(std::string name) :
NAME::BasePSD<dim>(name) 
{
}

//---------------------------------------------------------------------------
template<int dim>
NAME::DualPSD<dim>::~DualPSD() 
{
}

//---------------------------------------------------------------------------

template <int dim>
void
NAME::DualPSD<dim>::declare_parameters (ParameterHandler &param) const
{
  FuelCellShop::MicroScale::BasePSD<dim>::declare_parameters(param);
  
  psd_hi->declare_parameters(param);
  
  psd_ho->declare_parameters(param);
  
  param.enter_subsection("PSD parameters");
  {
      param.enter_subsection("BasePSD");
      {
          param.enter_subsection(concrete_name); 
          {
              
              
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
NAME::DualPSD<dim>::initialize (ParameterHandler &param)
{
    
    psd_hi->initialize(param);
    
    psd_ho->initialize(param);
    
    param.enter_subsection("PSD parameters");
    {
        param.enter_subsection("BasePSD");
        {
            
            param.enter_subsection(concrete_name); 
            {
                
                
            }
            param.leave_subsection(); 
        }
        param.leave_subsection();
    }
  param.leave_subsection();
  
}


// ---              ---
// --- get_critical_radius ---
// ---              ---


template<int dim>
const double
NAME::DualPSD<dim>::get_critical_radius(double& p_c) const
{
    return 1e6 * 2.0 * BasePSD<dim>::gamma * (std::cos( BasePSD<dim>::contact_angle )/p_c); 
}

template<int dim>
void
NAME::DualPSD<dim>::get_critical_radius(std::vector<double>& critical_radius) const
{
    std::vector<double> p_c;
    p_c.clear();
    critical_radius.clear();
    p_c.resize(this->Capillary_pressure_vector.size());
    
    for (unsigned int q = 0; q < p_c.size(); ++q)
        p_c[q] = this->Capillary_pressure_vector[q];

    
    for(unsigned int q = 0; q < p_c.size(); ++q)
    {
        p_c[q] = this->Capillary_pressure_vector[q];
        
        critical_radius[q] =  1e6 * 2.0 * BasePSD<dim>::gamma*(std::cos(BasePSD<dim>::contact_angle)/p_c[q]); 
    }
}

// ---              ---
// --- get_saturation ---
// ---              ---

template<int dim>
void
NAME::DualPSD<dim>::get_saturation(std::vector<double>& S) const
{
    std::vector<double> S_HO;
    psd_ho->get_saturation(S_HO);
    
    std::vector<double> S_HI;
    psd_hi->get_saturation(S_HI); 
    
    S.clear();
    S.resize(S_HI.size());
    
    for(unsigned int i = 0; i < S_HO.size(); ++i)
        
        S[i] = S_HI[i]+S_HO[i];
    
    
}

// ---              ---
// --- get_pore_saturated_permeability ---
// ---              ---

template<int dim>
void
NAME::DualPSD<dim>::get_global_saturated_permeability(double& saturated_permeability) const
{
    for(unsigned int i = 0; i < BasePSD<dim>::f_k.size(); ++i)
        
        saturated_permeability = 1/8 
                                 * std::pow(this->get_porosity()/BasePSD<dim>::lamda,2.0) 
                                 * (std::exp ((-2) * BasePSD<dim>::s_k[i] * BasePSD<dim>::s_k[i]) * BasePSD<dim>::r_k[i] * BasePSD<dim>::r_k[i] * BasePSD<dim>::f_k[i]);
    
}
// ---              ---
// --- get_pore_liquid_permeability ---
// ---              ---

template<int dim>
void
NAME::DualPSD<dim>::get_relative_liquid_permeability(std::vector<double>& liquid_permeability) const
{
    std::vector<double> liquid_HI_permeability;
    psd_hi->get_pore_HI_liquid_saturated_permeability(liquid_HI_permeability);
    
    std::vector<double> liquid_HO_permeability;
    psd_ho->get_pore_HO_liquid_saturated_permeability(liquid_HO_permeability);
    
    double saturated_permeability;
    get_global_saturated_permeability(saturated_permeability);
    
    liquid_permeability.clear();
    liquid_permeability.resize(liquid_HI_permeability.size());
    
    for(unsigned int q = 0; q < liquid_HI_permeability.size(); ++q)
        
        liquid_permeability[q] = (liquid_HI_permeability[q] + liquid_HO_permeability[q]) / saturated_permeability;
    
}
// ---              ---
// --- get_pore_gas_permeability ---
// ---              ---

template<int dim>
void
NAME::DualPSD<dim>::get_relative_gas_permeability(std::vector<double>& gas_permeability) const
{
    std::vector<double> HI_gas_permeability;
    psd_hi->get_pore_HI_liquid_saturated_permeability(HI_gas_permeability);
    
    std::vector<double> HO_gas_permeability;
    psd_ho->get_pore_HO_gas_saturated_permeability(HO_gas_permeability);
    
    double saturated_permeability;
    get_global_saturated_permeability(saturated_permeability);
    
    gas_permeability.clear();
    gas_permeability.resize(HI_gas_permeability.size());
    
    for(unsigned int q = 0; q < HI_gas_permeability.size(); ++q)
        
        gas_permeability[q] = (HI_gas_permeability[q] + HO_gas_permeability[q]) / saturated_permeability;
}
// ---              ---
// --- get_pore_liquid_gas_interfacial_surface ---
// ---              ---

template<int dim>
void
NAME::DualPSD<dim>::get_liquid_gas_interfacial_surface(std::vector<double>& liquid_gas_interfacial_surface) const
{
    std::vector<double> HI_liquid_gas_interfacial_surface;
    std::vector<double> HO_liquid_gas_interfacial_surface;
    
    psd_hi->get_liquid_gas_interfacial_surface (HI_liquid_gas_interfacial_surface);
    psd_ho->get_liquid_gas_interfacial_surface (HO_liquid_gas_interfacial_surface);
    
    liquid_gas_interfacial_surface.clear();
    liquid_gas_interfacial_surface.resize(HI_liquid_gas_interfacial_surface.size());
    
    for(unsigned int i = 0; i < HI_liquid_gas_interfacial_surface.size(); ++i)
        
        liquid_gas_interfacial_surface[i] = HI_liquid_gas_interfacial_surface[i] + HO_liquid_gas_interfacial_surface[i];
    
}
// ---              ---
// --- get_pore_wetted_wall ---
// ---              ---

template<int dim>
void
NAME::DualPSD<dim>::get_wetted_wall_surface_area(std::vector<double>& wetted_wall_surface_area) const
{
    std::vector<double> HI_wetted_wall_surface_area;
    std::vector<double> HO_wetted_wall_surface_area;
    
    psd_hi->get_pore_HI_wetted_wall_surface_area(HI_wetted_wall_surface_area);
    psd_ho->get_pore_HO_wetted_wall_surface_area(HO_wetted_wall_surface_area);
    
    wetted_wall_surface_area.clear();
    wetted_wall_surface_area.resize(HI_wetted_wall_surface_area.size());
    
    for(unsigned int i = 0; i < HI_wetted_wall_surface_area.size(); ++i)
        
        wetted_wall_surface_area[i] = HI_wetted_wall_surface_area[i] + HO_wetted_wall_surface_area[i];
    
}
// ---              ---
// --- get_pore_knudsen_radius ---
// ---              ---
template<int dim>
void
NAME::DualPSD<dim>::get_knudsen_radius(std::vector<double>& knudsen_radius) const
{
    std::vector<double> knudsen_radius_C1;
    std::vector<double> knudsen_radius_C2;
    std::vector<double> knudsen_radius_C3;
    std::vector<double> knudsen_radius_C4;
    
    psd_hi->get_pore_knudsen_radius_C1(knudsen_radius_C1);
    psd_ho->get_pore_knudsen_radius_C2(knudsen_radius_C2);
    psd_hi->get_pore_knudsen_radius_C3(knudsen_radius_C3);
    psd_ho->get_pore_knudsen_radius_C4(knudsen_radius_C4);
    
    knudsen_radius.clear();
    knudsen_radius.resize(knudsen_radius_C1.size());
    
    for(unsigned int i = 0; i < knudsen_radius_C1.size(); ++i)
        
        knudsen_radius [i] = (knudsen_radius_C1[i] + knudsen_radius_C2[i]) / (knudsen_radius_C3[i] + knudsen_radius_C4[i]);
    
    
}
// ---              ---
// --- get_pore_diffusivity ---
// ---              ---

template<int dim>
void
NAME::DualPSD<dim>::get_diffusivity() const
{
    
    
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::DualPSD<deal_II_dimension>;