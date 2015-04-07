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

#include "PSD_HI.h"

namespace NAME = FuelCellShop::MicroScale;

template<int dim>
const std::string NAME::HIPSD<dim>::concrete_name("HIPSD");

template<int dim>
NAME::HIPSD<dim> const* NAME::HIPSD<dim>::PROTOTYPE = new NAME::HIPSD<dim>();

//---------------------------------------------------------------------------
template<int dim>
NAME::HIPSD<dim>::HIPSD() 
:
FuelCellShop::MicroScale::BasePSD<dim>() 
{
    this->get_mapFactory()->insert(std::pair<std::string, FuelCellShop::MicroScale::BasePSD<dim>* >(concrete_name, this));
} 

//---------------------------------------------------------------------------
template<int dim>
NAME::HIPSD<dim>::HIPSD(std::string name) 
:
NAME::BasePSD<dim>(name) 
{}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::HIPSD<dim>::declare_parameters (ParameterHandler &param) const
{
  FuelCellShop::MicroScale::BasePSD<dim>::declare_parameters(param);
  
  param.enter_subsection("PSD parameters");
  {
      param.enter_subsection("BasePSD");
      {
          param.enter_subsection(concrete_name); 
          {
            
            param.declare_entry("Hydrophilic Mode probability global",
                                "1.0",
                                Patterns::List(Patterns::Double(0.0)),
                                "Contribution of the distribution mode into the PSD ");
            
            param.declare_entry("Hydrophilic Mode characteristic radius global",
                                "1.0",
                                Patterns::List(Patterns::Double(0.0)),
                                "Characteristic pore size of the distribution mode into the PSD ");
            
            param.declare_entry("Hydrophilic Mode width global",
                                "1.0",
                                Patterns::List(Patterns::Double(0.0)),
                                "Characteristic pore size of the distribution mode into the PSD ");
            
              param.declare_entry("capillay pressure",
                                  "1.0",
                                  Patterns::Double(0.0),
                                  "capillay pressure in Pascal");

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
NAME::HIPSD<dim>::initialize (ParameterHandler &param)
{
    
    FuelCellShop::MicroScale::BasePSD<dim>::initialize(param);
    
    param.enter_subsection("PSD parameters");
    {
        param.enter_subsection("BasePSD");
        {
            param.enter_subsection(concrete_name); 
            {
                pressure_c = param.get_double("capillay pressure");
                
                fHI_k = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Hydrophilic Mode probability global") ) );
                
                rHI_k = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Hydrophilic Mode characteristic radius global") ) );
                
                sHI_k = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Hydrophilic Mode width global") ) );
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
//---------------------------------------------------------------------------

template<int dim>
void
NAME::HIPSD<dim>::get_critical_radius(std::vector<double>& critical_radius) const
{
    std::vector<double> p_c;
    p_c.clear();
    critical_radius.clear();
    
    if ( critical_radius_is_initialized )
    {
        critical_radius = critical_radius_computed;
    }
    else 
    {
        if (Capillary_pressure_vector.is_initialized ())
        {
            p_c.resize(this->Capillary_pressure_vector.size());
            
            for (unsigned int q = 0; q < p_c.size(); ++q)
                p_c[q] = this->Capillary_pressure_vector[q];
        }
        else
        {
            p_c.push_back (pressure_c);
        }
        
        critical_radius.resize(p_c.size());
        
        for(unsigned int q = 0; q < p_c.size() ; ++q)
        {
            
            critical_radius[q] =  1e6 * 2.0 * BasePSD<dim>::gamma * (std::cos(BasePSD<dim>::contact_angle)/p_c[q]); 
        }
    }
}

//---------------------------------------------------------------------------
template<int dim>
const double
NAME::HIPSD<dim>::get_maximum_cross_sectional_areas() const
{
    double a_max (0.0);
    
    for(unsigned int i = 0; i < this->r_k.size(); ++i)
        
        a_max += BasePSD<dim>::f_k[i] * std::exp (BasePSD<dim>::s_k[i] * BasePSD<dim>::s_k[i] /2.0 )/(4*BasePSD<dim>::r_k[i]);
    
    return a_max; 
}




// ---              ---
// --- get_saturation ---
// ---              ---

template<int dim>
void
NAME::HIPSD<dim>::get_saturation(std::vector<double>& S) const
{
    std::vector<double> critical_radius;
    S.clear();
    critical_radius.clear();
    
    if (saturation_is_initialized)
    {
        S = saturation_computed;
    }
    else
    {
        get_critical_radius(critical_radius);
        S.resize(critical_radius.size());
        
        for(unsigned int q = 0; q < critical_radius.size(); ++q)
            
            for(unsigned int i = 0; i < rHI_k.size(); ++i)
            {
                S [q] += BasePSD<dim>::F_HI 
                * 0.5 
                * fHI_k[i] * ( 1 + std::erf (   (std::log(critical_radius[q]) - std::log(rHI_k[i])  )/(  sHI_k[i] * std::pow(2,0.5)  ))); 
            }   
    }
}

// ---              ---
// --- get_pore_saturated_permeability ---
// ---              ---

template<int dim>
void
NAME::HIPSD<dim>::get_global_saturated_permeability(double& saturated_permeability) const
{
    for(unsigned int i = 0; i < fHI_k.size(); ++i)

        saturated_permeability +=  std::pow(this->get_porosity() / BasePSD<dim>::lamda , 2.0) 
                                   *( std::exp( (-2) * sHI_k[i] * sHI_k[i])      
                                   *rHI_k[i] * rHI_k[i] * fHI_k[i])   
                                   *1/8;

    
}
// ---              ---
// --- get_pore_liquid_permeability ---
// ---              ---

template<int dim>
void
NAME::HIPSD<dim>::get_pore_HI_liquid_saturated_permeability(std::vector<double>& saturated_HI_permeability) const
{
    std::vector<double> S;
    S.clear();
    get_saturation(S);
    
    std::vector<double> critical_radius;
    critical_radius.clear();
    get_critical_radius(critical_radius);
    
    saturated_HI_permeability.clear();
    
    saturated_HI_permeability.resize(S.size());
    
    for(unsigned int q = 0; q < S.size(); ++q)
        
        for(unsigned int i = 0; i < sHI_k.size(); ++i)
            
            saturated_HI_permeability [q] +=  BasePSD<dim>::F_HI 
                                              *std::pow(this->get_porosity()* S[q] /BasePSD<dim>::lamda,2.0) 
                                              *std::exp ( (-2) * sHI_k[i] * sHI_k[i] ) 
                                              * rHI_k[i] * rHI_k[i] * fHI_k[i] 
                                              *( std::erf ( (std::log(critical_radius[q]) - std::log(rHI_k[i])) / (sHI_k[i] * std::pow(2,0.5)) - sHI_k[i] * std::pow(2,0.5) ) + 1)  
                                              *1/16 ;
        
}

template<int dim>
void
NAME::HIPSD<dim>::get_relative_liquid_permeability(std::vector<double>& liquid_permeability) const
{
    std::vector<double> liquid_HI_permeability;
    liquid_HI_permeability.clear();
    get_pore_HI_liquid_saturated_permeability(liquid_HI_permeability);
    
    double saturated_permeability;
    get_global_saturated_permeability(saturated_permeability);
    
    liquid_permeability.clear();
    
    liquid_permeability.resize(liquid_HI_permeability.size());
    
    for(unsigned int q = 0; q < liquid_HI_permeability.size(); ++q)
        
        liquid_permeability[q] = liquid_HI_permeability[q] / saturated_permeability;
    
}

// ---              ---
// --- get_pore_gas_permeability ---
// ---              ---

template<int dim>
void
NAME::HIPSD<dim>::get_pore_HI_gas_saturated_permeability(std::vector<double>& saturated_HI_permeability) const
{
    std::vector<double> S;
    S.clear();
    get_saturation(S);
    
    std::vector<double> critical_radius;
    critical_radius.clear();
    get_critical_radius(critical_radius);
    
    saturated_HI_permeability.clear();
    saturated_HI_permeability.resize(S.size());
    
    for(unsigned int q = 0; q < S.size(); ++q)
        
        for(unsigned int i = 0; i < sHI_k.size(); ++i)
            
            saturated_HI_permeability[q] += BasePSD<dim>::F_HI
                                            * std::pow(this->get_porosity() * (1 - S[q])/BasePSD<dim>::lamda,2.0) 
                                            * (std::exp ((-2) * sHI_k[i] * sHI_k[i]) * rHI_k[i] * rHI_k[i] * fHI_k[i]) 
                                            * 1/16 
                                            * ( -std::erf ( (std::log(critical_radius[q]) - std::log(rHI_k[i]))/(sHI_k[i] * std::pow(2,0.5)) - sHI_k[i] * std::pow(2,0.5) ) + 1 );
        
        
}

template<int dim>
void
NAME::HIPSD<dim>::get_relative_gas_permeability(std::vector<double>& gas_permeability) const
{
    std::vector<double> HI_gas_permeability;
    get_pore_HI_gas_saturated_permeability(HI_gas_permeability);
    
    double saturated_permeability;
    get_global_saturated_permeability(saturated_permeability);
    
    gas_permeability.clear();
    gas_permeability.resize(HI_gas_permeability.size());
    
    for(unsigned int q = 0; q < HI_gas_permeability.size(); ++q)
        
        gas_permeability[q] = HI_gas_permeability[q]  / saturated_permeability;
}

// ---              ---
// --- get_pore_liquid_gas_interfacial_surface ---
// ---              ---

template<int dim>
void
NAME::HIPSD<dim>::get_liquid_gas_interfacial_surface(std::vector<double>& HI_liquid_gas_interfacial_surface) const
{
    std::vector<double> critical_radius;
    get_critical_radius(critical_radius);
    
    std::vector<double> HI_liquid_gas_interfacial_surface_a;
    HI_liquid_gas_interfacial_surface_a.clear();
    HI_liquid_gas_interfacial_surface_a.resize(critical_radius.size());
    
    HI_liquid_gas_interfacial_surface.clear();
    HI_liquid_gas_interfacial_surface.resize(critical_radius.size());
    
    for(unsigned int q = 0; q <critical_radius.size(); ++q)
    {
        
        for(unsigned int i = 0; i < sHI_k.size(); ++i)
        {
            
            HI_liquid_gas_interfacial_surface_a[q] +=   BasePSD<dim>::P_b 
                                                        * BasePSD<dim>::F_HI 
                                                        * fHI_k[i] 
                                                        * std::exp(sHI_k[i] * sHI_k[i] /2 )  
                                                        /  rHI_k[i]  
                                                        * ( std::erf ( (std::log(critical_radius[q]) - std::log(rHI_k[i])) / (sHI_k[i] * std::pow(2,0.5)) + sHI_k[i] * std::pow(2,0.5) / 2) + 1 )
                                                        / 8;
        }
        
        HI_liquid_gas_interfacial_surface [q] = HI_liquid_gas_interfacial_surface_a[q]/get_maximum_cross_sectional_areas() 
                                                * (1 - HI_liquid_gas_interfacial_surface_a[q]/get_maximum_cross_sectional_areas() ) 
                                                * HI_liquid_gas_interfacial_surface_a[q];
    }
}

// ---              ---
// --- get_pore_wetted_wall ---
// ---              ---

template<int dim>
void
NAME::HIPSD<dim>::get_pore_HI_wetted_wall_surface_area(std::vector<double>& HI_wetted_wall_surface_area) const
{
    std::vector<double> critical_radius;
    get_critical_radius(critical_radius);
    
    HI_wetted_wall_surface_area.clear();
    HI_wetted_wall_surface_area.resize(critical_radius.size());
    
    for(unsigned int i = 0; i < critical_radius.size(); ++i)
        
        for(unsigned int q = 0; q < fHI_k.size(); ++q)
            
            HI_wetted_wall_surface_area[i] = BasePSD<dim>::F_HI 
                                             * (fHI_k[q] * std::exp (sHI_k[q]*sHI_k[q]/2) / rHI_k[q]) 
                                             * (1 + std::erf ( (std::log(critical_radius[i]) - std::log (rHI_k[q])) / (sHI_k[q] * std::sqrt(2.0)) + sHI_k[q] / std::sqrt(2.0)));
        
}

template<int dim>
void
NAME::HIPSD<dim>::get_wetted_wall_surface_area(std::vector<double>& wetted_wall_surface_area) const
{
    std::vector<double> HI_wetted_wall_surface_area;
    get_pore_HI_wetted_wall_surface_area(HI_wetted_wall_surface_area);
    
    wetted_wall_surface_area.clear();
    wetted_wall_surface_area.resize(HI_wetted_wall_surface_area.size());
    
    for(unsigned int i = 0; i < HI_wetted_wall_surface_area.size(); ++i)
        
        wetted_wall_surface_area[i] = HI_wetted_wall_surface_area[i];
    
}
// ---              ---
// --- get_pore_knudsen_radius ---
// ---              ---

template<int dim>
void
NAME::HIPSD<dim>::get_pore_knudsen_radius_C1(std::vector<double>& knudsen_radius_C1) const
{
    
    std::vector<double> critical_radius;
    get_critical_radius(critical_radius);
    
    knudsen_radius_C1.clear();
    knudsen_radius_C1.resize(critical_radius.size());
    
    for(unsigned int i = 0; i < critical_radius.size(); ++i)
        
        for(unsigned int q = 0; q < fHI_k.size(); ++q)
        {
            knudsen_radius_C1 [i] += BasePSD<dim>::F_HI 
                                    * 0.5 
                                    * fHI_k[q] 
                                    *( 1 - std::erf ((std::log(critical_radius[i]) - std::log(rHI_k[q]))/(sHI_k[q] * std::pow(2,0.5)))); 
        }  
        
}

template<int dim>
void
NAME::HIPSD<dim>::get_pore_knudsen_radius_C3(std::vector<double>& knudsen_radius_C3) const
{
    std::vector<double> critical_radius;
    get_critical_radius(critical_radius);
    
    knudsen_radius_C3.clear();
    knudsen_radius_C3.resize(critical_radius.size());
    
    for(unsigned int i = 0; i < critical_radius.size(); ++i)
        
        for(unsigned int q = 0; q < fHI_k.size(); ++q)
            
            knudsen_radius_C3[i] += BasePSD<dim>::F_HI 
                                   * (fHI_k[q] * std::exp (sHI_k[q]*sHI_k[q]/2) / (rHI_k[q])) 
                                   * ( 1 - std::erf ( (std::log(critical_radius[i]) - std::log (rHI_k[q]) ) / (sHI_k[q] * std::sqrt(2.0)) + sHI_k[q] / std::sqrt(2.0)));
        
        
}

template<int dim>
void
NAME::HIPSD<dim>::get_knudsen_radius(std::vector<double>& knudsen_radius) const
{
    std::vector<double> knudsen_radius_C1;
    std::vector<double> knudsen_radius_C3;
    
    get_pore_knudsen_radius_C1(knudsen_radius_C1);
    get_pore_knudsen_radius_C3(knudsen_radius_C3);
    
    knudsen_radius.clear();
    knudsen_radius.resize(knudsen_radius_C1.size());
    
    for(unsigned int i = 0; i < knudsen_radius_C1.size(); ++i)
        
        knudsen_radius [i] = 1e-9 * knudsen_radius_C1[i] / knudsen_radius_C3[i];
    
    
}
// ---              ---
// --- get_pore_diffusivity ---
// ---              ---

template<int dim>
void
NAME::HIPSD<dim>::get_diffusivity() const
{
    
    
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations.
template class NAME::HIPSD<deal_II_dimension>;