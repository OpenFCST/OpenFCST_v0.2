//---------------------------------------------------------------------------
//    $Id: operating_conditions.cc 2605 2014-08-15 03:36:44Z secanell $
//
//    Copyright (C) 2009 by Marc Secanell
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include "operating_conditions.h"

//---------------------------------------------------------------------------

FuelCell::OperatingConditions::OperatingConditions()
{}

//---------------------------------------------------------------------------
FuelCell::OperatingConditions::~OperatingConditions()
{}

//---------------------------------------------------------------------------
void
FuelCell::OperatingConditions::declare_parameters (ParameterHandler &param) const
{
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Operating conditions");
        {
            param.declare_entry("Adjust initial solution and boundary conditions",
                                "false",
                                Patterns::Bool(),
                                "Use the parameters in Operating conditions to create an initial solution"
                                "and overwrite the boundary conditions for the problem specified in Equations>>Initial Data"
                                "using the parameters in this section");
            param.declare_entry ("Temperature cell",
                                "353", // K, or 80 Celsius
                                Patterns::Double());
            param.declare_entry ("Cathode pressure",
                                 "101325", // Pa, or 1 atm
                                 Patterns::Double(),
                                 "Pressure at the cathode channel in Pa");
            param.declare_entry ("Cathode initial oxygen mole fraction (prior to humidification)",
                                 "0.21", // 
                                 Patterns::Double());
            param.declare_entry ("Cathode relative humidity",
                                 "0.7",
                                 Patterns::Double(),
                                 "Relative humidity (fraction) in the cathode channel");
            param.declare_entry ("Anode pressure",
                                 "101325", // Pa, or 1 atm
                                 Patterns::Double(),
                                 "Pressure at the anode channel in Pa");
            param.declare_entry ("Anode relative humidity",
                                 "0.7",
                                 Patterns::Double(),
                                 "Relative humidity (fraction) in the anode channel");
            param.declare_entry ("Voltage cell",
                                 "0.6", 
                                 Patterns::Double());
            param.declare_entry("Voltage drop in the anode",
                                "0.015",
                                Patterns::Double());
            param.declare_entry ("Open circuit voltage",
                                 "1.23",
                                 Patterns::Double());
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
void
FuelCell::OperatingConditions::initialize (ParameterHandler& param)
{
    R = Constants::R();
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Operating conditions");
        {
            adjust_BC = param.get_bool("Adjust initial solution and boundary conditions");
            T_cell = param.get_double("Temperature cell");
            // Cathode
            p_c = param.get_double("Cathode pressure");
            c_c = p_c/(R*T_cell)*1E-6; //mol/cm^3
            RH_c = param.get_double("Cathode relative humidity");
            channel_oxygen_mole_fraction = param.get_double("Cathode initial oxygen mole fraction (prior to humidification)");
            // Anode
            p_a = param.get_double("Anode pressure");
            c_a = p_a/(R*T_cell)*1E-6; //mol/cm^3
            RH_a = param.get_double("Anode relative humidity");
            // Cell voltage
            V_cell = std::fabs(param.get_double("Voltage cell"));
            dV_a = param.get_double("Voltage drop in the anode");
            OCV = param.get_double("Open circuit voltage");
            
            if (OCV > voltage_cell_th()) // NOTE: This also initializes E_th
                OCV = E_th;
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}
//---------------------------------------------------------------------------
void
FuelCell::OperatingConditions::adjust_initial_solution(std::vector< component_materialID_value_map >& maps,
                                                       const boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid) const
{
    for(unsigned int i = 0; i < maps.size(); ++i)
    {
        for(component_materialID_value_map::iterator iter  = maps[i].begin(); iter != maps[i].end(); ++iter)
        {                       
            if ((iter->first.compare("oxygen_molar_fraction") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Cathode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_o2();
                
                material_ID = grid->get_material_id("Cathode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_o2();
                
                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_o2();
            }
            //
            else if ((iter->first.compare("water_molar_fraction") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Cathode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                   iter->second[material_ID[ind]] = this->get_x_wv();

                material_ID = grid->get_material_id("Cathode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_x_wv();

                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)            
                    iter->second[material_ID[ind]] = this->get_x_wv();
                
                material_ID = grid->get_material_id("Anode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = 1 - this->get_x_h2();
                
                material_ID = grid->get_material_id("Anode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = 1 - this->get_x_h2();
                
                material_ID = grid->get_material_id("Anode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = 1 - this->get_x_h2();

            }
            //
            else if ((iter->first.compare("electronic_electrical_potential") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Cathode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                   iter->second[material_ID[ind]] = this->get_V();      
                
                material_ID = grid->get_material_id("Cathode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_V(); 
                
                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_V(); 

            }
            //
            else if ((iter->first.compare("protonic_electrical_potential") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Membrane");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = (this->get_V() - OCV)/4.0;
                
                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = (this->get_V() - OCV)/2.0; 

                material_ID = grid->get_material_id("Anode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = -0.0001;

            }
            //
            else if ((iter->first.compare("membrane_water_content") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Membrane");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = 5.0; 
                
                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = 4.0; 
                
               material_ID = grid->get_material_id("Anode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = 7.0;
            }
            //
            else if ((iter->first.compare("temperature_of_REV") == 0) && adjust_BC)
            {
                std::vector<unsigned int> material_ID = grid->get_material_id("Cathode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                   iter->second[material_ID[ind]] = this->get_T();

                material_ID = grid->get_material_id("Cathode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_T();

                material_ID = grid->get_material_id("Cathode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)            
                    iter->second[material_ID[ind]] = this->get_T();
                
                material_ID = grid->get_material_id("Membrane");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_T();
                
                material_ID = grid->get_material_id("Anode CL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_T();
                
                material_ID = grid->get_material_id("Anode MPL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_T();
                
                material_ID = grid->get_material_id("Anode GDL");
                for (unsigned int ind = 0; ind<material_ID.size(); ind++)
                    iter->second[material_ID[ind]] = this->get_T();

            }
        }
    }
}
//---------------------------------------------------------------------------
void
FuelCell::OperatingConditions::adjust_boundary_conditions(std::vector< component_boundaryID_value_map >& maps,
                                                          const boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid) const
{
    for(unsigned int i = 0; i < maps.size(); ++i)
        for(component_boundaryID_value_map::iterator iter  = maps[i].begin(); iter != maps[i].end(); ++iter)
        {                       
            if ((iter->first.compare("oxygen_molar_fraction") == 0) && adjust_BC)
            {
                unsigned int boundary_ID = grid->get_boundary_id("c_Ch/GDL");
                iter->second[boundary_ID] = this->get_x_o2();
                boundary_ID = grid->get_boundary_id("a_GDL/Ch");
                iter->second[boundary_ID] = 0.0;
            }
            else if ((iter->first.compare("water_molar_fraction") == 0) && adjust_BC)
            {
                unsigned int boundary_ID = grid->get_boundary_id("c_Ch/GDL");
                iter->second[boundary_ID] = this->get_x_wv();
                boundary_ID = grid->get_boundary_id("a_GDL/Ch");
                iter->second[boundary_ID] = 1 - this->get_x_h2();
            }
            else if ((iter->first.compare("electronic_electrical_potential") == 0) && adjust_BC)
            {
                unsigned int boundary_ID = grid->get_boundary_id("c_BPP/GDL");
                iter->second[boundary_ID] = this->get_V();  
                boundary_ID = grid->get_boundary_id("a_GDL/BPP");
                iter->second[boundary_ID] = 0.0;
            }
            else if ((iter->first.compare("temperature_of_REV") == 0) && adjust_BC)
            {
                unsigned int boundary_ID = grid->get_boundary_id("c_BPP/GDL");
                iter->second[boundary_ID] = this->get_T();  
                boundary_ID = grid->get_boundary_id("a_GDL/BPP");
                iter->second[boundary_ID] = this->get_T();
            }
        }
}

//---------------------------------------------------------------------------
double
FuelCell::OperatingConditions::saturation_pressure() const
{
  double T_celsius = T_cell - 273;
  return pow(10,(-2.1794+0.02953*T_celsius-0.000091837*pow(T_celsius,2)+0.00000014454*pow(T_celsius,3)));
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_x_wv() const 
{
  return RH_c*saturation_pressure()*101325/p_c; 
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_x_o2() const
{
  return channel_oxygen_mole_fraction*(p_c-saturation_pressure()*RH_c*101325)/p_c; 
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::get_x_h2() const
{
  return (p_a-saturation_pressure()*RH_a*101325)/p_a; 
}

//---------------------------------------------------------------------------
double 
FuelCell::OperatingConditions::voltage_cell_th()
{
  double p_h2 = p_a*get_x_h2()/101325;
  double p_o2 = p_c*get_x_o2()/101325;
  E_th  = 1.229 - (T_cell-298.15)*8.456E-4 + T_cell*4.31E-5*(log(p_h2)+0.5*log(p_o2));
  return E_th;
}


//---------------------------------------------------------------------------
void 
FuelCell::OperatingConditions::print_operating_conditions() const
{
	FcstUtilities::log<<"========= OPERATING CONDITIONS ========"<<std::endl;
	FcstUtilities::log<<"Temperature: "<< T_cell <<std::endl;
	FcstUtilities::log<<"Cathode Pressure: "<< p_c <<std::endl;
	FcstUtilities::log<<"Cathode RH: "<< RH_c <<std::endl;
	FcstUtilities::log<<"Anode Pressure: "<< p_a <<std::endl;
	FcstUtilities::log<<"Anode RH: "<< RH_a <<std::endl;
	FcstUtilities::log<<"Open Circuit Voltage: "<< OCV <<std::endl;
	FcstUtilities::log<<"Cell Voltage: "<< V_cell <<std::endl;
	FcstUtilities::log<<"======================================="<<std::endl;
}