//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2012, 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: dummy_GDL.cc
//    - Description: Implementation of a GDL class that setup us all properties from file
//    - Developers: M. Secanell
//    - Id: $Id: dummy_GDL.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <dummy_GDL.h>

namespace NAME = FuelCellShop::Layer; 

template <int dim>
const std::string NAME::DummyGDL<dim>::concrete_name ("DummyGDL");

template <int dim>
NAME::DummyGDL<dim> const* NAME::DummyGDL<dim>::PROTOTYPE = new NAME::DummyGDL<dim>();


//---------------------------------------------------------------------------
template <int dim>
NAME::DummyGDL<dim>::DummyGDL()
: FuelCellShop::Layer::GasDiffusionLayer<dim> ()
{
    //FcstUtilities::log<<" Register DummyGDL GDL to FactoryMap"<<std::endl;
    this->get_mapFactory()->insert(std::pair<std::string, FuelCellShop::Layer::GasDiffusionLayer<dim>* >(concrete_name, this));
}
            
            
//---------------------------------------------------------------------------
template <int dim>
NAME::DummyGDL<dim>::DummyGDL(std::string name)
  : NAME::GasDiffusionLayer<dim> (name)
{
  FcstUtilities::log<<" of type DummyGDL"<<std::endl;
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::declare_parameters (const std::string& name, 
                                         ParameterHandler &param) const
{
    
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_parameters(name, param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection(name);
        {
            param.enter_subsection(concrete_name); //-- Transport for the anisotropic case:
            {
                
                // ISOTROPIC PROPERTIES:
                param.declare_entry ("Oxygen diffusion coefficient, [cm^2/s]",
                                     "0.2741", //1atm, 353K
                                     Patterns::Double(),
                                     "Oxygen diffusion coefficient given by experiment");
                param.declare_entry ("Water vapour diffusion coefficient, [cm^2/s]",
                                     "0.29646", 
                                     Patterns::Double(),
                                     "Water vapour diffusion coefficient given by experiment");
                param.declare_entry ("Electrical conductivity, [S/cm]",
                                     "100", // [S/cm]
                                     Patterns::Double(),
                                     "Effective cond. if given is used, otherwise conducitivity of the raw material. Units [S/cm]");
                param.declare_entry ("Thermal conductivity, [W/(cm K)]",
                                     "16", // [S/cm]
                                     Patterns::Double(),
                                     "Effective thermal cond. if given is used, otherwise conducitivity of the raw material. Units [W/(cm K)]");
                // ANISOTROPIC PROPERTIES IF NEEDED:    
                param.declare_entry ("Anisotropic transport",
                                     "false", // [S/cm]
                                     Patterns::Bool(),
                                     "Boolean variable. Set to true if we want to account for anisotropy of the GDL");
                //--- XX
                param.declare_entry ("Oxygen diffusion coefficient X, [cm^2/s]",
                                     "0.2741", //1atm, 353K
                                     Patterns::Double(),
                                     "Oxygen diffusion coefficient given by experiment");
                param.declare_entry ("Water vapour diffusion coefficient X, [cm^2/s]",
                                     "0.29646", 
                                     Patterns::Double(),
                                     "Water vapour diffusion coefficient given by experiment");
                param.declare_entry ("Electrical conductivity X, [S/cm]",
                                     "100", // [S/cm]
                                     Patterns::Double(),
                                     "Component X of the electrical conductivity tensor. Units [S/cm]");
                param.declare_entry ("Thermal conductivity X, [W/(cm K)]",
                                     "16", // [S/cm]
                                     Patterns::Double(),
                                     "Effective thermal cond. if given is used, otherwise conducitivity of the raw material. Units [W/(cm K)]");
                // YY
                param.declare_entry ("Oxygen diffusion coefficient Y, [cm^2/s]",
                                     "0.2741", //1atm, 353K
                                     Patterns::Double(),
                                     "Oxygen diffusion coefficient given by experiment");
                param.declare_entry ("Water vapour diffusion coefficient Y, [cm^2/s]",
                                     "0.29646", 
                                     Patterns::Double(),
                                     "Water vapour diffusion coefficient given by experiment");
                param.declare_entry ("Electrical conductivity Y, [S/cm]",
                                     "100", // [S/cm]
                                     Patterns::Double(),
                                     "Component Y of the electrical conductivity tensor. Units [S/cm]");
                param.declare_entry ("Thermal conductivity Y, [W/(cm K)]",
                                     "16", // [S/cm]
                                     Patterns::Double(),
                                     "Effective thermal cond. if given is used, otherwise conducitivity of the raw material. Units [W/(cm K)]");
                // ZZ
                param.declare_entry ("Oxygen diffusion coefficient Z, [cm^2/s]",
                                     "0.2741", //1atm, 353K
                                     Patterns::Double(),
                                     "Oxygen diffusion coefficient given by experiment");
                param.declare_entry ("Water vapour diffusion coefficient Z, [cm^2/s]",
                                     "0.29646", 
                                     Patterns::Double(),
                                     "Water vapour diffusion coefficient given by experiment");
                param.declare_entry ("Electrical conductivity Z, [S/cm]",
                                     "100", 
                                     Patterns::Double(),
                                     "Component Z of the electrical conductivity tensor. Units [S/cm]");
                param.declare_entry ("Thermal conductivity Z, [W/(cm K)]",
                                     "16", // [S/cm]
                                     Patterns::Double(),
                                     "Effective thermal cond. if given is used, otherwise conducitivity of the raw material. Units [W/(cm K)]");
                
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
NAME::DummyGDL<dim>::initialize (ParameterHandler &param)
{
    
    NAME::GasDiffusionLayer<dim>::initialize(param);
    
    param.enter_subsection("Fuel cell data"); 
    {
        param.enter_subsection(this->name); 
        { 
            param.enter_subsection(concrete_name); 
            { 
                // Anisotropy
                anisotropy = param.get_bool ("Anisotropic transport");
                
                D_O2.resize(dim);
                D_wv.resize(dim);
                sigma_e.resize(dim);
                k_T.resize(dim);
                
                if (anisotropy == true)
                {
                    switch (dim)
                    {
                    // X
                    case 1:
                        D_O2[0] = param.get_double("Oxygen diffusion coefficient X, [cm^2/s]");
                        D_wv[0] = param.get_double("Water vapour diffusion coefficient X, [cm^2/s]");
                        sigma_e[0] = param.get_double("Electrical conductivity X, [S/cm]");	  
                        k_T[0] = param.get_double("Thermal conductivity X, [W/(cm K)]");                        
                        break;

                    case 2:
                        D_O2[0] = param.get_double("Oxygen diffusion coefficient X, [cm^2/s]");
                        D_wv[0] = param.get_double("Water vapour diffusion coefficient X, [cm^2/s]");
                        sigma_e[0] = param.get_double("Electrical conductivity X, [S/cm]");   
                        k_T[0] = param.get_double("Thermal conductivity X, [W/(cm K)]");
                        // Y
                        D_O2[1] = param.get_double("Oxygen diffusion coefficient Y, [cm^2/s]");
                        D_wv[1] = param.get_double("Water vapour diffusion coefficient Y, [cm^2/s]");
                        sigma_e[1] = param.get_double("Electrical conductivity Y, [S/cm]");
                        k_T[1] = param.get_double("Thermal conductivity Y, [W/(cm K)]");
                        break;
                        
                    case 3:
                        D_O2[0] = param.get_double("Oxygen diffusion coefficient X, [cm^2/s]");
                        D_wv[0] = param.get_double("Water vapour diffusion coefficient X, [cm^2/s]");
                        sigma_e[0] = param.get_double("Electrical conductivity X, [S/cm]");   
                        k_T[0] = param.get_double("Thermal conductivity X, [W/(cm K)]");
                        // Y
                        D_O2[1] = param.get_double("Oxygen diffusion coefficient Y, [cm^2/s]");
                        D_wv[1] = param.get_double("Water vapour diffusion coefficient Y, [cm^2/s]");
                        sigma_e[1] = param.get_double("Electrical conductivity Y, [S/cm]");
                        k_T[1] = param.get_double("Thermal conductivity Y, [W/(cm K)]");
                        // Z
                        D_O2[2] = param.get_double("Oxygen diffusion coefficient Z, [cm^2/s]");	  
                        D_wv[2] = param.get_double("Water vapour diffusion coefficient Z, [cm^2/s]");
                        sigma_e[2] = param.get_double("Electrical conductivity Z, [S/cm]");
                        k_T[2] = param.get_double("Thermal conductivity Z, [W/(cm K)]");
                        break;
                    }
                }
                else
                {
                    for (unsigned int i =0; i < D_O2.size(); i++)
                    {
                        D_O2[i] = param.get_double("Oxygen diffusion coefficient, [cm^2/s]");
                        D_wv[i] = param.get_double("Water vapour diffusion coefficient, [cm^2/s]");
                        sigma_e[i] = param.get_double("Electrical conductivity, [S/cm]");
                        k_T[i] = param.get_double("Thermal conductivity, [W/(cm K)]");
                    }
                }
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
NAME::DummyGDL<dim>::effective_gas_diffusivity(Table< 2, double> & prop_eff) const
{
    prop_eff.reinit(this->gases.size(),this->gases.size());
    
    unsigned int dimension = 0;

    for (unsigned int i = 0; i<this->gases.size(); i++)
    {
        std::string solute_name(this->gases[i]->name_material());

        for (unsigned int j = i+1; j<this->gases.size(); j++)
        {  
            std::string solvent_name(this->gases[j]->name_material());
            
            if (solute_name.compare("oxygen") == 0)
            {
                
                if (solvent_name.compare("nitrogen") == 0)
                {
                    prop_eff(i,j) = Units::convert(D_O2[dimension], Units::UNIT2, Units::C_UNIT2);
                    prop_eff(j,i) = prop_eff(i,j);
                    break;
                }
                
            }
            else if (solute_name.compare("water") == 0)
            {
                
                if (solvent_name.compare("nitrogen") == 0)
                {
                    prop_eff(i,j) = Units::convert(D_wv[dimension], Units::UNIT2, Units::C_UNIT2);
                    prop_eff(j,i) = prop_eff(i,j);
                    break;
                }
            }
        }
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::effective_gas_diffusivity(Table< 2, Tensor< 2, dim > > &prop_eff) const
{  
    
    prop_eff.reinit(this->gases.size(),this->gases.size());
    
    for (unsigned int i = 0; i<this->gases.size(); i++)
    {
        std::string name_i;
        name_i = this->gases[i]->name_material();
        
        for (unsigned int j = i+1; j<this->gases.size(); j++)
        {  
            std::string name_j;
            name_j = this->gases[j]->name_material();
            
            if (name_i.compare("oxygen") == 0)
            {
                
                if (name_j.compare("nitrogen") == 0)
                {
                    for (unsigned int d = 0; d<dim; d++)
                    {
                        prop_eff(i,j)[d][d] = Units::convert(D_O2[d], Units::UNIT2, Units::C_UNIT2);
                        prop_eff(j,i)[d][d] = prop_eff(i,j)[d][d];                       
                    }
                }
                else
                {
                    FcstUtilities::log << "Species "<< name_j.c_str() << " diffusivity requested in"
                    << "DummyGDL<dim>::effective_gas_diffusivity"
                    <<"not implemented"<< std::endl;
                    exit(1);
                } 
                
            }
            else if (name_i.compare("water") == 0)
            {
                
                if (name_j.compare("nitrogen") == 0)
                {
                    for (unsigned int d = 0; d<dim; d++)
                    {
                        prop_eff(i,j)[d][d] = Units::convert(D_wv[d], Units::UNIT2, Units::C_UNIT2);
                        prop_eff(j,i)[d][d] = prop_eff(i,j)[d][d];
                        break;
                    }
                }
                else
                {
                    FcstUtilities::log << "Species "<< name_j.c_str() << " diffusivity requested in"
                    << "DummyGDL<dim>::effective_gas_diffusivity"
                    <<"not implemented"<< std::endl;
                    exit(1);
                } 
            }
            else
            {
                FcstUtilities::log << "Species "<< name_i.c_str() << " diffusivity requested in"
                << "DummyGDL<dim>::effective_gas_diffusivity"
                <<"not implemented"<< std::endl;
                exit(1);
            } 
        }
    }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::effective_electron_conductivity(double& prop_eff) const
{
  if (anisotropy == false)
  {
    prop_eff = sigma_e[0];
  }
  else
  {
    FcstUtilities::log << "The member function " << __FUNCTION__
		  << " called in Class DummyGDL"
		  << "can only be used for isotropic materials. Set anisotropy to false"
		  <<std::endl;
  }
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::effective_electron_conductivity(Tensor<2,dim>& prop_eff) const
{
     for (unsigned int i=0; i<dim; i++)
       prop_eff[i][i] = this->sigma_e[i]; // Include in declare paramters and initialize.
}

//---------------------------------------------------------------------------
template <int dim>
void
NAME::DummyGDL<dim>::effective_thermal_conductivity(double& prop_eff) const
{
   if (anisotropy == false)
  {
    prop_eff = k_T[0]; //Note that in the case on isotropic parameters, it is all stored in the same location.
  }
  else
  {
    FcstUtilities::log << "The member function " << __FUNCTION__
		  << " called in Class DummyGDL"
		  << "can only be used for isotropic materials. Set anisotropy to false"
		  <<std::endl;
  }
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// Explicit instantiations. 
template class NAME::DummyGDL<deal_II_dimension>;