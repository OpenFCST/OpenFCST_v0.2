//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PSD_test.cc
//    - Description: Unit testing class for PSD
//    - Developers: Prafful Mangal
//    - Id: $Id: PSD_HI_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include "PSD_HI_test.h"
#include <boost/concept_check.hpp>

//-------------------------------------------------------------
void PSD_HI_Test::setup()
{
    ParameterHandler param;
    
    boost::shared_ptr<FuelCellShop::Layer::GasDiffusionLayer<dim> > CGDL;
    
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Cathode gas diffusion layer", param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Cathode gas diffusion layer");
        {
            
            psd_object.declare_parameters(param);
            
            param.enter_subsection("PSD parameters");
            {
                param.enter_subsection("BasePSD");
                {
                    
                    //param.set("porosity", "0.84");
                    param.set("Gamma", "0.0728");
                    param.set("Contact angle", "1.396");
                    param.set("lambda", "1.0");
                    param.set("Volume fraction Hydrophilic", "0.7");
                    param.set("probability P_b", "1.0");
                    param.set("Mode probability global", "0.72, 0.28");
                    param.set("Mode characteristic radius global", "34.0, 14.2");
                    param.set("Mode width global", "0.35, 1.0");
                    param.set("psd type", "HIPSD"); 
                    
                    param.enter_subsection("HIPSD");
                    {
                        param.set("Hydrophilic Mode probability global", "0.72, 0.28");
                        param.set("Hydrophilic Mode characteristic radius global", "34.0, 14.2");
                        param.set("Hydrophilic Mode width global", "0.35, 1.0");
                        param.set("capillay pressure", "10100.0");   
                    }
                    param.leave_subsection(); 
                    
                }
                param.leave_subsection();
            }
            param.leave_subsection();
            
            param.enter_subsection("Generic data");
            {
                param.set("Porosity", "0.84");
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Cathode gas diffusion layer");
        {
            psd_object.initialize(param);
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    CGDL = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Cathode gas diffusion layer",param);
    
    psd_object.set_porosity (CGDL->get_porosity());
    
}

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_rc_HI()
{
    
    std::vector<double> answer(0.0);
    
    psd_object.set_critical_radius();
    
    psd_object.get_critical_radius(answer);
    
    double expectedAnswer = 2.507024;
  
    std::ostringstream streamOut;
    streamOut <<"The value of the rc_HI (microns) is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
    TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-6, streamOut.str().c_str()); 
}

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_k_sat_HI()
{
    double answer(0.0);
    
    psd_object.set_saturation();
  
    psd_object.get_global_saturated_permeability(answer); 
    
    double expectedAnswer = 58.1327578;
  
    std::ostringstream streamOut;
    streamOut <<"The value of the k_sat_HI (microns^2) is: "<<answer<<". The expected value is: "<<expectedAnswer<<std::endl;
    TEST_ASSERT_DELTA_MSG(expectedAnswer, answer, 1e-7, streamOut.str().c_str());
}

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_sat_HI()
{
    std::vector<double> answer(0.0);
  
    psd_object.get_saturation(answer); 
    
    double expectedAnswer = 0.0081234135;
  
    std::ostringstream streamOut;
    streamOut <<"The value of the Saturation_HI is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
    TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-7, streamOut.str().c_str());
}

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_k_L_HI()
{
    std::vector<double> answer(0.0);
    psd_object.get_pore_HI_liquid_saturated_permeability(answer);  
    
    double expectedAnswer = 2.9317742175e-9;
  
    std::ostringstream streamOut;
    streamOut <<"The value of the K_L_HI(microns^2) is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
    TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-13, streamOut.str().c_str());
}

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_kr_L_HI()
{
    std::vector<double> answer(0.0);
    
    psd_object.get_relative_liquid_permeability(answer);  
    
    double expectedAnswer = 5.04323951277e-11;
  
    std::ostringstream streamOut;
    streamOut <<"The value of the Kr_L_HI is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
    TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-15, streamOut.str().c_str());
}

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_k_G_HI()
{
    std::vector<double> answer(0.0);
  
    psd_object.get_pore_HI_gas_saturated_permeability(answer);  
    
    double expectedAnswer = 40.0344411166;
  
    std::ostringstream streamOut;
    streamOut <<"The value of the K_G_HI(microns^2) is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
    TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-5, streamOut.str().c_str());
}

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_kr_G_HI()
{
    std::vector<double> answer(0.0);
  
    psd_object.get_relative_gas_permeability(answer);  
    
    double expectedAnswer = 0.6886726621;
  
    std::ostringstream streamOut;
    streamOut <<"The value of the Kr_G_HI is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
    TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-7, streamOut.str().c_str());
}

//-------------------------------------------------------------
void PSD_HI_Test::testcompute_interfacial_area_per_volume_HI()
{
    std::vector<double> answer(0.0);
  
    psd_object.get_liquid_gas_interfacial_surface(answer);  
    double expectedAnswer = 0.000113962;
  
    std::ostringstream streamOut;
    streamOut <<"The value of the interfacial area per unit volume for hydrophilic pores is: "<<answer[0]<<". The expected value is: "<<expectedAnswer<<std::endl;
    TEST_ASSERT_DELTA_MSG(expectedAnswer, answer[0], 1e-9, streamOut.str().c_str());
}

//-------------------------------------------------------------
//-------------------------------------------------------------