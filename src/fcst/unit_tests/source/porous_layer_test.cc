//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: 
//    - Description: 
//    - Developers: 
//    - $Id: porous_layer_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include <porous_layer_test.h>
#include <boost/config/no_tr1/utility.hpp>

//---------------------------------------------------------------------------
void PorousLayerTest::setup(){
    
    ParameterHandler param;
    // Note since I cannot create a parent, I will test the functions of the parent in one of the children:
    FuelCellShop::Layer::GasDiffusionLayer<dim>::declare_GasDiffusionLayer_parameters("Sample Gas Diffusion Layer", param);
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Sample Gas Diffusion Layer");
        {
            param.set("Gas diffusion layer type", "DesignFibrousGDL");
            {
                param.enter_subsection("Generic data");
                { 
                    param.set("Knudsen pore radius [um]",5e-8);
                }
                param.leave_subsection();
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    
    layer = FuelCellShop::Layer::GasDiffusionLayer<dim>::create_GasDiffusionLayer("Sample Gas Diffusion Layer", param);

}

//---------------------------------------------------------------------------
  /*
void PorousLayerTest::molecular_diffusion_test(){
    

  double pressure(1.0); //atm
  double temperature(273); //K
  FuelCellShop::Layer::PorousLayer<dim> test("test");
  FuelCellShop::Material::Oxygen oxygen;
  FuelCellShop::Material::Nitrogen nitrogen;
  FcstUtilities::log<<"Computing diffusion coefficient of: "<<oxygen.name_material()<<" in "<<nitrogen.name_material()<<" at "<<pressure<<" atm and "<<temperature<<" Kelvin"<<std::endl;
  std::vector<FuelCellShop::Material::PureGas*> gases;
  gases.push_back(&oxygen);
  gases.push_back(&nitrogen);
  std::vector<VariableNames> deriv_flags;
  deriv_flags.push_back(total_pressure);
  deriv_flags.push_back(temperature_of_REV);
  test.set_derivative_flags(deriv_flags);
  test.set_gases_and_compute(gases, pressure, temperature);
  Table<2, double> D_eff;
  test.gas_diffusion_coefficients(D_eff);
  std::vector<Table< 2, double > > dD_dx;
  test.derivative_gas_diffusion_coefficients(dD_dx);

  FcstUtilities::log<<"The diffusion coefficients are: D_11="<<D_eff(0,0)<<" D_12 = "<<D_eff(0,1)<<"m2/s"<<std::endl;
  FcstUtilities::log<<"The diffusion coefficients are: D_21="<<D_eff(1,0)<<" D_22 = "<<D_eff(1,1)<<"m2/s"<<std::endl;
  FcstUtilities::log<<"The experimental value at 273K and 1atm is: 0.181 cm2/s"<<std::endl;

  FcstUtilities::log<<"The diffusion coefficients are: dD_11_dp="<<dD_dx[0](0,0)<<" dD_12_dp = "<<dD_dx[0](0,1)<<"m2/s"<<std::endl;
  FcstUtilities::log<<"The diffusion coefficients are: dD_21_dp="<<dD_dx[0](1,0)<<" dD_22_dp = "<<dD_dx[0](1,1)<<"m2/s"<<std::endl;

  FcstUtilities::log<<"The diffusion coefficients are: dD_11_dT="<<dD_dx[1](0,0)<<" dD_12_dT = "<<dD_dx[1](0,1)<<"m2/s"<<std::endl;
  FcstUtilities::log<<"The diffusion coefficients are: dD_21_dT="<<dD_dx[1](1,0)<<" dD_22_dT = "<<dD_dx[1](1,1)<<"m2/s"<<std::endl;

}
 */  
//---------------------------------------------------------------------------
void PorousLayerTest::Knudsen_diffusion(){
    
    std::vector<double> T_dummy(8, 298.0);
    FuelCellShop::SolutionVariable sols(T_dummy,temperature_of_REV);
    
    double D_k_check(1.480150e-5);
    
    std::vector<double> D_k;
    FuelCellShop::Material::Oxygen oxygen;
    
    layer->Knudsen_diffusion(&oxygen, sols, D_k);
    
    //Final check
    std::ostringstream streamOut;
    
    for (unsigned int i=0; i<T_dummy.size(); ++i){
        streamOut <<"The value of the D_k is: "<<D_k[i]<<". The expected value is: "<<D_k_check<<std::endl;
        streamOut <<"Knudsen_diffusion!"<<std::endl;
        TEST_ASSERT_MSG(std::fabs((D_k_check - D_k[i])/D_k_check)<1e-6, streamOut.str().c_str());
        streamOut.flush();
    }
}

//---------------------------------------------------------------------------
void PorousLayerTest::Knudsen_diffusion_derivatives(){
    
    double T(298.0);
    std::vector<double> T_dummy(8, T);
    FuelCellShop::SolutionVariable sols(T_dummy,temperature_of_REV);
    
    double D_k_check(1.480150e-5);
    double dD_k_dT_check(0.5*D_k_check*(1/T));
    
    std::vector<double> D_k;
    std::vector<double> dD_k_dT;
    FuelCellShop::Material::Oxygen oxygen;
    
    layer->Knudsen_diffusion(&oxygen, sols, D_k, dD_k_dT);
    
    //Final check
    std::ostringstream streamOut;
    streamOut <<"The value of the D_k is: "<<D_k[0]<<". The expected value is: "<<D_k_check<<std::endl;
    streamOut <<"Knudsen_diffusion_derivatives!"<<std::endl;
    TEST_ASSERT_MSG(std::fabs((D_k_check - D_k[0])/D_k_check)<1e-6, streamOut.str().c_str());    
    streamOut.flush();
    
    streamOut <<"The value of the dD_k_dT is: "<<dD_k_dT[0]<<". The expected value is: "<<dD_k_dT_check<<std::endl;
    streamOut <<"Knudsen_diffusion_derivatives!"<<std::endl;
    TEST_ASSERT_MSG(std::fabs((dD_k_dT_check - dD_k_dT[0])/dD_k_dT_check)<1e-6, streamOut.str().c_str());
    streamOut.flush();
    
}
