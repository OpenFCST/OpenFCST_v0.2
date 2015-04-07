#include <water_agglomerate_test.h>
using namespace FuelCellShop;


void  WaterAgglomerateTest::testO2CurrentDensity()
{
    ParameterHandler param;
    //kinetics.declare_parameters(param);
    FuelCellShop::Layer::CatalystLayer<dim>::declare_CatalystLayer_parameters("CathodeCL", param);
    FuelCellShop::Material::PolymerElectrolyteBase::declare_PolymerElectrolyte_parameters(param);
    FuelCellShop::Material::CatalystSupportBase::declare_CatalystSupport_parameters(param);
    FuelCellShop::Material::CatalystBase::declare_Catalyst_parameters(param);
    OC.declare_parameters(param);

    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("CathodeCL");
        {
            param.set("Material id", "4");
            param.set("Catalyst layer type", "MultiScaleCL");
            
            param.enter_subsection("ConventionalCL");
            {
                param.set("Platinum loading on support (%wt)", "4:0.46");
                param.set("Platinum loading per unit volume (mg/cm3)", "4:400");
                param.set("Electrolyte loading (%wt)", "4:0.30");
                param.set("Active area [cm^2/cm^3]", "4:2.0e5");
            }
            param.leave_subsection();

            param.enter_subsection("MultiScaleCL");
            {
                param.enter_subsection("MicroScale");
                {
                    param.set("Microscale type", "WaterAgglomerateNumerical");
                }
                param.leave_subsection();
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }   
    param.leave_subsection();
    
    CCL = FuelCellShop::Layer::CatalystLayer<dim>::create_CatalystLayer("CathodeCL", param);
    CCL->set_local_material_id(4);

    OC.initialize(param);

    CCL->set_reaction_kinetics(ORR);
    CCL->set_constant_solution(OC.get_pc_Pa(), total_pressure);
    CCL->set_constant_solution(OC.get_T(), temperature_of_REV);


    std::vector<FuelCellShop::SolutionVariable> sols;
    sols.push_back(FuelCellShop::SolutionVariable(0.1,1,oxygen_molar_fraction));
    sols.push_back(FuelCellShop::SolutionVariable(0.7,1,electronic_electrical_potential));
    sols.push_back(FuelCellShop::SolutionVariable(-0.1,1,protonic_electrical_potential));
    sols.push_back(FuelCellShop::SolutionVariable(8,1,membrane_water_content));
    CCL->set_solution(sols);

    std::vector<double> current(1, 0.0);
    CCL->current_density(current);

    //water numerical: 189.5 match within 5%
    double expected = 25.0;
    double match = expected*0.05;
    std::string msg = "Current from model ionomer analytical does not match expected results (" + std::to_string(std::abs(100*(expected-current[0])/expected)) + "% wrong)!";

    TEST_ASSERT_DELTA_MSG(expected, current[0], match, msg.c_str());
    //TEST_ASSERT_MSG(false, "Test not yet implemented");


}


void WaterAgglomerateTest::test02CurrentDerivative(){

    std::map< VariableNames, std::vector<double> > derivatives;
    std::vector<VariableNames> sol_names;
    sol_names.push_back(oxygen_molar_fraction);
    sol_names.push_back(protonic_electrical_potential);
    sol_names.push_back(electronic_electrical_potential);

    CCL->set_derivative_flags(sol_names);
    CCL->derivative_current_density(derivatives);
    double dCurrentdO2 = derivatives[oxygen_molar_fraction][0];    //ionomer numerical : 2062.3 //ionomer analytical : 2062.4  //water numerical: 1885.3
    double dCurrentdPhi_m = derivatives[protonic_electrical_potential][0]; //ionomer numerical : 6082.1 //ionomer analytical : 6074.5  //water numerical: 6090.75
    double dCurrentdPhi_s = derivatives[electronic_electrical_potential][0]; //ionomer numerical :-6066   //ionomer analytical : -6074.5 //water numerical: -6027.7

    double expDCurrentDO2 = 255.0;     double diffDCurrentDO2    = expDCurrentDO2*0.05;
    double expDCurrentPhi_m = 835.0;  double diffDCurrentPhi_m  = expDCurrentPhi_m*0.05;
    double expDCurrentPhi_s = -835.0;  double diffDCurrentPhi_s  = -expDCurrentPhi_s*0.05; //must be positive

    TEST_ASSERT_DELTA_MSG(dCurrentdO2, expDCurrentDO2, diffDCurrentDO2, "Current derivative w.r.t. O2 is incorrect");
    TEST_ASSERT_DELTA_MSG(dCurrentdPhi_m, expDCurrentPhi_m, diffDCurrentPhi_m, "Current derivative w.r.t. Phi_m is incorrect");
    TEST_ASSERT_DELTA_MSG(dCurrentdPhi_s, expDCurrentPhi_s, diffDCurrentPhi_s, "Current derivative w.r.t. Phi_s is incorrect");

}

