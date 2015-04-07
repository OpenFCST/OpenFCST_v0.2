//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: puregas_oxygen_test.cc
//    - Description: Unit testing class for GasMixture
//    - Developers: Jie Zhou and Marc Secanell
//    - $Id: GasMixture_test.cc 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

#include "GasMixture_test.h"

namespace NAME = FuelCellShop::Material;

//---------------------------------------------------------------------------
void GasMixtureTest::setup()
{
    fluid.declare_parameters(param);  
    
    param.enter_subsection("Fuel cell data");
    {
        param.enter_subsection("Materials");
        {
            
            param.enter_subsection("fluid");
            {
                
                param.set("Total pressure of gas mixture",
                          "100000.0");
                param.set("Temperature of gas mixture",
                          "293.0");  
            }
            param.leave_subsection();
        }
        param.leave_subsection();
    }
    param.leave_subsection();
    
    fluid.initialize(param);
       
}

//---------------------------------------------------------------------------
void GasMixtureTest::testChapmanEnskog_isobaric_diffusion_coefficient(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_ChapmanEnskog_isobaric_diffusion_coefficient(); 
    
    double expectedAnswer = 7.8840471941;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testChapmanEnskog_isobaric_diffusion_coefficient failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_isobaric_diffusion_coefficient(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> answer(3);
    
    fluid.get_ChapmanEnskog_isobaric_diffusion_coefficient(answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 7.8840471941;
    expectedAnswer[1] = 7.8840471941;
    expectedAnswer[2] = 7.8840471941;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficient failed!");
    }
}

//---------------------------------------------------------------------------
void GasMixtureTest::testtemp_ChapmanEnskog_isobaric_diffusion_coefficient(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_ChapmanEnskog_isobaric_diffusion_coefficient(303.0); 
    
    double expectedAnswer = 8.3512973681;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testtemp_ChapmanEnskog_isobaric_diffusion_coefficient failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testtemp_diff_ChapmanEnskog_isobaric_diffusion_coefficient(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> temperature(3);
    
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 293.0;
    
    std::vector<double> answer(3);
    
    fluid.get_ChapmanEnskog_isobaric_diffusion_coefficient(temperature,answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 7.8840471941;
    expectedAnswer[1] = 8.3512973681;
    expectedAnswer[2] = 7.8840471941;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testtemp_diff_ChapmanEnskog_isobaric_diffusion_coefficient failed!");
    }
}

///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficient. Binary gas mixture only.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(293.0); 
    
    double expectedAnswer = 0.046204176679;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperaturevector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> temperature(3);
    
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 293.0;
    
    std::vector<double> answer(3);
    
    fluid.get_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature(temperature,answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 0.046204176679;
    expectedAnswer[1] = 0.047244210264;
    expectedAnswer[2] = 0.046204176679;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperaturevector failed!");
    }
}

///@name Service functions. Chapman Enskog diffusion coefficient. Binary gas mixture only.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficient(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_ChapmanEnskog_diffusion_coefficient(); 
    
    double expectedAnswer = 7.8840471941/100000;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficientvector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> answer(3);
    
    fluid.get_ChapmanEnskog_diffusion_coefficient(answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 7.8840471941/100000;
    expectedAnswer[1] = 7.8840471941/100000;
    expectedAnswer[2] = 7.8840471941/100000;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_ChapmanEnskog_diffusion_coefficientvector failed!");
    }
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(293.0); 
    
    double expectedAnswer = 7.8840471941/100000;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressure failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressurevector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> temperature(3);
    
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 293.0;
    
    std::vector<double> answer(3);
    
    fluid.get_ChapmanEnskog_diffusion_coefficient_at_constant_pressure(temperature,answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 7.8840471941/100000;
    expectedAnswer[1] = 8.3512973681/100000;
    expectedAnswer[2] = 7.8840471941/100000;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressurevector failed!");
    }
}

//---------------------------------------------------------------------------
///@name Service functions. Derivatives of Chapman Enskog diffusion coefficient. Binary gas mixture only.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dpressure(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_DChapmanEnskog_diffusion_coefficient_Dpressure(100000.0); 
    
    double expectedAnswer = -7.8840471941/100000/100000;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testget_DChapmanEnskog_diffusion_coefficient_Dpressure failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dpressurevector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> temperature(3);
    
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 293.0;
    
    std::vector<double> pressure(3);
    
    pressure[0] = 100000.0;
    pressure[1] = 200000.0;
    pressure[2] = 100000.0;
    
    std::vector<double> answer(3);
    
    fluid.get_DChapmanEnskog_diffusion_coefficient_Dpressure(pressure,temperature,answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = -7.8840471941/100000/100000;
    expectedAnswer[1] = -8.3512973681/200000/200000;
    expectedAnswer[2] = -7.8840471941/100000/100000;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_DChapmanEnskog_diffusion_coefficient_Dpressurevector failed!");
    }
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dtemperature(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    double answer = fluid.get_DChapmanEnskog_diffusion_coefficient_Dtemperature(293.0); 
    
    double expectedAnswer = 0.046204176679/100000;
    
    TEST_ASSERT_MSG(std::fabs((answer-expectedAnswer)/expectedAnswer)<10e-5, "testget_DChapmanEnskog_diffusion_coefficient_Dtemperature failed!");
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dtemperaturevector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    fluid.set_gases(gases);
    
    std::vector<double> temperature(3);
    
    temperature[0] = 293.0;
    temperature[1] = 303.0;
    temperature[2] = 293.0;
    
    std::vector<double> pressure(3);
    
    pressure[0] = 100000.0;
    pressure[1] = 200000.0;
    pressure[2] = 100000.0;
    
    std::vector<double> answer(3);
    
    fluid.get_DChapmanEnskog_diffusion_coefficient_Dtemperature(pressure,temperature,answer); 
    
    std::vector<double> expectedAnswer (answer.size());
    
    expectedAnswer[0] = 0.046204176679/100000;
    expectedAnswer[1] = 0.047244210264/200000;
    expectedAnswer[2] = 0.046204176679/100000;
    
    for (int i=0; i<answer.size(); i++)
    {
        TEST_ASSERT_MSG(std::fabs((answer[i]-expectedAnswer[i])/expectedAnswer[i])<10e-5, "testget_DChapmanEnskog_diffusion_coefficient_Dtemperaturevector failed!");
    }
}

//---------------------------------------------------------------------------
///@name Service functions. Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.clear();
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    gases.push_back(&nitrogen);
    
    fluid.set_gases(gases);
    
    Table< 2, double > answer(gases.size(), gases.size());
    
    answer = fluid.get_ChapmanEnskog_isobaric_diffusion_coefficients();
    
    Table< 2, double > expectedAnswer(gases.size(), gases.size());
    expectedAnswer(0,1) = 7.8840471941;
    expectedAnswer(1,0) = 7.8840471941;
    expectedAnswer(0,2) = 1.9944813200;
    expectedAnswer(2,0) = 1.9944813200;
    expectedAnswer(1,2) = 7.4629860961;
    expectedAnswer(2,1) = 7.4629860961;
    
    for(unsigned int i = 0; i < gases.size(); ++i)
    {
        for(unsigned int j = 0; j < gases.size(); ++j)
        {
            if (i != j)
            {
                TEST_ASSERT_MSG(std::fabs((answer(i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_DChapmanEnskog_diffusion_coefficient_Dtemperature failed!");
            }
        }
        
    }
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.clear();
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    gases.push_back(&nitrogen);
    
    fluid.set_gases(gases);
    
    Table< 2, double > answer(gases.size(), gases.size());
    
    double temperature = 293.0;
    
    answer = fluid.get_ChapmanEnskog_isobaric_diffusion_coefficients(temperature);
    
    Table< 2, double > expectedAnswer(gases.size(), gases.size());
    expectedAnswer(0,1) = 7.8840471941;
    expectedAnswer(1,0) = 7.8840471941;
    expectedAnswer(0,2) = 1.9944813200;
    expectedAnswer(2,0) = 1.9944813200;
    expectedAnswer(1,2) = 7.4629860961;
    expectedAnswer(2,1) = 7.4629860961;
    
    for(unsigned int i = 0; i < gases.size(); ++i)
    {
        for(unsigned int j = 0; j < gases.size(); ++j)
        {
            if (i != j)
            {
                TEST_ASSERT_MSG(std::fabs((answer(i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector failed!");
            }
        }
        
    }
}

//---------------------------------------------------------------------------
///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature_table(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.clear();
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    gases.push_back(&nitrogen);
    
    fluid.set_gases(gases);
    
    std::vector< Table< 2, double > > answer (3);
    
    std::vector <double> atemperature(3);
    atemperature[0] = 293.0;
    atemperature[1] = 293.0;
    atemperature[2] = 293.0;
    
    fluid.get_DChapmanEnskog_diffusion_coefficients_Dtemperature(atemperature,answer);
    
    Table< 2, double > expectedAnswer(gases.size(), gases.size());
    expectedAnswer(0,1) = 0.046204176679/100000;
    expectedAnswer(1,0) = 0.046204176679/100000;
    expectedAnswer(0,2) = 0.012129471075/100000;
    expectedAnswer(2,0) = 0.012129471075/100000;
    expectedAnswer(1,2) = 0.043597154276/100000;
    expectedAnswer(2,1) = 0.043597154276/100000;
    
    unsigned int q = 2;
    
    for(unsigned int i = 0; i < gases.size(); ++i)
    {
        for(unsigned int j = 0; j < gases.size(); ++j)
        {
            if (i != j)
            {
                TEST_ASSERT_MSG(std::fabs((answer[q](i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector failed!");
            }
        }
        
    }
}

//---------------------------------------------------------------------------
///@name Service functions. Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficients_table(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.clear();
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    gases.push_back(&nitrogen);
    
    fluid.set_gases(gases);
    
    std::vector< Table< 2, double > > answer (3);
    
    fluid.get_ChapmanEnskog_diffusion_coefficients(answer);
    
    Table< 2, double > expectedAnswer(gases.size(), gases.size());
    expectedAnswer(0,1) = 7.8840471941/100000;
    expectedAnswer(1,0) = 7.8840471941/100000;
    expectedAnswer(0,2) = 1.9944813200/100000;
    expectedAnswer(2,0) = 1.9944813200/100000;
    expectedAnswer(1,2) = 7.4629860961/100000;
    expectedAnswer(2,1) = 7.4629860961/100000;
    
    
    unsigned int q = 2;
    
    for(unsigned int i = 0; i < gases.size(); ++i)
    {
        for(unsigned int j = 0; j < gases.size(); ++j)
        {
            if (i != j)
            {
                TEST_ASSERT_MSG(std::fabs((answer[q](i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector failed!");
            }
        }
        
    }
}

//---------------------------------------------------------------------------
void GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficients_table_tp(){
	
	std::vector< NAME::PureGas* > gases;
	
	gases.clear();
	
	gases.push_back(&oxygen);
		
	gases.push_back(&hydrogen);
	
	gases.push_back(&nitrogen);
	
	fluid.set_gases(gases);
	
	std::vector< Table< 2, double > > answer (3);
	
	std::vector <double> tpressure(3);
	tpressure[0] = 500000.0;
	tpressure[1] = 500000.0;
	tpressure[2] = 500000.0;
	
	std::vector <double> atemperature(3);
	atemperature[0] = 500.0;
	atemperature[1] = 500.0;
	atemperature[2] = 500.0;
	
    fluid.get_ChapmanEnskog_diffusion_coefficients(tpressure,atemperature,answer);
	
 	Table< 2, double > expectedAnswer(gases.size(), gases.size());
	expectedAnswer(0,1) = 19.52273379/500000;
	expectedAnswer(1,0) = 19.52273379/500000;
	expectedAnswer(0,2) = 5.0620398945/500000;
	expectedAnswer(2,0) = 5.0620398945/500000;
	expectedAnswer(1,2) = 18.4376329159709/500000;
	expectedAnswer(2,1) = 18.4376329159709/500000;
 	
	unsigned int q = 2;
	
  	for(unsigned int i = 0; i < gases.size(); ++i)
  	{
  		for(unsigned int j = 0; j < gases.size(); ++j)
  		{
  			if (i != j)
  			{
  			TEST_ASSERT_MSG(std::fabs((answer[q](i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector failed!");
			
  			}
  		}
  		
  	}
}

//---------------------------------------------------------------------------
///@name Service functions. Derivatives of Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.
//---------------------------------------------------------------------------
void GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficients_Dtemperature(){
    
    std::vector< NAME::PureGas* > gases;
    
    gases.clear();
    
    gases.push_back(&oxygen);
    
    gases.push_back(&hydrogen);
    
    gases.push_back(&nitrogen);
    
    fluid.set_gases(gases);
    
    std::vector< Table< 2, double > > answer (3);
    
    std::vector <double> tpressure(3);
	tpressure[0] = 500000.0;
	tpressure[1] = 500000.0;
	tpressure[2] = 500000.0;
    
    std::vector <double> atemperature(3);
	atemperature[0] = 500.0;
	atemperature[1] = 500.0;
	atemperature[2] = 500.0;
    
    fluid.get_DChapmanEnskog_diffusion_coefficients_Dtemperature(tpressure,atemperature,answer);
    
    Table< 2, double > expectedAnswer(gases.size(), gases.size());
 	expectedAnswer(0,1) = 0.0655584367221128/500000;
 	expectedAnswer(1,0) = 0.0655584367221128/500000;
 	expectedAnswer(0,2) = 0.0173378420615506/500000;
 	expectedAnswer(2,0) = 0.0173378420615506/500000;
 	expectedAnswer(1,2) = 0.0617900564100835/500000;
 	expectedAnswer(2,1) = 0.0617900564100835/500000;
    
    unsigned int q = 2;
    
    for(unsigned int i = 0; i < gases.size(); ++i)
    {
        for(unsigned int j = 0; j < gases.size(); ++j)
        {
            if (i != j)
            {
                TEST_ASSERT_MSG(std::fabs((answer[q](i,j)-expectedAnswer(i,j))/expectedAnswer(i,j))<10e-5, "testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector failed!");
                
            }
        }
        
    }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------