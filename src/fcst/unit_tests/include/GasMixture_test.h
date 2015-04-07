//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: puregas_oxygen_test.h
//    - Description: Unit testing class for GasMixture
//    - Developers: Jie Zhou and Marc Secanell
//    - Id: $Id: GasMixture_test.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------

/**
 * A unit test class that tests the GasMixture class. 
 * 
 */

#ifndef _FCST_GasMixtureTest_TESTSUITE
#define _FCST_GasMixtureTest_TESTSUITE

// STD:
#include <cpptest.h>
#include <string.h>
#include <stdexcept>
#include <vector>
#include <iostream>

// openFCST
#include "GasMixture.h"

namespace NAME = FuelCellShop::Material;

class GasMixtureTest: public Test::Suite
{
    
    
public:
    
    GasMixtureTest()
    :
    fluid("fluid")
    {
        
        TEST_ADD(GasMixtureTest::testChapmanEnskog_isobaric_diffusion_coefficient);
        TEST_ADD(GasMixtureTest::testget_ChapmanEnskog_isobaric_diffusion_coefficient);
        TEST_ADD(GasMixtureTest::testtemp_ChapmanEnskog_isobaric_diffusion_coefficient);
        TEST_ADD(GasMixtureTest::testtemp_diff_ChapmanEnskog_isobaric_diffusion_coefficient);
        
        ///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficient. Binary gas mixture only.
        
        TEST_ADD(GasMixtureTest::testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature);
        TEST_ADD(GasMixtureTest::testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperaturevector);
        
        ///@name Service functions. Chapman Enskog diffusion coefficient. Binary gas mixture only.
        
        TEST_ADD(GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficient);
        TEST_ADD(GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficientvector);
        TEST_ADD(GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressure);
        TEST_ADD(GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressurevector);
        
        ///@name Service functions. Derivatives of Chapman Enskog diffusion coefficient. Binary gas mixture only.
        
        TEST_ADD(GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dpressure);
        TEST_ADD(GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dpressurevector);
        TEST_ADD(GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dtemperature);
        TEST_ADD(GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficient_Dtemperaturevector);
        
        ///@name Service functions. Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
        
        TEST_ADD(GasMixtureTest::testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table);
        TEST_ADD(GasMixtureTest::testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector);
        
        ///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
        
        TEST_ADD(GasMixtureTest::testget_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature_table);
        
        ///@name Service functions. Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.
        
        TEST_ADD(GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficients_table);
        TEST_ADD(GasMixtureTest::testget_ChapmanEnskog_diffusion_coefficients_table_tp);
        
        ///@name Service functions. Derivatives of Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.
        
        TEST_ADD(GasMixtureTest::testget_DChapmanEnskog_diffusion_coefficients_Dtemperature);
        
        
    }
    
    
protected:
    virtual void setup(); // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:
    
    FuelCellShop::Material::Oxygen oxygen;
    
    FuelCellShop::Material::Hydrogen hydrogen;
    
    FuelCellShop::Material::Nitrogen nitrogen;
    
    FuelCellShop::Material::GasMixture fluid;
    
    ParameterHandler param;
    
    void testChapmanEnskog_isobaric_diffusion_coefficient();
    void testget_ChapmanEnskog_isobaric_diffusion_coefficient();
    void testtemp_ChapmanEnskog_isobaric_diffusion_coefficient();
    void testtemp_diff_ChapmanEnskog_isobaric_diffusion_coefficient();
    
    ///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficient. Binary gas mixture only.
    
    void testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperature();
    void testget_DChapmanEnskog_isobaric_diffusion_coefficient_Dtemperaturevector();
    
    ///@name Service functions. Chapman Enskog diffusion coefficient. Binary gas mixture only.
    
    void testget_ChapmanEnskog_diffusion_coefficient();
    void testget_ChapmanEnskog_diffusion_coefficientvector();
    void testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressure();
    void testget_ChapmanEnskog_diffusion_coefficient_at_constant_pressurevector();
    
    ///@name Service functions. Derivatives of Chapman Enskog diffusion coefficient. Binary gas mixture only.
    
    void testget_DChapmanEnskog_diffusion_coefficient_Dpressure();
    void testget_DChapmanEnskog_diffusion_coefficient_Dpressurevector();
    void testget_DChapmanEnskog_diffusion_coefficient_Dtemperature();
    void testget_DChapmanEnskog_diffusion_coefficient_Dtemperaturevector();
    
    ///@name Service functions. Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
    
    void testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table();
    void testget_ChapmanEnskog_isobaric_diffusion_coefficients_for_table_vector();
    
    ///@name Service functions. Derivatives of Chapman Enskog isobaric diffusion coefficients. Ternary and more complicated gas mixtures.
    
    void testget_DChapmanEnskog_isobaric_diffusion_coefficients_Dtemperature_table();
    
    ///@name Service functions. Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.
    
    void testget_ChapmanEnskog_diffusion_coefficients_table();
    void testget_ChapmanEnskog_diffusion_coefficients_table_tp();
    
    ///@name Service functions. Derivatives of Chapman Enskog diffusion coefficients. Ternary and more complicated gas mixtures.	  
    
    void testget_DChapmanEnskog_diffusion_coefficients_Dtemperature();
    
    
};

#endif
