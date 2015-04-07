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
//    - $Id: porous_layer_test.h 2605 2014-08-15 03:36:44Z secanell $ 
//
//---------------------------------------------------------------------------
#ifndef _Porous_Layer_Test
#define _Porous_Layer_Test

#include <cpptest.h>
#include <boost/lexical_cast.hpp>

#include <layers/porous_layer.h>
#include<layers/design_fibrous_GDL.h>
#include <catalyst_layer.h>
/**
 * A unit test class that tests the FCST enumeration provided in system managment.
 *
 */
class PorousLayerTest: public Test::Suite
{
    public:
    PorousLayerTest()
    {
       //Add a number of tests that will be called during Test::Suite.run()
       //Generic cases
        TEST_ADD(PorousLayerTest::Knudsen_diffusion);
        TEST_ADD(PorousLayerTest::Knudsen_diffusion_derivatives);

    }
    protected:
        virtual void setup(); // setup resources... called before Test::Suite.run() ..not implemented for this test suite
        virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
        
        
    private:
        
        void Knudsen_diffusion();
        void Knudsen_diffusion_derivatives();
        
        boost::shared_ptr< FuelCellShop::Layer::GasDiffusionLayer< dim > > layer;
       
        
};

#endif
