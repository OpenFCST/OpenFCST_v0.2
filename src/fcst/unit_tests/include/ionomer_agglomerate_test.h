/**
 * A unit test class that tests the fcst_units class. 
 * 
 */



#ifndef _FCST_IonomerAgglomerate_TESTSUITE
#define _FCST_IonomerAgglomerate_TESTSUITE

#include <map>
#include <base_layer.h>
#include <cpptest.h>
#include <boost/lexical_cast.hpp>
#include <string>
#include <catalyst_layer.h>
#include <tafel_kinetics.h>
#include <operating_conditions.h>
#include "boost/shared_ptr.hpp"
#include <chrono>
#include <fcst_utilities.h>

class IonomerAgglomerateTest: public Test::Suite
{
public:
	IonomerAgglomerateTest()
    {

        TEST_ADD(IonomerAgglomerateTest::testO2CurrentDensity);
        TEST_ADD(IonomerAgglomerateTest::testO2CurrentDerivative);
        //TEST_ADD(IonomerAgglomerateTest::testH2CurrentDensity);
    }
protected:
    virtual void setup() {} // setup resources... called before Test::Suite.run() ..not implemented for this test suite
    virtual void tear_down() {} // remove resources...called after Test::Suite.run()  ..not implemented for this test suite
private:


    void testO2CurrentDensity();
    void testO2CurrentDerivative();
    void testH2CurrentDensity();
    //Layer
    boost::shared_ptr<FuelCellShop::Layer::CatalystLayer<dim> > CCL;
    FuelCell::OperatingConditions OC;

};

#endif
