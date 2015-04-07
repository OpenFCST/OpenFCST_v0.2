//---------------------------------------------------------------------------
//    $Id: operating_conditions.h 2605 2014-08-15 03:36:44Z secanell $
//
//    Copyright (C) 2009 by Marc Secanell
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef _FUELCELL_OPERATING_CONDITIONS__H
#define _FUELCELL_OPERATING_CONDITIONS__H

// Include deal.II classes
#include <base/parameter_handler.h>


#include <boost/lexical_cast.hpp>

#include <logging.h>
#include "fcst_constants.h"
#include <application_core/initial_and_boundary_data.h>
#include <grid/geometry.h>

using namespace dealii;

namespace FuelCell
{
    
    /**
     * Class used to store, read from file and define the operating conditions for a fuel cell.
     * It stores the following:
     * - cell temperature
     * - cell voltage
     * - oxygen mole fraction at the inlet of the channel
     * - relative humidity at anode and cathode
     * - total pressure at anode and cathode
     * 
     * This information, in conjunction with the water saturation equation is used to compute the
     * mole fractions for each component in the mixture.
     * 
     * If the inlet is humidified air, please set the oxygen mole fraction to 0.21 (default).
     * 
     * In the input file, the following parameters can be specified (see declare_parameters ):
     * @code
     * subsection Fuel cell data
     * (...)
     *   subsection Operating conditions
     *     set Adjust initial solution and boundary conditions = false
     *     set Temperature cell = 353
     *     set Cathode pressure = 101325
     *     set Cathode initial oxygen mole fraction (prior to humidification) = 0.21
     *     set Cathode relative humidity = 0.7
     *     set Anode pressure = 101325
     *     set Anode relative humidity = 0.7
     *     set Voltage cell = 0.6
     *     set Voltage drop in the anode = 0.015
     *     set Open circuit voltage = 1.23
     *   end
     * end
     * @endcode 
     * 
     * <h3>Usage Details:</h3>
     * 
     * In order to use this class, first an object of the class needs to be created. Usually, one such objects exists in every application. To create the object,
     * include the .h file in the include application file and in the application data member section add the object. For example:
     * @code
     * #include "operating_conditions.h"
     * 
     * // Then in the data member declaration (usually a private member)
     * FuelCell::OperatingConditions OC;
     * @endcode
     * 
     * Once the object is created, the section where the input data will be specified in the input file needs to be delcared. To do so, in the declare_parameters section 
     * of your application call the following:
     * 
     * @code
     * //--------- IN DECLARE_PARAMETERS ------------------------------------------------------
     * template <int dim>
     * void 
     * NAME::AppCathode<dim>::declare_parameters(ParameterHandler& param)
     * {
     *   (...)
     *   OC.declare_parameters(param);
     *   (...)
     * }
     * @endcode
     *          
     * 
     * Finally, once the input file has been read by our application, your class needs to be initialized. This is achieved using the function initialize()
     * @code
     * //--------- IN INITIALIZE ------------------------------------------------------
     * template <int dim>
     * void
     * NAME::AppCathode<dim>::_initialize(ParameterHandler& param)
     * {   
     *  (...) 
     *  OC.initialize(param);
     * }
     * @endcode
     * 
     * You are now ready to use your OperatingConditions object!.
     * 
     * 
     * @author M. Secanell, 2009-2013
     * 
     */
    class OperatingConditions
    {
    public:
        /**
         * Constructor
         */
        OperatingConditions();
        
        /**
         * Destructor
         */
        ~OperatingConditions();
        
        /**
         * Declare all necessary parameters in order to compute the coefficients
         * 
         * The parameters that can be specified in the input file are as follows:
         * 
         * @code
         * subsection Fuel cell data
         * (...)
         *   subsection Operating conditions
         *     set Adjust initial solution and boundary conditions = false  
         *     set Temperature cell = 353
         *     set Cathode pressure = 101325
         *     set Cathode initial oxygen mole fraction (prior to humidification) = 0.21
         *     set Cathode relative humidity = 0.7
         *     set Anode pressure = 101325
         *     set Anode relative humidity = 0.7
         *     set Voltage cell = 0.6
         *     set Voltage drop in the anode = 0.015
         *     set Open circuit voltage = 1.23
         *   end
         * end
         * @endcode 
         */
        void declare_parameters (ParameterHandler &param) const;

        /**
         * Class used to read in data and initialize the necessary data
         * to compute the coefficients.
         */
        void initialize (ParameterHandler& param);
        
        
        /**
         * Get the water vapour saturation pressure in atmospheres (atm) using cell temperature. 
         */
        double saturation_pressure() const;
        
        /**
         * Get the mole fraction of water at the cathode channel from pressure, temperature, and relative humidity
         *   note that the anode mole fraction of water is not required, as there are only two species
         *   so x_wv_a = 1 - x_h2 
         */
        double get_x_wv() const;
        
        /**
         * Get the mole fraction of oxygen at the cathode channel from pressure, temperature, and relative humidity
         */
        double get_x_o2() const;
        
        /**
         * Get the mole fraction of hydrogen at the anode channel from pressure, temperature, and relative humidity
         */       
        double get_x_h2() const;
        
        /**
         * NOTE: This function is redefined in base_kinetics class, considering variable temperature and gas pressures. The function here should only be used for defining Initial condition values.
         * Get the theoretical cell voltage using the cell temperature, and the reactant gas pressures
         */
        double voltage_cell_th();
        
        /**
         * Output operating conditions to screen
         */
        void print_operating_conditions() const;
        
        /**
         * 
         */
        void adjust_initial_solution(std::vector< component_materialID_value_map >& maps,
                                     const boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid) const;
        /**
         * Routine used in order to adjust boundary conditions for O2 and cell voltage. This routine should be called after
         * the component_boundaryID_value_map has been initialized.
         * 
         * In order for the application to work well, the boundary IDs for the channel/GDL and land/GDL interfaces
         * need to be properly identified in the grid section of the input file, i.e.
         * @code
         * subsection Grid generation
         *   subsection Internal mesh generator parameters
         *     subsection Boundary ID
         *       set c_Ch/GDL = 100
         *       set c_BPP/GDL = 200
         *       set a_Ch/GDL = 100
         *       set a_BPP/GDL = 200
         *     end  
         *   end
         * end
         * @endcode
         */
        void adjust_boundary_conditions(std::vector< component_boundaryID_value_map >& maps,
                                        const boost::shared_ptr< FuelCellShop::Geometry::GridBase<dim> > grid) const;
        /**
         * Get the total gas concentration in the cathode. 
         * NOTE: This is a constant value. Use only if the total gas concentration
         * is assumed to be constant
         */
        inline double get_c_c() const
        {	return c_c; }
        /** Return cell temperture as input in Operating Conditions subsection */
        inline double get_T() const
        { return T_cell; }
        /** Return cell voltage as input in Operating Conditions subsection */
        inline double get_V() const
        { return V_cell; }
        /** Return the voltage drop in the anode 
         * 
         * \note Can be used as boundary condition for anode model or initial condition
         * 		in a full MEA model
         */
        inline double get_dV_a() const
        { return dV_a; }
        /** Return cathode pressure as input in Operating Conditions subsection */
        inline double get_pc_Pa() const
        { return p_c; }
        /** Return cathode pressure as input in Operating Conditions subsection */
        inline double get_pc_atm() const
        { return p_c/101325; }
        /** Get the total gas concentration in the anode. */
        inline double get_c_a() const
        {	return c_a; }
        /** Return anode pressure as input in Operating Conditions subsection */
        inline double get_pa_Pa() const
        { return p_a; }
        /** Return anode pressure as input in Operating Conditions subsection */
        inline double get_pa_atm() const
        { return p_a/101325; }
        /** Return anode relative humidity as input in Operating Conditions subsection */
        inline double get_RH_a() const
        {return RH_a;}
        /** Return cathode relative humidity as input in Operating Conditions subsection */
        inline double get_RH_c() const
        {return RH_c;}
        /** Get the open circuit voltage for the cell */
        inline double get_OCV() const
        {return OCV;}
        
    private:
        //------------ BOUNDARY CONDITIONS -------------------------------
        /** Bool set to true if you want to modify boundary conditions */
        bool adjust_BC;
        //------------ CONSTANTS -----------------------------------------
        double R; //Gas constant 8.3144 J / (mol K)
        
        /** Initial amount of oxygen in channel prior to humidification */
        double channel_oxygen_mole_fraction;
        
        //------------ CELL DATA -----------------------------------------
        /** Operating temperature of the cell */
        double T_cell;
        /** Operating voltage of the cell */
        double V_cell;
        /** Voltage drop in the anode */
        double dV_a;
        /** Open circuit voltage for the cell */
        double OCV;
        /** Theoretical voltage for the cell */
        double E_th;
        
        //------------- ANODE DATA -------------------------------------
        /** Pressure of the gas mixture in the anode B.C. */
        double p_a; 
        /** Concentration of the gas mixture in the anode B.C. */
        double c_a;
        /** Relative humidity of the gas mixture in the anode B.C */
        double RH_a;
        
        //------------- CATHODE DATA -------------------------------------
        /** Pressure of the gas mixture in the cathode B.C. */
        double p_c;
        /** Concentration of the gas mixture in the cathode B.C. */
        double c_c;
        /** Relative humidity of the gas mixture in the anode B.C */
        double RH_c;
        
    };
    
}

#endif // _FUELCELL_OPERATING_CONDITIONS__H
