//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: PSD_base.h
//    - Description: Base class for pore size distribution model.
//    - Developers: 2009-13 by Marc Secanell, University of Alberta
//                  2013-14 by Jie Zhou, University of Alberta
//    - $ $
//
//---------------------------------------------------------------------------
#ifndef _FUELCELLSHOP__NONE__PSD_H
#define _FUELCELLSHOP__NONE__PSD_H

// Include deal.II classes
#include <base/parameter_handler.h>
#include <base/point.h> 
#include <base/function.h>
#include <lac/vector.h>
#include <fe/fe_values.h>  

//Include STL
#include <cmath>
#include <iostream>

// Include OpenFCST routines:
#include "application_core/fcst_variables.h"
#include "application_core/system_management.h"
#include "utils/fcst_utilities.h"
#include "utils/fcst_constants.h"
#include "PSD_base.h"

using namespace dealii;

namespace FuelCellShop
{
    
    
    namespace MicroScale
    {
    /**
     *  @brief This is a class for the layer does not contain PSD information
     * 
     * 
     * @author J. Zhou
     * 
     *   Marc Secanell
     * @date 2014
     */
        template <int dim>
        class NonePSD 
        : 
        public BasePSD<dim>
        {
        public:
            ///@name Initalization
            //@{
            /** 
             * Consturctor
             */
            NonePSD();
                                    
            /**
             * Constructor.
             */
            NonePSD (std::string name);

            /**
             * Destructor.
             */
            ~NonePSD() {}

            /**
             * Concrete name used for objects of this class. This name is used when
             * setting up the subsection where the data is stored in the input file.
             * 
             * The data will be store under
             * \code
             *    set psd type = NonePSD # <-here I select the type of object of type psd
             *    subsection NonePSD # <- this is the concrete_name for this class
             *       set all info relevant to this object
             *    end
             * end
             * \endcode
             */
            static const std::string concrete_name;
            
            //@}
            
            ///@name Accessors and info
            //@{
            /**
             * This function is used to compute the saturation by using PSD
             * 
             *  This function contains no information, does not return anything.
             */
            
            virtual inline void get_saturation(std::vector<double>& S) const {}
            
            /**
             * This function is used to compute the saturated_permeability by using PSD
             * 
             * This function contains no information, does not return anything.
             */
            
            virtual inline void get_global_saturated_permeability(double& saturated_permeability) const {}
            
            /**
             * This function is used to compute the liquid_permeability by using PSD
             * 
             * This function contains no information, does not return anything.
             */
            
            virtual inline void get_relative_liquid_permeability(std::vector<double>& liquid_permeability) const {}
            
            /**
             * This function is used to compute the gas_permeability by using PSD
             * 
             * This function contains no information, does not return anything.
             */
            
            virtual inline void get_relative_gas_permeability(std::vector<double>& gas_permeability) const {}
            
            /**
             * This function is used to compute the liquid_gas_interfacial_surface by using PSD
             * 
             * This function contains no information, does not return anything.
             */
                                                                            
            virtual inline void get_liquid_gas_interfacial_surface(std::vector<double>& HI_liquid_gas_interfacial_surface) const {}
            
            /**
             * This function is used to compute the pore_wetted_wall by using PSD
             * 
             * This function contains no information, does not return anything.
             */
            
            virtual inline void get_wetted_wall_surface_area(std::vector<double>& wetted_wall_surface_area) const {}
            
            
            /**
             * This function is used to compute the knudsen_radius by using PSD
             * 
             * This function contains no information, does not return anything.
             */
            
            virtual inline void get_knudsen_radius(std::vector<double>& knudsen_radius) const {}
            
            /**
             * This function is used to compute the diffusivity by using PSD
             * 
             * This function contains no information, does not return anything.
             */
            virtual inline void get_diffusivity() const {}
            //@}
            
        protected:
            
            ///@name Instance Delivery (Replica creator)
            //@{
            /**
             * This member function is used to create an object of type psd
             */
            virtual boost::shared_ptr<FuelCellShop::MicroScale::BasePSD <dim>> create_replica (const std::string &psd_section_name)
            {
                return boost::shared_ptr<FuelCellShop::MicroScale::BasePSD <dim>> (new FuelCellShop::MicroScale::NonePSD<dim> (psd_section_name));
            }   
            //@}           
            ///@name Instance Delivery (Prototype)
            //@{
            /**
             * PROTOTYPE is the pointer is the dynamic pointer pointing to the NonePSD class itself.
             */            
            static NonePSD<dim> const* PROTOTYPE;
            //@} 
        };
        
    } // PSD
    
}  // FuelCellShop

#endif