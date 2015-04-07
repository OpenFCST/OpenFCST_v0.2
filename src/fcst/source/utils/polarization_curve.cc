//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: polarization_curve.cc
//    - Description: Child of ApplicationWrapper used to implement IV curve analysis
//    - Developers: M. Secanell
//    - $Id: polarization_curve.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include "utils/polarization_curve.h"

//---------------------------------------------------------------------------
template <int dim>
FuelCell::PolarizationCurve<dim>::PolarizationCurve()
{}

//---------------------------------------------------------------------------
template <int dim>
FuelCell::PolarizationCurve<dim>::~PolarizationCurve()
{}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::PolarizationCurve<dim>::declare_parameters(ParameterHandler& param) const
{
    // Initialize to parameters needed for evaluate (note: this is provisional)
    param.enter_subsection("Simulator");
    {
        param.enter_subsection("Polarization Curve");
        {
            param.declare_entry("Polarization curve file output",
                                "polarization_curve.dat",
                                Patterns::Anything(),
                                "File where the polarization curve results should be stored");
            param.declare_entry ("Initial voltage [V]",
                                 "1",
                                 Patterns::Double(),
                                 "Voltage at which the first point in the polarization curve will be evaluated");
            param.declare_entry ("Final voltage [V]",
                                 "0.1",
                                 Patterns::Double(),
                                 "Voltage at which the polarization curve will be terminated");
            param.declare_entry ("Increment [V]",
                                 "0.1",
                                 Patterns::Double(),
                                 "Spacing between points in the polarization curve");
            param.declare_entry ("Adaptive Increment",
                                 "true",
                                 Patterns::Bool(),
                                "Set to true if you would like to reduce the voltage increment adaptively"
                                "if convergence could not be achieved with the larger value");
            param.declare_entry ("Min. Increment [V]",
                                 "0.025",
                                 Patterns::Double(),
                                 "If Adaptive Increment? is set to true, this value controls the "
                                 "minimum change in cell voltage before the polarization curve gives up"
                                 "and the voltage is updated again. Note that this value has to be more "
                                 "than 0.01 V as a value of zero would lead to an infinite loop.");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::PolarizationCurve<dim>::initialize(ParameterHandler& param)
{
    param.enter_subsection("Simulator");
    {
        param.enter_subsection("Polarization Curve");
        {
            this->parameter_filename = param.get("Polarization curve file output");
            this->p_init = param.get_double("Initial voltage [V]");
            this->p_end = param.get_double("Final voltage [V]");
            this->dp =  param.get_double("Increment [V]");
            this->adaptive = param.get_bool("Adaptive Increment");
            this->min_dp = param.get_double("Min. Increment [V]");
            if (this->min_dp < 0.01)
                this->min_dp = 0.01;                
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    this->n_dpPts = floor( (this->p_init - this->p_end)/this->dp ) + 1;
    
    // Initialize convergence flag
    this->convergence = false;
    
    // Initialize coarse solution w/ empty vector:
    this->coarse_solution = FuelCell::ApplicationCore::FEVector();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// PRIVATE:
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::PolarizationCurve<dim>::set_parameters(ParameterHandler& param,
                                                 const shared_ptr<FuelCell::ApplicationCore::AdaptiveRefinement<dim> >& solver,
                                                 const int iteration,
                                                 const std::string parameter_name,
                                                 const double param_value)
{
    //-- Set design variables to the values given by DAKOTA:
    // Make sure that boundary conditions are modified for cell voltage:
    FcstUtilities::modify_parameter_file("Fuel cell data>>Operating conditions>>Adjust initial solution and boundary conditions", true, param);
    // Modify cell voltage:
    FcstUtilities::modify_parameter_file("Fuel cell data>>Operating conditions>>Voltage cell", param_value, param); 
    // Make sure a solution is not read from file after the first iteration:
    if (iteration != 0)
        FcstUtilities::modify_parameter_file("Initial Solution>>Read in initial solution from file", false, param);
    
    param.enter_subsection("Output Variables");
    param.set("num_output_vars", 1.0);
    param.set("Output_var_0", "current");
    param.leave_subsection();

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::PolarizationCurve<dim>::print_parameters() const
{
    FcstUtilities::log<<"== Polarization curve parameters: =="<<std::endl;
    FcstUtilities::log<<"Initial voltage [V] : "<<this->p_init<<std::endl;
    FcstUtilities::log<<"Final voltage [V] : "<<this->p_end<<std::endl;
    FcstUtilities::log<<"Increment [V] : "<<this->dp<<std::endl;
    FcstUtilities::log<<"Adaptive Increment? "<<this->adaptive<<std::endl;
    FcstUtilities::log<<"Min. Increment [V] : "<<this->min_dp<<std::endl;
    FcstUtilities::log<<"==  =="<<std::endl;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::PolarizationCurve<dim>::print_parameteric_study(const std::vector<std::vector<double> >& curve) const
{
    std::ofstream myfile;
    myfile.open (this->parameter_filename);
    myfile<< "# ===================================================="<<std::endl;
    myfile<< "# OpenFCST: Fuel cell simulation toolbox "<<std::endl;
    myfile<< "# ===================================================="<<std::endl;
    myfile<< "# Polarization curve data :"<<this->parameter_filename<<std::endl;
    myfile<< "# "<<std::endl;
    myfile<< "Cell voltage [V]"<<"\t"<<"Cathode current [A/cm2]"<<std::endl;
    for (unsigned int i = 0; i < curve.size(); i++)
    {
        myfile<<curve[i][0]<<"\t"<<curve[i][1]<<std::endl;
    }
    myfile.close();
}

//---------------------------------------------------------------------------
// Explicit instantations
template class FuelCell::PolarizationCurve<deal_II_dimension>;