// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: simulation_selector.cc
// - Description: This class selects an openFCST application which will run
// - Developers: P. Dobson,
//               M. Secanell,
//               A. Koupaei,
//               V. Zingan,
//               M. Bhaiya
// - Id: $Id$
//
// ----------------------------------------------------------------------------

#include "simulation_selector.h"

// ----------------------------------------------------------------------------
template<int dim>
SimulationSelector<dim>::SimulationSelector()
{ }

// ----------------------------------------------------------------------------
template<int dim>
SimulationSelector<dim>::~SimulationSelector()
{ }

// ----------------------------------------------------------------------------
template<int dim>
void
SimulationSelector<dim>::declare_parameters(ParameterHandler& param) const
{
     param.enter_subsection("Simulator");
     {
          param.declare_entry("simulator name",
                              "cathode",
                               Patterns::Selection( get_simulator_names() ),
                              "Name of the application that you would like to run ");

          param.declare_entry("simulator specification",
                              "Poiseuille",
                               Patterns::Selection( get_simulator_specifications() ),
                              "Select a sub-application ");

          param.declare_entry("solver name",
                              "NewtonLineSearch",
                               Patterns::Selection( get_solver_names() ),
                              "Select type of linear or non-linear solver ");

          param.declare_entry("solver method",
                              "AdaptiveRefinement",
                               Patterns::Selection( get_solver_methods() ),
                              " ");
     }
     param.leave_subsection();
}

// ----------------------------------------------------------------------------
template<int dim>
void
SimulationSelector<dim>::initialize(ParameterHandler& param)
{
     param.enter_subsection("Simulator");
     {
          name_application  = param.get("simulator name");
          app_specification = param.get("simulator specification");
          name_solver       = param.get("solver name");
          name_solve_method = param.get("solver method");
     }
     param.leave_subsection();
}

// ----------------------------------------------------------------------------
template<int dim>
shared_ptr< FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim> >
SimulationSelector<dim>::select_application()
{
    boost::shared_ptr< FuelCell::ApplicationCore::ApplicationData > data(new FuelCell::ApplicationCore::ApplicationData);
    data->enter_flag(app_specification,
                     true);

    //-- Cathode without convection, pseudo-homogeneous
    if (name_application.compare("cathode") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY SOLVING A CATHODE MODEL" << std::endl;
        return shared_ptr<FuelCell::Application::AppCathode<dim> > (new FuelCell::Application::AppCathode<dim>);
    }
    //-- Complete MEA without convection
    else if (name_application.compare("MEA") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY SOLVING AN MEA MODEL" << std::endl;
        return shared_ptr<FuelCell::Application::AppPemfc<dim> > (new FuelCell::Application::AppPemfc<dim>);

    }
    //Run application defined in app_thermal_test.h
    else if (name_application.compare("thermalTest") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY SOLVING A Thermal Transport Equation TEST CASE" << std::endl;
        return shared_ptr<FuelCell::Application::AppThermalTest<dim>> (new FuelCell::Application::AppThermalTest<dim>);
    }
    //Run application defined in app_pemfc_nonisothermal.h
    else if (name_application.compare("meaNIT") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY SOLVING A NON-ISOTHERMAL MEA MODEL" << std::endl;
        return shared_ptr<FuelCell::Application::AppPemfcNIThermal<dim>> (new FuelCell::Application::AppPemfcNIThermal<dim>);
    }
    //run mesh test application
    else if (name_application.compare("test_mesh") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY RUNNING THE TEST MESH APPLICATION" << std::endl;
        return shared_ptr<FuelCell::Application::AppReadMesh<dim> > (new FuelCell::Application::AppReadMesh<dim>);
    }
    else if (name_application.compare("diffusion") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY RUNNING THE DIFFUSION APPLICATION" << std::endl;
        return shared_ptr<FuelCell::Application::AppDiffusion<dim> > (new FuelCell::Application::AppDiffusion<dim>);
    }
    else if (name_application.compare("ohmic") == 0)
    {
        FcstUtilities::log << "YOU ARE CURRENTLY RUNNING THE OHMIC ELECTRON TRANSPORT APPLICATION" << std::endl;
        return shared_ptr<FuelCell::Application::AppOhmic<dim> > (new FuelCell::Application::AppOhmic<dim>);
    }
    else
    {
        FcstUtilities::log << "Application not found. See " << __FILE__ << " and " << __FUNCTION__ << std::endl;
        Assert(false, ExcNotImplemented());
        FcstUtilities::log << "Application not implemented. See " << __FILE__ << std::endl;
    }
}

// ----------------------------------------------------------------------------
template<int dim>
shared_ptr< FuelCell::ApplicationCore::ApplicationWrapper >
SimulationSelector<dim>::select_solver(FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>* app_lin)
{
     if  (name_solver.compare("Linear") == 0)
     {
        FcstUtilities::log << "YOUR PROBLEM IS LINEAR" << std::endl;
        return shared_ptr<FuelCell::ApplicationCore::ApplicationWrapper> (new FuelCell::ApplicationCore::ApplicationWrapper(*app_lin));
     }
     else if( name_solver.compare("NewtonBasic") == 0 )
     {
        FcstUtilities::log << "YOU ARE USING NewtonBasic NEWTON SOLVER" << std::endl;
        return shared_ptr<FuelCell::ApplicationCore::NewtonBasic> ( new FuelCell::ApplicationCore::NewtonBasic(*app_lin) );
     }
     else if (name_solver.compare("Newton3pp") == 0)
     {
        FcstUtilities::log << "YOU ARE USING Newton3pp NEWTON SOLVER" << std::endl;
        return shared_ptr<FuelCell::ApplicationCore::Newton3pp> (new FuelCell::ApplicationCore::Newton3pp(*app_lin));
     }
     else if (name_solver.compare("NewtonLineSearch") == 0)
     {
        FcstUtilities::log << "YOU ARE USING NewtonLineSearch NEWTON SOLVER" << std::endl;
        return shared_ptr<FuelCell::ApplicationCore::NewtonLineSearch> (new FuelCell::ApplicationCore::NewtonLineSearch(*app_lin));
     }
     else
     {
        Assert(false, ExcNotImplemented());
        FcstUtilities::log << "This newton solver is not implemented yet! See " << __FILE__ << std::endl;
     }
}

// ----------------------------------------------------------------------------
template<int dim>
shared_ptr< FuelCell::ApplicationCore::AdaptiveRefinement<dim> >
SimulationSelector<dim>::select_solver_method(FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim>* app_lin,
                                              FuelCell::ApplicationCore::ApplicationWrapper*                      newton_solver,
                                              const FuelCell::ApplicationCore::FEVector&                          solution)
{
    if (name_solve_method.compare("AdaptiveRefinement") == 0)
    {
        return shared_ptr<FuelCell::ApplicationCore::AdaptiveRefinement<dim> > (new FuelCell::ApplicationCore::AdaptiveRefinement<dim> (*app_lin, *newton_solver, solution));
    }
    else
    {
        Assert(false, ExcNotImplemented());
        FcstUtilities::log << "This solver is not implemented yet! See " << __FILE__ << std::endl;
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
template class SimulationSelector<deal_II_dimension>;

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
