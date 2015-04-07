//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2009-13 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: parametric_study.cc
//    - Description: Child of ApplicationWrapper used to implement IV curve analysis
//    - Developers: M. Secanell
//    - $Id$
//
//---------------------------------------------------------------------------

#include "utils/parametric_study.h"

//---------------------------------------------------------------------------
template <int dim>
FuelCell::ParametricStudy<dim>::ParametricStudy()
{}

//---------------------------------------------------------------------------
template <int dim>
FuelCell::ParametricStudy<dim>::~ParametricStudy()
{}

//---------------------------------------------------------------------------

template <int dim>
void
FuelCell::ParametricStudy<dim>::declare_parameters(ParameterHandler& param) const
{
    // Initialize to parameters needed for evaluate (note: this is provisional)
    param.enter_subsection("Simulator");
    {
        param.enter_subsection("Parametric study");
        {
            param.declare_entry("Parameter file output",
                                "parameteric_study.dat",
                                Patterns::Anything(),
                                "File where the parametric study results should be stored");
            param.declare_entry("Parameter name",
                                "no_name",
                                Patterns::Anything(),
                                "Enter the name of the parameter you would like to study. Use one of the following formats: \n"
                                "For normal parameter: Subsection_1>>Subsection_2>>Value \n"
                                "For boundary value or graded: Subsection_1>>Subsection_2>>Material_id:Value \n"
                                "where Subsection_1 and Subsection_2 would be the sections where the parameter is found in the data file");
            param.declare_entry ("Initial value",
                                 "1",
                                 Patterns::Double(),
                                 "Enter the value you would like to start the parametric study from.");
            param.declare_entry ("Final value",
                                 "0.1",
                                 Patterns::Double(),
                                 "Enter the final value for the parametric study.");
            param.declare_entry ("Increment",
                                 "0.1",
                                 Patterns::Double(),
                                 "Spacing between points in the polarization curve");
            param.declare_entry ("Adaptive Increment",
                                 "true",
                                 Patterns::Bool(),
                                 "Set to true if you would like to reduce the voltage increment adaptively"
                                 "if convergence could not be achieved with the larger value");
            param.declare_entry ("Min. Increment",
                                 "0.025",
                                 Patterns::Double(),
                                 "If Adaptive Increment? is set to true, this value controls the "
                                 "minimum change in cell voltage before the polarization curve gives up"
                                 "and the voltage is updated again. Note that this value has to be more "
                                 "than 0.01 V as a value of zero would lead to an infinite loop.");
            param.declare_entry("Parameter values",
                                """",
                                Patterns::List( Patterns::Double() ),
                                "This list contains the discrete values of a parameter of study.");
        }
        param.leave_subsection();
    }
    param.leave_subsection();
}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ParametricStudy<dim>::initialize(ParameterHandler& param)
{
    param.enter_subsection("Simulator");
    {
        param.enter_subsection("Parametric study");
        {
            parameter_filename = param.get("Parameter file output");
            parameter_name = param.get("Parameter name");
            p_init = param.get_double("Initial value");
            p_end = param.get_double("Final value");
            dp =  param.get_double("Increment");
            adaptive = param.get_bool("Adaptive Increment");
            min_dp = param.get_double("Min. Increment");
            if(!param.get("Parameter values").empty())
            {
              p_values = FcstUtilities::string_to_number<double>( Utilities::split_string_list( param.get("Parameter values") ) );
            }
        }
        param.leave_subsection();
    }
    param.leave_subsection();

    n_dpPts = floor( (p_init-p_end)/dp ) + 1;
}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ParametricStudy<dim>::run(ParameterHandler& param,
                                      std::string simulator_parameter_file_name,
                                      boost::shared_ptr<SimulationSelector<dim> > sim_selector)
{
    std::map<std::string, double> functionals;
    double param_value;
    unsigned int iteration = 0;
    
    FcstUtilities::log<<"============================================================================="<<std::endl;
    print_parameters();
    FcstUtilities::log<<"============================================================================="<<std::endl;
        
    if(p_values.size() == 0) // Marc
    {
        param_value = p_init;
        double param_value_step = dp;

        while( (p_init > p_end) ? (param_value >= p_end) : (param_value <= p_end) )
        {
            // --
            this->print_iteration_info(iteration, param_value);
            
            // -- Solve problem:
            this->run_point(param, simulator_parameter_file_name, sim_selector, iteration, parameter_name, param_value, functionals);

            // -- If convergence achieved, register solution:
            if (convergence)
            {
                register_data(param_value, functionals);
            }
            //-- If convergence not achieved:
            else
            {
                unsigned int step_size = 0;
                double param_value_conv = param_value;

                while (!convergence && adaptive && (param_value_step/pow(2,double(step_size)) > min_dp) )
                {
                    FcstUtilities::log<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
                    FcstUtilities::log<<"Convergence not achieved"<<std::endl;
                    FcstUtilities::log<<"Reducing step size: "<<param_value_step/pow(2,double(step_size))<<std::endl;

                    // Set param_value to previous value
                    param_value_conv += param_value_step/pow(2,double(step_size));
                    // Now, scale
                    step_size += 1;
                    param_value_conv -= param_value_step/pow(2,double(step_size));
                    FcstUtilities::log<<"New param_value: "<<param_value_conv<<std::endl;
                    FcstUtilities::log<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
                    this->run_point(param, simulator_parameter_file_name, sim_selector, iteration, parameter_name, param_value_conv, functionals);

                    if (convergence)
                    {
                        register_data(param_value_conv, functionals);
                        param_value = param_value_conv;
                    }
                }
            }
            //--

            // Update to next param_value:
            param_value -= param_value_step;

            ++iteration;

            FcstUtilities::log<<"============================================================================="<<std::endl;
            FcstUtilities::log<<"============================================================================="<<std::endl;
        }
    }
    else // V
    {
        
        for(unsigned int i = 0; i < p_values.size(); i++)
        {
            param_value = p_values[i];

            //
            this->print_iteration_info(i, param_value);

            // -- Solve problem:
            this->run_point(param, simulator_parameter_file_name, sim_selector, i, parameter_name, param_value, functionals);

            // -- If convergence achieved, register solution:
            if (convergence)
            {
                register_data(param_value, functionals);
            }
            //-- If convergence not achieved:
            else
            {
                FcstUtilities::log << "Convergence can not be achieved at Parameter values [" << i+1 << "] = " << p_values[i] << std::endl;
                AssertThrow(false, ExcInternalError());
            }

            FcstUtilities::log<<"============================================================================="<<std::endl;
            FcstUtilities::log<<"============================================================================="<<std::endl;
        }
    }

    this->print_parameteric_study(curve);
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// PRIVATE:
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ParametricStudy<dim>::set_parameters(ParameterHandler& param,
                                               const shared_ptr<FuelCell::ApplicationCore::AdaptiveRefinement<dim> >& solver,
                                               const int iteration,
                                               const std::string parameter_name,
                                               const double param_value)
{
    //-- Set design variables to the values given by DAKOTA:
    // Modify cell param_value:
    FcstUtilities::modify_parameter_file(parameter_name, param_value, param);
    // Make sure a solution is not read from file after the first iteration:
    if (iteration != 0)
        FcstUtilities::modify_parameter_file("Initial Solution>>Read in initial solution from file", false, param);

    param.enter_subsection("Output Variables");
    param.set("num_output_vars", 1.0);
    param.set("Output_var_0", "current");
    param.leave_subsection();


}

//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ParametricStudy<dim>::run_point(ParameterHandler& param,
                                          const std::string simulator_parameter_file_name,
                                          const boost::shared_ptr<SimulationSelector<dim> > sim_selector,
                                          const int iteration,
                                          const std::string parameter_name,
                                          const double param_value,
                                          std::map<std::string, double>& functionals)
{

    shared_ptr<FuelCell::ApplicationCore::OptimizationBlockMatrixApplication<dim> > app_linear;
    shared_ptr<FuelCell::ApplicationCore::ApplicationWrapper> newton;
    shared_ptr<FuelCell::ApplicationCore::AdaptiveRefinement<dim> > solver;

    // Create linear, nonlinear and adaptive refinement objects:
    app_linear = sim_selector->select_application();
    newton = sim_selector->select_solver(app_linear.get());
    solver = sim_selector->select_solver_method(app_linear.get(), newton.get(), coarse_solution);

    // Declare parameter file:
    solver->declare_parameters(param);
    FcstUtilities::read_parameter_files(param,
                                        simulator_parameter_file_name);

    // Set new parameter in the parameter file:
    this->set_parameters(param, solver, iteration, parameter_name, param_value);

    // Initialize used to create objects, etc.
    solver->initialize(param);

    // -- Allocate space for responses:
    int n_resp = app_linear->get_n_resp();
    std::vector<double> f(n_resp);

    // -- Run simulation
    solver->run_Newton(f);

    // -- Store desired response in functionals:
    std::vector<std::string> name_responses;
    name_responses = app_linear->get_name_responses();
    for (unsigned int i=0; i<n_resp; i++)
    {
        functionals[name_responses[i]] = f[i];
    }

    // -- Was convergence achieved
    convergence = newton->get_data()->flag("Newton convergence");

    // -- If convergence achieved, store coarse mesh solution for future use:
    if (convergence)
        coarse_solution = solver->get_coarse_solution();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ParametricStudy<dim>::print_parameters() const
{
    FcstUtilities::log<<"== Parameteric study parameters: =="<<std::endl;

    if( p_values.size() == 0 )
    {
        FcstUtilities::log<<"Initial value : "<<p_init<<std::endl;
        FcstUtilities::log<<"Final value : "<<p_end<<std::endl;
        FcstUtilities::log<<"Increment : "<<dp<<std::endl;
        FcstUtilities::log<<"Adaptive Increment? "<<adaptive<<std::endl;
        FcstUtilities::log<<"Min. Increment : "<<min_dp<<std::endl;
        FcstUtilities::log<<"==  =="<<std::endl;
    }
    else
    {
        for(unsigned int i = 0; i < p_values.size(); i++)
          FcstUtilities::log << "Parameter values [" << i+1 << "] = " << p_values[i] << std::endl;
    }
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ParametricStudy<dim>::print_parameteric_study(const std::vector<std::vector<double> >& curve) const
{
    std::ofstream myfile;
    myfile.open (parameter_filename);
    myfile<< "# ===================================================="<<std::endl;
    myfile<< "# OpenFCST: Fuel cell simulation toolbox "<<std::endl;
    myfile<< "# ===================================================="<<std::endl;
    myfile<< "# Parametric study data :"<<parameter_filename<<std::endl;
    myfile<< "# Parametric name :"<<parameter_name<<std::endl;
    myfile<< " Parameter value "<<"\t"<<"Cathode current [A/cm2] "<<std::endl;
    for (unsigned int i = 0; i < curve.size(); i++)
    {
        myfile<<curve[i][0]<<"\t"<<curve[i][1]<<std::endl;
    }
    myfile.close();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int dim>
void
FuelCell::ParametricStudy<dim>::register_data(const double param_value,
                                              std::map<std::string, double>& functionals)
{
    std::vector<double> aux;
    aux.push_back(param_value);
    aux.push_back((-1.0)*functionals.find("current")->second);
    curve.push_back(aux);
    FcstUtilities::log<<"Current density: "<<aux[1]<<" at cell param_value "<<aux[0]<<std::endl;
}

//---------------------------------------------------------------------------
// Explicit instantations
template class FuelCell::ParametricStudy<deal_II_dimension>;