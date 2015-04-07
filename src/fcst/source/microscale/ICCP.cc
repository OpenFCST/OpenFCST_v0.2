//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2011-14 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Class: dummy_micro_scale.cc
//    - Description: A solid spherical catalyst particle surrounded by thin ionomer film
//    - Developers: Philip Wardlaw
//    - $Id: dummy_micro_scale.cc 2605 2014-08-15 03:36:44Z secanell $
//
//---------------------------------------------------------------------------

#include <ICCP.h>
namespace NAME = FuelCellShop::MicroScale;

const std::string NAME::ICCP::concrete_name("ICCP");
NAME::ICCP const* NAME::ICCP::PROTOTYPE = new NAME::ICCP(concrete_name);

//---------------------------------------------------------------------------
NAME::ICCP::ICCP(std::string name): F(Constants::F()),pi(Constants::Pi())
{
    this->get_mapFactory()->insert(
            std::pair<std::string, NAME::ICCP*>(
                    concrete_name, this));
}

//---------------------------------------------------------------------------
void
NAME::ICCP::set_structure(){

    //Convert to cm
    core_radius *=1e-7;
    film_thickness *=1e-7;

    std::map<FuelCellShop::Layer::MultiScaleCL<deal_II_dimension>::Properties, double> props = this->layer->get_properties();
    kinetics = this->layer->get_resource<FuelCellShop::Kinetics::BaseKinetics>();

    //Get active area and scale to surface
    ActiveArea = props[CLPropNames::active_area_scaled]*((4.0/3.0)*pi*std::pow(core_radius,3.0));

    P = this->layer->get_properties()[CLPropNames::pressure];

    boost::shared_ptr<FuelCellShop::Material::PolymerElectrolyteBase> electrolyte;
    electrolyte = this->layer->get_resource<FuelCellShop::Material::PolymerElectrolyteBase>();
    HO2N = electrolyte->get_H_O2();
}


//---------------------------------------------------------------------------
void
NAME::ICCP::declare_parameters (ParameterHandler &param) const{

    param.enter_subsection(concrete_name);{
        param.declare_entry("Radius [nm]", "100.0", Patterns::Double(0.0,std::numeric_limits<double>::max()));
        param.declare_entry("Film Thickness [nm]", "5.0", Patterns::Double(0.0,std::numeric_limits<double>::max()));
        param.declare_entry("Non Equilibrium BC Rate constant", "0.13", Patterns::Double());
        param.declare_entry("Use non equilibrium BC", "false", Patterns::Bool());
    }
    param.leave_subsection();

}


//---------------------------------------------------------------------------
void
NAME::ICCP::initialize (ParameterHandler &param){

    param.enter_subsection(concrete_name);{
        core_radius  = param.get_double("Radius [nm]");
        film_thickness  = param.get_double("Film Thickness [nm]");
        non_eq_BC_coeff = 100*param.get_double("Non Equilibrium BC Rate constant"); //convert to cm/s
        non_eq_BC  = param.get_bool("Use non equilibrium BC");
    }
    param.leave_subsection();
}


//---------------------------------------------------------------------------
void
NAME::ICCP::set_solution(const std::map<VariableNames,SolutionVariable>& sols,const VariableNames& react, const int& index){


    reactants.clear();

    AssertThrow((react == oxygen_molar_fraction) or (react == hydrogen_molar_fraction) ,
            ExcMessage("Agglomerate cannot solve for this type of reactant."));

    AssertThrow(sols.find(protonic_electrical_potential) != sols.end(),
            ExcMessage("Solution is missing protonic potential!"));

    AssertThrow((sols.find(electronic_electrical_potential) != sols.end()),
            ExcMessage("Solution is missing solid potential!"));



    P = this->layer->get_properties()[CLPropNames::pressure];

    reactants.push_back( SolutionVariable(sols.at(oxygen_molar_fraction)[index] * (P/HO2N), 1, oxygen_concentration));
    proton_pot = SolutionVariable(sols.at(protonic_electrical_potential)[index], 1, protonic_electrical_potential);
    electron_pot = SolutionVariable(sols.at(electronic_electrical_potential)[index], 1, electronic_electrical_potential);

    double temp = sols.at(temperature_of_REV)[index];
    electrolyte = this->layer->get_resource<FuelCellShop::Material::PolymerElectrolyteBase>();
    electrolyte->set_T(temp);
    kinetics->set_temperature(SolutionVariable(temp,1,temperature_of_REV));
    electrolyte->oxygen_diffusivity(DO2N);

}


//---------------------------------------------------------------------------
FuelCellShop::SolutionMap
NAME::ICCP::compute_current(){
    //Same as the Macro-Hom current term, just with a Macro:Micro effectiveness factor



    const double c_outer = reactants[0][0];
    
    if(c_outer <= 0.0){
        SolutionMap sols;
        sols.push_back(SolutionVariable(0, 1, current_density));
        sols.push_back(SolutionVariable(0,1, CL_effectiveness));

        return sols; //To stabalize FEM solution
    }

    double c_inner = c_outer;

    //Get surface reaction rate
    std::vector<double> J(1,0.0);
    kinetics->set_reactant_concentrations(reactants);
    kinetics->set_electrolyte_potential(proton_pot);
    kinetics->set_solid_potential(electron_pot);
    kinetics->current_density(J);
    double J_ideal = J[0];

    if(film_thickness > 0.0){
        //Do Newton loop
        double f_x = residual(c_inner, c_outer);
        double f_x_; //First derivative
        unsigned int loops = 0;

        while(std::abs(f_x) > 1e-15){

            f_x_ = (residual(c_inner*1.000001, c_outer) - residual(c_inner*0.999999, c_outer))/(0.000002*c_inner);

            double step_factor = 1.0;
            while(c_inner - step_factor*(f_x/f_x_) < 0.0){
                //if the step will take c_inner below 0 then reduce it
                step_factor /=2.0;

                if(step_factor < 1e-5) //Something has gone really wrong...
                    throw std::runtime_error("ICCP inner CO2 value failed to converge (Step reduction loop).");
            }


            c_inner -= step_factor*(f_x/f_x_);
            f_x = residual(c_inner, c_outer);

            if(loops++ > 1000)
                throw std::runtime_error("ICCP inner CO2 value failed to converge.");


        }

        //CO_2 inner solved, get current density
        std::vector<SolutionVariable> temp_react;
        temp_react.push_back(SolutionVariable(c_inner, 1, oxygen_concentration));
        kinetics->set_reactant_concentrations(temp_react);
        kinetics->current_density(J);
    }


    SolutionMap sols;
    sols.push_back(SolutionVariable(J[0]*ActiveArea/((4.0/3.0)*pi*std::pow(core_radius,3.0)), 1, current_density));
    sols.push_back(SolutionVariable(J[0]/J_ideal,1, CL_effectiveness));

    if(this->kinetics->has_coverage(OH_coverage)){
        std::vector<double> OH_c;
        this->kinetics->OH_coverage(OH_c);
        sols.push_back(SolutionVariable(OH_c,OH_coverage));
    }

    if(this->kinetics->has_coverage(O_coverage)){
        std::vector<double> O_c;
        this->kinetics->O_coverage(O_c);
        sols.push_back(SolutionVariable(O_c,O_coverage));
    }


    return sols;

}


//---------------------------------------------------------------------------
double
NAME::ICCP::residual(const double & c_inner, const double & c_outer){
    double answer = 0.0;

    std::vector<double> J(1,0.0);
    std::vector<SolutionVariable> temp_react;
    temp_react.push_back(SolutionVariable(c_inner, 1, oxygen_concentration));
    kinetics->set_reactant_concentrations(temp_react);
    kinetics->set_electrolyte_potential(proton_pot);
    kinetics->set_solid_potential(electron_pot);
    kinetics->current_density(J);


    if (non_eq_BC)
    {
        double a = ((J[0]*ActiveArea)/(16*F*pi*std::pow(core_radius + film_thickness,2.0)*non_eq_BC_coeff));
        double b= (film_thickness/(core_radius*(core_radius+film_thickness)))*((J[0]*ActiveArea)/(16*F*pi*DO2N));
        answer = c_inner + a + b - c_outer;
    }
    else
    {
        answer = c_inner - c_outer + (film_thickness/(core_radius*(core_radius+film_thickness)))*((J[0]*ActiveArea)/(16*F*pi*DO2N));
    }

    return answer;
}


//---------------------------------------------------------------------------
void
NAME::ICCP::print_properties(){
    FcstUtilities::log << "=========== CL MICROSTRUCTURE =========" << std::endl;
    FcstUtilities::log <<  "Agglomerate Type: " << get_name() << std::endl;
    FcstUtilities::log <<  "Core radius    [nm]: " << core_radius*1e7 << std::endl;
    FcstUtilities::log <<  "Film thickness [nm]: " << film_thickness*1e7 << std::endl;
    FcstUtilities::log <<  "Boundary conditions: " << (non_eq_BC ? "Non Equilibrium": "Equilibrium") << std::endl;
    FcstUtilities::log << "=======================================" << std::endl;
}




//---------------------------------------------------------------------------
void
NAME::ICCP::make_thread_safe(ParameterHandler &param, unsigned int thread_index){

    catalyst = FuelCellShop::Material::CatalystBase::create_Catalyst(param, param.get("Catalyst type"));
    electrolyte = FuelCellShop::Material::PolymerElectrolyteBase::create_PolymerElectrolyte(param, param.get("Electrolyte type"));
    kinetics = FuelCellShop::Kinetics::BaseKinetics::create_Kinetics(param, param.get("Kinetics type"));
    kinetics->set_catalyst(catalyst.get());
    kinetics->set_electrolyte(electrolyte.get());

    ReactionNames rxn_name;

    if((typeid(*kinetics.get()) == typeid(FuelCellShop::Kinetics::TafelKinetics)) or (typeid(*kinetics.get()) == typeid(FuelCellShop::Kinetics::DoubleTrapKinetics))){
        rxn_name = ORR;
    }
    else if(typeid(*kinetics.get()) == typeid(FuelCellShop::Kinetics::DualPathKinetics)){
        rxn_name = HOR;

    }


    kinetics->set_reaction_kinetics(rxn_name);
    catalyst->set_reaction_kinetics(rxn_name);
}
