//---------------------------------------------------------------------------
// C++ Interface: fcst_units.cc
//
// Description: A class consisting of static members that can be used by fcst 
// 			programmers when manipulating units in order to preserve
//			standards and prevent confusion.
//
// Author: Philip Wardlaw 2013
//
// Copyright: See COPYING file that comes with this distribution
//
//---------------------------------------------------------------------------


#include "fcst_units.h"


double Units::PER_UNIT =  1;
double Units::PER_C_UNIT =  1E-2;
double Units::PER_MILLI_UNIT =  1E-3;
double Units::PER_MICRO_UNIT =  1E-6;
double Units::PER_N_UNIT =  1E-9;
double Units::PER_P_UNIT =  1E-12;
	
double Units::PER_UNIT2 =  1;
double Units::PER_C_UNIT2 =  1E-4;
double Units::PER_MILLI_UNIT2 =  1E-6;
double Units::PER_MICRO_UNIT2 =  1E-12;
double Units::PER_N_UNIT2 =  1E-18;
double Units::PER_P_UNIT2 =  1E-24;
	
double Units::PER_UNIT3 =  1;
double Units::PER_C_UNIT3 =  1E-6;
double Units::PER_MILLI_UNIT3 =  1E-9;
double Units::PER_MICRO_UNIT3 =  1E-18;
double Units::PER_N_UNIT3 =  1E-27;
double Units::PER_P_UNIT3 =  1E-36;

double Units::UNIT =  -1;
double Units::C_UNIT = - 1E-2;
double Units::MILLI_UNIT =  -1E-3;
double Units::MICRO_UNIT = - 1E-6;
double Units::N_UNIT =-  1E-9;
double Units::P_UNIT =  -1E-12;
	
double Units::UNIT2 = - 1;
double Units::C_UNIT2 =  -1E-4;
double Units::MILLI_UNIT2 =  -1E-6;
double Units::MICRO_UNIT2 = - 1E-12;
double Units::N_UNIT2 =  -1E-18;
double Units::P_UNIT2 =  -1E-24;
	
double Units::UNIT3 = - 1;
double Units::C_UNIT3 =  -1E-6;
double Units::MILLI_UNIT3 = - 1E-9;
double Units::MICRO_UNIT3 = - 1E-18;
double Units::N_UNIT3 = - 1E-27;
double Units::P_UNIT3 = - 1E-36;