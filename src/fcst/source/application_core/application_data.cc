// ----------------------------------------------------------------------------
//
// FCST: Fuel Cell Simulation Toolbox
//
// Copyright (C) 2006-2009 by Guido Kanschat
// Copyright (C) 2006-2015 by Energy Systems Design Laboratory, University of Alberta
//
// This software is distributed under the MIT license
// For more information, see the README file in /doc/LICENSE
//
// - Class: application_data.cc
// - Description: This class implements general data of applications
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
// - Id: $Id$
//
// ----------------------------------------------------------------------------

#include "application_data.h"

namespace NAME = FuelCell::ApplicationCore;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NAME::ApplicationData::~ApplicationData()
{ }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationData::enter_flag(std::string name,
                                  const bool& s)
{
  named_flags.insert(std::make_pair(name, s));
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationData::enter(std::string   name,
                             const double& s)
{
  named_scalars.insert(std::make_pair(name, &s));
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationData::enter(std::string           name,
                             const NAME::FEVector& v)
{
  named_vectors.insert(std::make_pair(name, &v));
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationData::enter(std::string                name,
                             const std::vector<double>& v)
{
  named_std_vectors.insert(std::make_pair(name, &v));
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationData::erase_flag(std::string name)
{
  flag_map::iterator iter = named_flags.find(name);
  Assert(iter != named_flags.end(), ExcNotInitialized());
  named_flags.erase(iter);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationData::erase_scalar(std::string name)
{
  scalar_map::iterator iter = named_scalars.find(name);
  Assert(iter != named_scalars.end(), ExcNotInitialized());
  named_scalars.erase(iter);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationData::erase_vector(std::string name)
{
  vector_map::iterator iter = named_vectors.find(name);
  Assert(iter != named_vectors.end(), ExcNotInitialized());
  named_vectors.erase(iter);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationData::erase_std_vector(std::string name)
{
  std_vector_map::iterator iter = named_std_vectors.find(name);
  Assert(iter != named_std_vectors.end(), ExcNotInitialized());
  named_std_vectors.erase(iter);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool
NAME::ApplicationData::flag_exists(const std::string& name) const
{
  flag_map::const_iterator iter = named_flags.find(name);

  if(iter != named_flags.end())
    return true;
  else
    return false;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool
NAME::ApplicationData::flag(std::string name) const
{
  flag_map::const_iterator iter = named_flags.find(name);
  AssertThrow(iter != named_flags.end(), ExcNotFound("boolean flag", name));
  return iter->second;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const double*
NAME::ApplicationData::scalar(std::string name) const
{
  scalar_map::const_iterator iter = named_scalars.find(name);
  AssertThrow(iter != named_scalars.end(), ExcNotFound("scalar", name));
  return iter->second;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const NAME::FEVector*
NAME::ApplicationData::vector(std::string name) const
{
  vector_map::const_iterator iter = named_vectors.find(name);
  AssertThrow(iter != named_vectors.end(), ExcNotFound("vector", name));
  return iter->second;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

const std::vector<double>*
NAME::ApplicationData::std_vector(std::string name) const
{
  std_vector_map::const_iterator iter = named_std_vectors.find(name);
  AssertThrow(iter != named_std_vectors.end(), ExcNotFound("std_vector", name));
  return iter->second;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void
NAME::ApplicationData::log() const
{
  FcstUtilities::log.push("boolean flags");
  for(flag_map::const_iterator iter = named_flags.begin(); iter != named_flags.end(); ++iter)
    FcstUtilities::log << iter->first << '\t' << iter->second << std::endl;
  FcstUtilities::log.pop();

  FcstUtilities::log.push("scalar");
  for(scalar_map::const_iterator iter = named_scalars.begin(); iter != named_scalars.end(); ++iter)
    FcstUtilities::log << iter->first << '\t' << iter->second << std::endl;
  FcstUtilities::log.pop();

  FcstUtilities::log.push("vector");
  for(vector_map::const_iterator iter = named_vectors.begin(); iter != named_vectors.end(); ++iter)
    FcstUtilities::log << iter->first << '\t' << iter->second->size() << std::endl;
  FcstUtilities::log.pop();

  FcstUtilities::log.push("std_vector");
  for(std_vector_map::const_iterator iter = named_std_vectors.begin(); iter != named_std_vectors.end(); ++iter)
    FcstUtilities::log << iter->first << '\t' << iter->second->size() << std::endl;
  FcstUtilities::log.pop();
}