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
// - Class: application_data.h
// - Description: This class implements general data of applications
// - Developers: Guido Kanschat,     Texas A&M University
//               Valentin N. Zingan, University of Alberta
//               Marc Secanell,      University of Alberta
// - Id: $Id$
//
// ----------------------------------------------------------------------------

#ifndef _FUEL_CELL_APPLICATION_CORE_APPLICATION_DATA_H_
#define _FUEL_CELL_APPLICATION_CORE_APPLICATION_DATA_H_

#include <lac/vector_memory.h>
#include <lac/block_vector.h>
#include "logging.h"
#include <algorithm>
#include <fstream>

using namespace dealii;

namespace FuelCell
{
namespace ApplicationCore
{

/**
 * The vector class used by applications.
 */
typedef BlockVector<double> FEVector;

/**
 * Here we handle general data of applications. It is the idea of this
 * class that all applications inheriting from ApplicationBase share
 * the <b>same</b> ApplicationData object, thus accessing the same
 * vector pools and using the object as a means for communicating
 * global information.
 *
 * First purpose is providing vector memory objects for Vector and
 * BlockVector objects for applications. These objects, #vector_pool
 * and #block_vector_pool can be accessed directly.
 *
 * Furthermore, it provides a map which allows applications to
 * exchange data in named registers. To this end, the functions
 * enter_flag(std::string,const bool&), enter(std::string,const double&),
 * enter(std::string,const FEVector&), and enter(std::string,const std::vector<double>&)
 * can be used to enter boolean flag, scalar and vector data, respectively.
 * These data are entered as references, so
 * the controlling application can change them at any time. All other
 * applications sharing the same ApplicationData can access them read
 * only through the functions flag(), scalar(), vector(), and std_vector().
 *
 * @author Guido Kanschat
 * @author Valentin N. Zingan
 */

class ApplicationData
{
public:

  /**
   * Destructor, releasing all data.
   */
 ~ApplicationData();

  /**
   * Register a named boolean flag.
   */
  void enter_flag(std::string name,
                  const bool& s);

  /**
   * Register a named scalar.
   */
  void enter(std::string   name,
             const double& s);

  /**
   * Register a named vector.
   */
  void enter(std::string     name,
             const FEVector& v);

  /**
   * Register a named std vector.
   */
  void enter(std::string                name,
             const std::vector<double>& v);

  /**
   * Delete a previously registered
   * boolean flag.
   */
  void erase_flag(std::string name);

  /**
   * Delete a previously registered
   * scalar.
   */
  void erase_scalar(std::string name);

  /**
   * Delete a previously registered
   * vector.
   */
  void erase_vector(std::string name);

  /**
   * Delete a previously registered
   * std vector.
   */
  void erase_std_vector(std::string name);

  /**
   * This function returns @p true
   * if a boolean flag exists (doesn't matter if the flag itself is @p true or @p false).
   * Otherwise returns @p false.
   */
  bool flag_exists(const std::string& name) const;

  /**
   * Get read-only access to a
   * registered boolean flag.
   */
  bool flag(std::string name) const;

  /**
   * Get read-only access to a
   * registered scalar.
   * It only offers read access and
   * returns a null pointer, if
   * the name has not been
   * registered.
   */
  const double* scalar(std::string name) const;

  /**
   * Get read-only access to a
   * registered vector.
   * It only offers read access and
   * returns a null pointer, if
   * the name has not been
   * registered.
   */
  const FEVector* vector(std::string name) const;

  /**
   * Get read-only access to a
   * registered std vector.
   * It only offers read access and
   * returns a null pointer, if
   * the name has not been
   * registered.
   */
  const std::vector<double>* std_vector(std::string name) const;

  /**
   * List all stored objects to FcstUtilities::log.
   */
  void log() const;

  /**
   * VectorMemory object for
   * simple vectors. All
   * applications should allocate
   * their vectors here, so they
   * can be reused and operating
   * system memory management can
   * be avoided.
   */
  GrowingVectorMemory< Vector<double> > vector_pool;

  /**
   * VectorMemory object for
   * block vectors. All
   * applications should allocate
   * their vectors here, so they
   * can be reused and operating
   * system memory management can
   * be avoided.
   */
  GrowingVectorMemory<FEVector> block_vector_pool;

  /**
   * The typedef for the map of
   * boolean flags.
   */
  typedef std::map< std::string, bool > flag_map;

  /**
   * The typedef for the map of
   * scalars.
   */
  typedef std::map< std::string, const double* > scalar_map;

  /**
   * The typedef for the map of
   * vectors.
   */
  typedef std::map< std::string, const FEVector* > vector_map;

  /**
   * The typedef for the map of
   * std vectors.
   */
  typedef std::map< std::string, const std::vector<double>* > std_vector_map;

  /**
   * Exception thrown when a
   * named scalar or vector was
   * searched but not found.
   */
  DeclException2(ExcNotFound,
                 char*,
                 std::string,
                 << "A " << arg1 << " with name " << arg2 << " was not stored in this data object");

private:

  /**
   * A map linking names of data
   * to actual boolean flags.
   */
  flag_map named_flags;

  /**
   * A map linking names of data
   * to actual scalars.
   */
  scalar_map named_scalars;

  /**
   * A map linking names of data
   * to actual vectors.
   */
  vector_map named_vectors;

  /**
   * A map linking names of data
   * to actual std vectors.
   */
  std_vector_map named_std_vectors;
};

} // ApplicationCore

} // FuelCell

#endif