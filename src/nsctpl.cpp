// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// $Author$
// $Id$
//
// nsctpl.cpp: These are the MASA class member functions and constructors for
// the Navier--Stokes compressible transient power law viscosity case.  The
// implementation is svn:external-ed in as nsctpl_fwd.hpp and nsctpl.hpp.
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h>

// External nsctpl::manufactured_solution implementation
namespace MASA {
#include "nsctpl.hpp"
}

template <typename Scalar>
struct registration_helper
{
  MASA::navierstokes_4d_compressible_powerlaw<Scalar> *p;

  registration_helper(MASA::navierstokes_4d_compressible_powerlaw<Scalar> *p) : p(p) {}

  void operator() (const std::string &name, Scalar &value) {
    p->register_var(name.c_str(), &value);
  }
};

template <typename Scalar>
MASA::navierstokes_4d_compressible_powerlaw<Scalar>::navierstokes_4d_compressible_powerlaw()
{
  this->mmsname   = "navierstokes_4d_compressible_powerlaw";
  this->dimension = 4; // x + y + z + t = 4

  // Register parameters using nsctpl::manufactured_solution::foreach_parameter
  registration_helper<Scalar> rh(this);
  this->foreach_parameter(rh);

  this->init_var();
}

template <typename Scalar>
struct set_helper
{
  MASA::navierstokes_4d_compressible_powerlaw<Scalar> *p;

  set_helper(MASA::navierstokes_4d_compressible_powerlaw<Scalar> *p) : p(p) {}

  void operator() (const std::string &name, Scalar &value) {
    p->set_var(name.c_str(), value);
  }
};

template <typename Scalar>
int MASA::navierstokes_4d_compressible_powerlaw<Scalar>::init_var()
{
  // Set parameter values directly on the instance
  MASA::nsctpl::isothermal_channel(*this);

  // Tell MASA about each parameter so it believes they are initialized
  set_helper<Scalar> sh(this);
  this->foreach_parameter(sh);

  return 0;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::navierstokes_4d_compressible_powerlaw);
