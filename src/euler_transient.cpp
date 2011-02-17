// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010 The PECOS Development Team
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
// $Author: nick $
// $Id: euler.cpp 17621 2011-02-14 16:53:09Z nick $
//
// euler.cpp: These are the MASA class member functions and constructors
//          For the Euler Equations
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         EULER EQUATION 1D
 *
 *
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::euler_transient_1d<Scalar>::euler_transient_1d()
{
  this->mmsname = "euler_1d";
  this->dimension=1;

}

template <typename Scalar>
int MASA::euler_transient_1d<Scalar>::init_var()
{


  return 0;
}


template <typename Scalar>
Scalar MASA::euler_transient_1d<Scalar>::eval_q_rho_u(Scalar x)
{

  return 0;

}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::euler_transient_1d);
//MASA_INSTANTIATE_ALL(MASA::euler_transient_2d);
//MASA_INSTANTIATE_ALL(MASA::euler_transient_3d);

