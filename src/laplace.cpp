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
// $Author: roystgnr $
// $Id: heat.cpp 20657 2011-05-03 19:05:52Z roystgnr $
//
// laplace.cpp: These are the MASA class member functions and constructors
//          For Laplace's equation
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         HEAT EQUATION
 *
 *   steady -- constant coeff
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::laplace_2d<Scalar>::laplace_2d()
{

  this->mmsname = "laplace_2d";
  this->dimension=2;
  
  this->register_var("Lx",&Lx);   
  this->register_var("Ly",&Ly);
  
  this->init_var();
  
}

template <typename Scalar>
int MASA::laplace_2d<Scalar>::init_var()
{

  int err = 0;

  err += this->set_var("Lx",1.4);
  err += this->set_var("Ly",.82);

  return err;
  
}

template <typename Scalar>
Scalar MASA::laplace_2d<Scalar>::eval_q_f(Scalar x,Scalar y)
{
  using std::pow;

  Scalar f;
  f  = 2*pow(Lx-x,2) - 8*(Lx-x)*(Lx+x) + 2*pow(Lx+x,2);
  f += 2*pow(Ly-y,2) - 8*(Ly-y)*(Ly+y) + 2*pow(Ly+y,2);
  return f;
}

template <typename Scalar>
Scalar MASA::laplace_2d<Scalar>::eval_exact_phi(Scalar x,Scalar y)
{
  using std::pow;

  Scalar phi;
  phi = pow(Ly-y,2)*pow(Ly+y,2) + pow(Lx-x,2)*pow(Lx+x,2);
  return phi;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::laplace_2d);
