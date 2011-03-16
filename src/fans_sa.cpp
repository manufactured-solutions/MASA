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
// $Id: heat.cpp 18162 2011-03-01 05:23:07Z nick $
//
// heat.cpp: These are the MASA class member functions and constructors
//          For Radiation
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *   radiation_integrated_intensity
 *
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::fans_sa_transient_d_finite<Scalar>::fans_sa_transient_d_finite()
{
    this->mmsname = "fans_sa_transient_d_finite";
    this->dimension=1;

    this->register_var("cb1",&cb1);   
    this->init_var();
  
}//done with constructor

template <typename Scalar>
int MASA::fans_sa_transient_d_finite<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("cb1",1.4);
  return err;

}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_q_u(Scalar x)
{
  Scalar Q_T;
  Q_T = (cb1 * x);
  return Q_T;
}

template <typename Scalar>
Scalar MASA::fans_sa_transient_d_finite<Scalar>::eval_exact_u(Scalar x)
{
  Scalar exact_t;
  exact_t = cos(cb1 * x);
  return exact_t;
}

//int x[5] = {2, 3, 5, 7, 11};
//vector<int> vector1(&x[0], &x[5]);

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::fans_sa_transient_d_finite);
