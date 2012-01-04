// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012 The PECOS Development Team
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
// $Id: rans_sa.cpp 18184 2011-03-01 20:09:57Z nick $
//
// fans_sa_transient_d_finite: this is an example of the API used for 
//                             calling the wall-bounded spelart alamaras 
//                             Favre Average Navier Stokes (FANS) model
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>

using namespace MASA;
typedef double Scalar;

int main()
{
  int err;

  // initialize the problem 
  err = masa_init<Scalar>("sa example","fans_sa_steady_wall_bounded");

  // call the sanity check routine 
  // (tests that all variables have been initialized)
  err = masa_sanity_check<Scalar>();

  // stay tuned: obviously this is not fully implimented yet


  return err;

}// end program
