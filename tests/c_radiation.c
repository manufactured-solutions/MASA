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
// $Id: heat-eq.cpp 18184 2011-03-01 20:09:57Z nick $
//
// c_radiation.cpp: regression testing radiation solution class
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include<tests.h>

int main()
{
  int err=0;
  //int n = 6;
  //double vec1[n];

  // initialize the problem
  err += masa_init("radiation","radiation_integrated_intensity");

  // initialize the default parameters
  err += masa_init_param();

  // intialize the various parameters required for Euler 2D
  // call the sanity check routine 
  // (tests that all variables have been initialized)
  err += masa_sanity_check();

  double x=12;
  double source = masa_eval_1d_source_u(x);
  double exact  = masa_eval_1d_exact_u(x);

  nancheck(source);
  nancheck(exact);

  return err;

}// end program
