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
// $Id: misc.cpp 13998 2010-09-23 18:31:43Z nick $
//
// pass_func.cpp: passes a function into the masa realm
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <iostream>
#include <string>

using namespace MASA;
using namespace std;

//typedef long double Scalar;
typedef double Scalar;

Scalar tester()
{
  
  double a = 4.4;
  return a;

}


int main()
{
  Scalar u_0;
  Scalar q = 0;

  // start problem
  masa_init<Scalar>("masa-test","euler_chem_1d");

  u_0 = 1.234567890123456789;
  pass_func<Scalar>(&tester,u_0);

  q = tester();

  cout << "q is: " << q << endl;

  // steady as she goes
  return 0;

}
