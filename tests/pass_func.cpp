// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011 The PECOS Development Team
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

#include<tests.h>

using namespace MASA;
using namespace std;

//typedef long double Scalar;
typedef double Scalar;


int main()
{

  const Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();

  // start problem
#ifndef portland_compiler
  Scalar u_0;
  masa_init<Scalar>("masa-test","euler_chem_1d");

  u_0 = 1.234567890123456789;
  Scalar out = pass_func<Scalar>(&tester,u_0);
  Scalar q = tester(u_0);
  
  if((out - q) > threshold) 
    {
      cout << "\nMASA ERROR::function passing not functioning properly in MASA!\n";
      return 1;
    }

#endif

  // steady as she goes
  return 0;

}
