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
// $Id: misc.cpp 13822 2010-09-21 18:14:25Z nick $
//
// fail.cpp : program that tests masa failures, of various functions
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <iostream>
#include <stdio.h>

using namespace MASA;
using namespace std;

typedef double Scalar;

int main()
{
  int test;

  masa_init<Scalar>("test","masa_test_function");

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.
  freopen("/dev/null","w",stdout);

  masa_set_param<Scalar>("demo_var_2",0);
  masa_set_param<Scalar>("A_1",0); // does not exist
  // we are intentionally not setting demo_var_3

  // test error on sanity check when user has not initialized variables
  int err = masa_sanity_check<Scalar>();
  if(err!=1)
    {
      cout << "masa_sanity_check FAILED\n";
      return 1;
    }

  // now try to 'get_var' that does not exist
  test = masa_get_param<Scalar>("A_1");
  if(test != -20)
    {
      cout << "REGRESSION TEST FAILED: masa_get_param error condition not triggered\n";
    }

}
