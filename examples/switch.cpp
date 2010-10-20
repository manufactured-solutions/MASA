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
// $Author$
// $Id$
//
// switch.cpp:
// C++ example demonstrating how to switch between two different solutions
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

// header contains all subroutine information
#include <masa.h>
#include <iostream>

// MASA was designed with a custom namespace
using namespace MASA;
using namespace std;

typedef double Scalar;

int main()
{
  Scalar solution;

  // print all solutions
  masa_printid<Scalar>();
  cout << endl << endl;

  // initialize first solution
  masa_init<Scalar>("alice","heateq_1d_steady_const");
  masa_init_param<Scalar>();

  // initialize 2nd solution
  masa_init<Scalar>("bob"  ,"euler_2d");
  masa_init_param<Scalar>();
  
  // list all initialized solutions
  masa_list_mms<Scalar>();

  // lets manipulate the alice set
  masa_select_mms<Scalar>("alice");
  masa_display_param<Scalar>();
  solution = masa_eval_t_source<Scalar>(1.2);
  cout << solution << endl;

  // now switch to and edit bob
  masa_select_mms<Scalar>("bob");
  masa_display_param<Scalar>();
  solution = masa_eval_u_source<Scalar>(1,1);
  cout << solution << endl;

}

