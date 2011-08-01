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
// $Author$
// $Id$
//
// switch.cpp: C++ example demonstrating how to switch 
//             between two different solutions
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

// header contains all subroutine information
#include <masa.h>

// MASA was designed with a custom namespace
using namespace MASA;
using namespace std;

typedef double Scalar;

int main()
{
  Scalar solution;

  // reroute stdout: comment to display to screen
  freopen("/dev/null","w",stdout);

  // print all solutions
  masa_printid<Scalar>();

  // initialize first solution
  masa_init<Scalar>("alice","heateq_1d_steady_const");

  // initialize 2nd solution
  masa_init<Scalar>("bob"  ,"euler_2d");

  // list all initialized solutions
  masa_list_mms<Scalar>();

  // lets manipulate the alice set
  masa_select_mms<Scalar>("alice");
  masa_display_param<Scalar>();
  solution = masa_eval_source_t<Scalar>(1.2);

  // now switch to and edit bob
  masa_select_mms<Scalar>("bob");
  masa_display_param<Scalar>();
  solution = masa_eval_source_rho_u<Scalar>(1,1);

}

