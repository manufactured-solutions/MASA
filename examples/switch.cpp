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

// MASA was designed with a custom namespace
using namespace MASA;
using namespace std;

int main()
{
  double solution;

  // print all solutions
  masa_printid();
  cout << endl << endl;

  // initialize first solution
  masa_init("alice","heateq_1d_steady_const");
  masa_init_param();


  // initialize 2nd solution
  masa_init("bob"  ,"euler_2d");
  masa_init_param();
  
  // list all initialized solutions
  masa_list_mms();

  // lets manipulate the alice set
  masa_select_mms("alice");
  masa_display_param();
  masa_eval_t_source(1.2,&solution);
  cout << solution << endl;

  // now switch to and edit bob
  masa_select_mms("bob");
  masa_display_param();
  masa_eval_u_source(1,1,&solution);
  cout << solution << endl;

}

