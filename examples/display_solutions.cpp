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
// $Id: verify.cpp 12695 2010-08-26 03:47:26Z nick $
//
//
// display_solutions.cpp: when run, will display all available 
//                        masa solutions in STDOUT
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
#include <iostream>
#include <math.h>
#include <masa.h>
#include <stdio.h>

using namespace MASA;
using namespace std;

int main()
{
  // reroute stdout: comment to display to screen
  freopen("/dev/null","w",stdout);

  // show available solutions
  return masa_printid<double>();

}
