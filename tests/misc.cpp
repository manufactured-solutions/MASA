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
// misc.cpp : program that tests masa helper functions
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace MASA;
using namespace std;

int main()
{
  int i;
  double x=0;
  double y=0;
  double z=0;

  string str;
  masa_init("euler-test","euler_1d");
  masa_get_dimension(&i);
  masa_get_name(&str);
  
  if(str.compare("euler_1d") != 0)
    {
      cout << "masa_get_name FAILED";
      return 1;
    }

  if(i!=1)
    {
      cout << "masa_get_dimension FAILED";
      return 1;
    }

  // test a few other functions
  masa_version_stdout();
  masa_get_numeric_version();

  return 0; // steady as she goes

}
