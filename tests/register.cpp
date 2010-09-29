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
// register.cpp : program that tests masa register var
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

using namespace MASA;
using namespace std;

int main()
{
  int i,err;
  double x=0;
  double y=0;
  double z=0;
  string str;

  // testing masa_init
  masa_init("masa-test","euler_1d");
  err = masa_init_param();
  if(err!=0) // test function returns correct error value for properly built function
    {
      cout << "masa_init success condition not triggered properly!\n";
      return 1;
    }

  err = masa_init("masa-test","masa_test");
  if(err!=0)
    {
      cout << "masa_init FAILED\n";
      return 1;
    }

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa.
  //freopen("/dev/null","w",stdout);

  err = masa_init_param();
  if(err!=1) // function designed to fail for masa_test mms
    {
      cout << "masa_init fail condition not triggered properly!\n";
      return 1;
    }
    
  return 0; // steady as she goes

}
