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
// register.cpp : program that tests masa register var
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <cstdio>
#include <iostream>
#include <string>

using namespace MASA;
using namespace std;

typedef double Scalar;

int main()
{
  int err = 0;
  string str;

  // testing masa_init
  masa_init<Scalar>("masa-test","euler_1d");
  err = masa_init_param<Scalar>();
  if(err!=0) // test function returns correct error value for properly built function
    {
      cout << "masa_init success condition not triggered properly!\n";
      return 1;
    }

  err = masa_init<Scalar>("masa-test","masa_test_function");
  if(err!=0) 
    {
      cout << "masa_init success condition not triggered properly!\n";
      return 1;
    }

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa.
  freopen("/dev/null","w",stdout);

  err = masa_init_param<Scalar>();
  if(err!=2) // function designed to fail (twice) for masa_test mms
    {
      cout << "masa_init fail condition not triggered properly!\n";
      return 1;
    }

#ifdef MASA_EXCEPTIONS
  // also will fail -- because user is mucking with vararray
  try
    {
      masa_sanity_check<Scalar>();
    }
  catch(int err2) // return one on fatal error
    {
      if(err2 != 1)
	{
	  cout << "regression test failed: masa sanity_check error condition not triggered!\n";
	  return 1; // fail
	}      
    } 
#endif

  return 0; // steady as she goes

}
