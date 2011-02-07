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
// catch_exception.cpp : program that tests fatal masa failures
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <iostream>
#include <stdio.h>

using namespace MASA;
using namespace std;

int main()
{

  // only test with masa exception handling turned on
#ifdef MASA_EXCEPTIONS
  
  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.
  freopen("/dev/null","w",stdout);

  // verify_pointer_sanity
  try
    {
      masa_eval_source_t<double>(1.0);
    }
  catch(int err) // return one on fatal error
    {
      if(err != 1)
	{
	  cout << "regression test failed: masa pointer sanity fatal error condition not triggered!\n";
	  return 1; // fail
	}      
    }
 
  // masa_init
  try
    {
      masa_init<double>("it is pitch black","grue"); // grue is not a valid masa mms, nor should it be
    }
  catch(int err) // return one on fatal error
    {
      if(err != 1)
	{
	  cout << "regression test failed: masa_init fatal error condition not triggered!\n";
	  return 1; // fail
	}      
    }
 
  // masa_select_mms
  try
    {
      masa_select_mms<double>("cthulhu"); // cthulhu is not a valid masa solution, nor has it been intialized
    }
  catch(int err)
    {
      if(err != 1)
	{
	  cout << "regression test failed: masa_select_mms fatal error condition not triggered!\n";
	  return 1; // fail
	}      
    }

#endif

  return 0; // success

}
