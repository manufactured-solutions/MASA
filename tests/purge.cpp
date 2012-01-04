// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012 The PECOS Development Team
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
// purge.cpp : program that tests masa purge function
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <tests.h>

#include <cmath>

using namespace MASA;
using namespace std;

typedef double Scalar;
Scalar MASA_DEFAULT = -12345.67;

int main()
{
  Scalar u_0;
  const Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();

  // start problem
  masa_init<Scalar>("masa-test","euler_1d");
  
  // values should be set to something other than default: checking
  u_0 = masa_get_param<Scalar>("u_0");
  if(std::abs((u_0 - MASA_DEFAULT)/MASA_DEFAULT) < threshold)
    {
      cout << "\nMASA ERROR:: Variables not being auto initalized!\n";
      cout << " value is : " << std::abs((u_0 - MASA_DEFAULT)/MASA_DEFAULT) << endl;
      cout << " value is : " << u_0 << endl;
      return 1;
    }

  // now purge values
  masa_purge_default_param<Scalar>();

  // values should be set to default: checking
  u_0 = masa_get_param<Scalar>("u_0");
  if((u_0 - MASA_DEFAULT) > threshold)
    {
      cout << "\nMASA ERROR:: Variables not being purged properly!\n";
      return 1;
    }

#ifdef MASA_EXCEPTIONS
  
  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.
  freopen("/dev/null","w",stdout);
  
  // check sanity check fails
  try
    {
      masa_sanity_check<Scalar>();
    }
  catch(int err) // return one on fatal error
    {
      if((err - 1) > threshold)
	{
	  cout << "regression test failed: masa sanity_check error condition not triggered!\n";
	  return 1; // fail
	}      
    }  
#endif

  // steady as she goes
  return 0;

}
