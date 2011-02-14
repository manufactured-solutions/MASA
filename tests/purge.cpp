// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// $Author: nick $
// $Id: misc.cpp 13998 2010-09-23 18:31:43Z nick $
//
// purge.cpp : program that tests masa purge function
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <iostream>
#include <string>

using namespace MASA;
using namespace std;

typedef double Scalar;
Scalar MASA_DEFAULT = -12345.67;

int main()
{
  Scalar u_0;

  // start problem
  masa_init<Scalar>("masa-test","euler_1d");
  
  // values should be set by default: checking
  u_0 = masa_get_param<Scalar>("u_0");
  if(u_0 == MASA_DEFAULT)
    {
      cout << "\nMASA ERROR:: Variables not being auto initalized!\n";
      return 1;
    }

  // now purge values
  masa_purge_default_param<Scalar>();

  // values should be set to default: checking
  u_0 = masa_get_param<Scalar>("u_0");
  if(u_0 != MASA_DEFAULT)
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
      if(err != 1)
	{
	  cout << "regression test failed: masa sanity_check error condition not triggered!\n";
	  return 1; // fail
	}      
    }  
#endif
  

}
