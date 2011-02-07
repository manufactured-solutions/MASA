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
// sod.cpp : test sod mms
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace MASA;
using namespace std;

template<typename Scalar>
int run_regression()
{
  int i,j,k;
  Scalar out;
  Scalar x;
  Scalar t;

  int nx = 200;  // number of points
  int lx=10;     // length
  int nt = 22;
  int lt=4;   
  Scalar dx=Scalar(lx)/Scalar(nx);
  Scalar dt=Scalar(lt)/Scalar(nt);

  masa_init<Scalar>("sod-test","sod_1d");
  masa_init_param<Scalar>();
  masa_sanity_check<Scalar>();
  
  for(i=0;i<nx;i++)
    for(j=0;j<nt;j++)
      {
	x=i*dx;
	t=j*dt;
	
	out = masa_eval_source_rho  <Scalar>(x,t);
	out = masa_eval_source_rho_u<Scalar>(x,t);
	//out = masa_eval_p_source    <Scalar>(x,t);
      } //done iterating

  // below are barely a regression test: 
  // mostly 'touching' code for coverage
  
  // test rarefaction wave before origin
  x = -1;
  t =  1;
  out = masa_eval_source_rho<Scalar>(x,t);
  out = masa_eval_source_rho_u<Scalar>(x,t);

  // touching sod error condition if root not bracketed in rtbis
  x = 1;
  t = 2;

  // nonphysical gamma, but useful to test coverage in rtbis
  masa_set_param<Scalar>("Gamma",-1.0);
  out = masa_eval_source_rho<Scalar>(x,t);

#ifdef MASA_EXCEPTIONS
  
  masa_set_param<Scalar>("Gamma",1.2);

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.
  freopen("/dev/null","w",stdout);

  try
    {
      out = masa_eval_source_t<Scalar>(x);
    }
  catch(int err) // return one on fatal error
    {
      if(err != 1)
	{
	  cout << "regression test failed: sod rtbis root not bracketed error condition not triggered!\n";
	  return 1; // fail
	}      
    }

  try
    {
      out = masa_eval_source_t<Scalar>(x,t);
    }
  catch(int err) // return one on fatal error
    {
      if(err != 1)
	{
	  cout << "regression test failed: sod rtbis `too many bisections` condition not triggered!\n";
	  return 1; // fail
	}      
    }

#endif 
  
  return 0;

}// end program

int main()
{
  int err=0;

  err += run_regression<double>();
  err += run_regression<long double>();

  return err;
}
