// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
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
// $Id: rans_sa.cpp 22699 2011-08-01 03:22:24Z nick $
//
// sod.cpp: this is an example of the API used for 
//          calling the sod shock tube solution
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>

using namespace MASA;
using namespace std;

typedef double Scalar;

int main()
{

  int i,j;
  int err = 0;
  Scalar x;
  Scalar t;
  Scalar out;

  int nx = 200;  // number of points
  int lx=10;     // length
  int nt = 22;
  int lt=4;   
  Scalar dx=Scalar(lx)/Scalar(nx);
  Scalar dt=Scalar(lt)/Scalar(nt);

  err  = masa_init<Scalar>("sod-test","sod_1d");
  err += masa_sanity_check<Scalar>();
  
  for(i=0;i<nx;i++)
    {
      for(j=0;j<nt;j++)
	{
	  x=i*dx;
	  t=j*dt;
	  
	  out = masa_eval_source_rho  <Scalar>(x,t);
	  out = masa_eval_source_rho_u<Scalar>(x,t);	
	  
	}
    } //done iterating

  return err;

}// end program
