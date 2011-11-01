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
// $Id: rans_sa.cpp 18184 2011-03-01 20:09:57Z nick $
//
// fans_sa_transient_d_finite: this is an example of the API used for 
//                             calling the spelart alamaras Favre Average 
//                             Navier Stokes (FANS) model
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>

using namespace MASA;
typedef double Scalar;

int main()
{

  Scalar ufield;
  Scalar vfield;
  Scalar efield;
  Scalar rho;
  Scalar nu;

  // declarations
  Scalar tempx;
  Scalar tempy;
  Scalar tempt;

  Scalar exact_u;
  Scalar exact_v;
  Scalar exact_p;
  Scalar exact_rho;
  Scalar exact_nu;

  //error handing
  int err = 0;

  //problem size

  //problem size
  Scalar lx,ly,lt;
  Scalar dx,dy,dt;
  int nx,ny,nt;

  // initialize
  nx = 10;  // number of points
  ny = 10;  
  nt = 10;  
  
  lx=1;     // length
  ly=1; 
  lt=3;

  dx=(Scalar)lx/(Scalar)nx;
  dy=(Scalar)ly/(Scalar)ny;
  dt=(Scalar)lt/(Scalar)nt;
  
  // initialize the problem 
  err = masa_init<Scalar>("sa example","fans_sa_steady_wall_bounded");

  // call the sanity check routine 
  // (tests that all variables have been initialized)
  err = masa_sanity_check<Scalar>();

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(int i=0;i<nx;i++)
    {
      for(int j=0;j<nx;j++)
	{  
	  tempx=i*dx;
	  tempy=j*dy;
	  
	  // evaluate source terms
	  ufield = masa_eval_source_rho_u<Scalar>  (tempx,tempy);
	  vfield = masa_eval_source_rho_v<Scalar>  (tempx,tempy);
	  efield = masa_eval_source_rho_e<Scalar>  (tempx,tempy);
	  rho    = masa_eval_source_rho  <Scalar>  (tempx,tempy);
	  nu     = masa_eval_source_nu   <Scalar>  (tempx,tempy);
	
	  //evaluate analytical solution
	  exact_u   = masa_eval_exact_u  <Scalar>   (tempx,tempy);
	  exact_v   = masa_eval_exact_v  <Scalar>   (tempx,tempy);
	  exact_p   = masa_eval_exact_p  <Scalar>   (tempx,tempy);
	  exact_rho = masa_eval_exact_rho<Scalar>   (tempx,tempy);
	  exact_nu  = masa_eval_exact_nu <Scalar>   (tempx,tempy);

	}
    } // done with spatial loop
  return err;

}// end program
