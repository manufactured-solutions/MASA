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
// $Author: 
// $Id$
//
// c_heat2dsc.c: program that tests heat equation steady, constant
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include<tests.h> 

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double SourceQ_t(double x,double y,double A_x,double B_y,double k_0)
{
  double Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * (A_x * A_x + B_y * B_y);
  return Q_T;
}

double Source_exact_t(double x,double y,double A_x,double B_y)
{
  double exact_t;
  exact_t = cos(A_x * x) * cos(B_y * y);
  return exact_t;
}

int main()
{
  int i;
  int j;

  double tfield,tfield2,tfield3;
  double exact_t,exact_t2,exact_t3;
  double x;
  double y;

  double A_x;
  double k_0;

  double B_y;

  //problem size
  int nx = 100;  // number of points
  int lx = 10;     // length
  double dx = (double)lx/(double)nx;

  int ny  = 100;  // number of points
  int ly  = 10;     // length
  double dy = (double)ly/(double)ny;

  // initalize everyone
  masa_init("temp-test-1d","heateq_2d_steady_const");
  masa_init_param();
  masa_sanity_check();

  A_x = masa_get_param("A_x");
  k_0 = masa_get_param("k_0"); 
  B_y = masa_get_param("B_y"); 

  // evaluate source terms (1D)
  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)      
      {
      x=i*dx;
      y=j*dy;
      
      //evalulate source terms
      tfield = masa_eval_2d_source_t(x,y);
      
      //evaluate analytical terms
      exact_t   = masa_eval_2d_exact_t(x,y);
	
      // get fundamental source term solution
      tfield2   = SourceQ_t(x,y,A_x,B_y,k_0);
      exact_t2     = Source_exact_t(x,y,A_x,B_y);

      // test the result is roughly zero
      // choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

      tfield3 = fabs(tfield-tfield2);
      exact_t3   = fabs(exact_t-exact_t2);

#else

      tfield3 = fabs(tfield-tfield2)/fabs(tfield2);
      exact_t3   = fabs(exact_t-exact_t2)/fabs(tfield2);

#endif

	if(tfield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Heat Equation Steady-2d\n");
	    printf("U Field Source Term\n");
	    printf("Threshold Exceeded: %g\n",tfield3);
	    printf("MASA:              %5.16f\n",tfield);
	    printf("Maple:              %5.16f\n",tfield2);
	    printf("@ x,y:              %5.16f %5.16f\n",x,y);
	    exit(1);
	  }

	if(exact_t3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Heat Equation Steady-2d\n");
	    printf("U Field Analytical Term\n");
	    printf("Threshold Exceeded: %g\n",exact_t3);
	    printf("MASA:              %5.16f\n",exact_t);
	    printf("Maple:              %5.16f\n",exact_t2);
	    printf("@ x,y:              %5.16f %5.16f\n",x,y);
	    exit(1);
	  }
      } // done iterating

  return 0; // steady as she goes

}
