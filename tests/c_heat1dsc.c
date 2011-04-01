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
// $Author$
// $Id$
//
// c_heat1dsc.c: program that tests heat equation steady, constant
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <tests.h>

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double SourceQ_t_1d(double x, double A_x, double k_0)
{
  double Q_T = A_x * A_x * k_0 * cos(A_x * x);
  return Q_T;
}

double Source_t_1d_exact(double A_x,double x)
{
  double exact_t;
  exact_t = cos(A_x * x);
  return exact_t;
}

double SourceQ_t_2d (
  double x,
  double y,
  double A_x,
  double B_y,
  double k_0)
{
  double Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * (A_x * A_x + B_y * B_y);
  return Q_T;
}

double SourceQ_t_3d (
  double x,
  double y,
  double z,
  double A_x,
  double B_y,
  double C_z,
  double k_0)
{
  double Q_T = k_0 * cos(A_x * x) * cos(B_y * y) * cos(C_z * z) * (A_x * A_x + B_y * B_y + C_z * C_z);
  return Q_T;
}


int main()
{
  int i;

  double tfield,tfield2,tfield3;
  double exact_t,exact_t2,exact_t3;
  double x;

  double A_x;
  double k_0;

  //problem size
  int nx = 200;  // number of points
  int lx=10;     // length
  double dx = (double)lx/(double)nx;

  // initalize everyone
  masa_init("temp-test-1d","heateq_1d_steady_const");
  masa_init_param();
  masa_sanity_check();

  A_x = masa_get_param("A_x");
  k_0 = masa_get_param("k_0"); 

  // evaluate source terms (1D)
  for(i=0;i<nx;i++)
    {
      x=i*dx;
      
      //evalulate source terms
      tfield = masa_eval_1d_source_t(x);
      
      //evaluate analytical terms
      exact_t   = masa_eval_1d_exact_t(x);
	
      // get fundamental source term solution
      tfield2   = SourceQ_t_1d  (x,A_x,k_0);
      exact_t2     = Source_t_1d_exact(x,A_x);

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
	    exit(1);
	  }

	if(exact_t3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Heat Equation Steady-2d\n");
	    printf("U Field Analytical Term\n");
	    printf("Threshold Exceeded: %g\n",exact_t3);
	    printf("MASA:              %5.16f\n",exact_t);
	    printf("Maple:              %5.16f\n",exact_t2);
	    exit(1);
	  }
    } // done iterating


  masa_init("temp-test-2d","heateq_2d_steady_const");
  masa_init_param();
  masa_sanity_check();  

  masa_init("temp-test-3d","heateq_3d_steady_const");
  masa_init_param();
  masa_sanity_check();

  return 0; // steady as she goes

}
