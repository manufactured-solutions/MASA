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
// c_euler1d.c :program that tests masa against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h> // for MASA_STRICT_REGRESSION
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double fsol_(double x)
{
  return 2*x;
}

double anQ_p (double x,double p_0,double p_x,double a_px,double L)
{
  const double pi = acos(-1);
  double exact_p = p_0 + p_x * cos(a_px * pi * x / L);
  return exact_p;
}
  
double anQ_u (double x,double u_0,double u_x,double a_ux,double L)
{
  const double pi = acos(-1); 
  double exact_u = u_0 + u_x * sin(a_ux * pi * x / L);
  return exact_u;
} 
 
double anQ_rho (double x,double rho_0,double rho_x,double a_rhox,double L)
{ 
  const double pi = acos(-1);  
  double exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  return exact_rho;
}

double SourceQ_rho (
  double x,
  double u_0,
  double u_x,
  double rho_0,
  double rho_x,
  double p_0,
  double p_x,
  double a_px,
  double a_rhox,
  double a_ux,
  double L)
{
  const double pi = acos(-1);  
  double Q_rho;
  double RHO;
  double U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);

  Q_rho = cos(a_ux * pi * x / L) * RHO * a_ux * pi * u_x / L + cos(a_rhox * pi * x / L) * U * a_rhox * pi * rho_x / L;

  return(Q_rho);
}

double SourceQ_u (
  double x,
  double u_0,
  double u_x,
  double rho_0,
  double rho_x,
  double p_0,
  double p_x,
  double a_px,
  double a_rhox,
  double a_ux,
  double L)
{
  const double pi = acos(-1); 
  double Q_u;
  double RHO;
  double U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);

  Q_u = 0.2e1 * cos(a_ux * pi * x / L) * RHO * U * a_ux * pi * u_x / L + cos(a_rhox * pi * x / L) * U * U * a_rhox * pi * rho_x / L - sin(a_px * pi * x / L) * a_px * pi * p_x / L;
  
  return(Q_u);
}

double SourceQ_e (
  double x,
  double u_0,
  double u_x,
  double rho_0,
  double rho_x,
  double p_0,
  double p_x,
  double a_px,
  double a_rhox,
  double a_ux,
  double Gamma,
  double mu,
  double L)
{
  const double pi = acos(-1);  
  double Q_e;
  double RHO;
  double U;
  double P;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  P = p_0 + p_x * cos(a_px * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);

  Q_e = cos(a_rhox * pi * x / L) * pow(U, 0.3e1) * a_rhox * pi * rho_x / L / 0.2e1 + cos(a_ux * pi * x / L) * P * a_ux * pi * u_x * Gamma / L / (Gamma - 0.1e1) + 0.3e1 / 0.2e1 * cos(a_ux * pi * x / L) * RHO * U * U * a_ux * pi * u_x / L - sin(a_px * pi * x / L) * U * a_px * pi * p_x * Gamma / L / (Gamma - 0.1e1);
  return(Q_e);
}

int main()
{
  //variables 
  double u_0;
  double u_x;
  double rho_0;
  double rho_x;
  double p_0;
  double p_x;
  double a_px;
  double a_rhox;
  double a_ux;
  double Gamma;
  double mu;
  double L;

  // parameters
  double x;
  int i;

  //problem size
  int nx    = 200;  // number of points
  int lx    = 10;     // length
  double dx = (double)(lx)/(double)(nx);

  // solutions
  double ufield,ufield2,ufield3;
  double vfield,vfield2,vfield3;
  double efield,efield2,efield3;
  double rho,rho2,rho3;

  double exact_u,exact_u2,exact_u3;
  double exact_v,exact_v2,exact_v3;
  double exact_p,exact_p2,exact_p3;
  double exact_rho,exact_rho2,exact_rho3;

  // initalize
  cmasa_init("euler-test","euler_1d");

  // initialize the default parameters
  cmasa_init_param();

  // get defaults for comparison to source terms
  // get vars
  u_0 = cmasa_get_param("u_0");
  u_x = cmasa_get_param("u_x");

  rho_0 = cmasa_get_param("rho_0");
  rho_x = cmasa_get_param("rho_x");

  p_0 = cmasa_get_param("p_0");
  p_x = cmasa_get_param("p_x");

  a_px = cmasa_get_param("a_px");

  a_rhox = cmasa_get_param("a_rhox");

  a_ux = cmasa_get_param("a_ux");

  Gamma = cmasa_get_param("Gamma");
  mu    = cmasa_get_param("mu");
  L     = cmasa_get_param("L");

  // check that all terms have been initialized
  cmasa_sanity_check();

  // evaluate source terms (1D)
  for(i=0;i<nx;i++)
    {
      x=i*dx;

      //evalulate source terms
      ufield = cmasa_eval_1d_source_rho_u(x);
      efield = cmasa_eval_1d_source_rho_e    (x);
      rho    = cmasa_eval_1d_source_rho  (x);
	
      //evaluate analytical terms
      exact_u   = cmasa_eval_1d_exact_u      (x);
      exact_p   = cmasa_eval_1d_exact_p      (x);
      exact_rho = cmasa_eval_1d_exact_rho     (x);
	
      // get fundamental source term solution
      ufield2   = SourceQ_u  (x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L);
      rho2      = SourceQ_rho(x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L);
      efield2   = SourceQ_e  (x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,Gamma,mu,L);
  
      exact_u2   = anQ_u   (x,u_0,u_x,a_ux,L);
      exact_rho2 = anQ_rho (x,rho_0,rho_x,a_rhox,L);
      exact_p2   = anQ_p   (x,p_0,p_x,a_px,L);

      // test the result is roughly zero
      // choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

      ufield3 = fabs(ufield-ufield2);
      efield3 = fabs(efield-efield2);
      rho3    = fabs(rho-rho2);
      
      exact_u3   = fabs(exact_u-exact_u2);
      exact_rho3 = fabs(exact_rho-exact_rho2);
      exact_p3   = fabs(exact_p-exact_p2);

#else

      ufield3 = fabs(ufield-ufield2)/fabs(ufield2);
      efield3 = fabs(efield-efield2)/fabs(efield2);
      rho3    = fabs(rho-rho2)/fabs(rho2);
      
      exact_u3   = fabs(exact_u-exact_u2)/fabs(exact_u2);
      exact_rho3 = fabs(exact_rho-exact_rho2)/fabs(exact_rho2);
      exact_p3   = fabs(exact_p-exact_p2)/fabs(exact_p2);

#endif

	if(ufield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	    printf("U Field Source Term\n");
	    printf("Threshold Exceeded: %g\n",ufield3);
	    printf("CMASA:              %5.16f\n",ufield);
	    printf("Maple:              %5.16f\n",ufield2);
	    exit(1);
	  }

	if(exact_u3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	    printf("U Field Analytical Term\n");
	    printf("Threshold Exceeded: %g\n",exact_u3);
	    printf("CMASA:              %5.16f\n",exact_u);
	    printf("Maple:              %5.16f\n",exact_u2);
	    exit(1);
	  }

	if(efield3 > threshold)
	  {

	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	    printf("Energy Source Term\n");
	    printf("Threshold Exceeded: %g\n",efield3);
	    printf("CMASA:              %5.16f\n",efield);
	    printf("Maple:              %5.16f\n",efield2);
	    exit(1);
	  }

	if(exact_p3 > threshold)
	  {
	    
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	    printf("P Field Analytical Term\n");
	    exit(1);
	  }

	if(rho3 > threshold)
	  {

	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	    printf("RHO Source Term\n");
	    exit(1);
	  }

	if(exact_rho3 > threshold)
	  {	    
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	    printf("RHO Analytical Term\n");
	    exit(1);
	  }

    } // done interating 

  // tests passed
  return 0;
}
