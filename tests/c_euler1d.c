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

#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double anQ_p (double x,double p_0,double p_x,double a_px,double L)
{
  const double pi = acos(-1);
  double p_an = p_0 + p_x * cos(a_px * pi * x / L);
  return p_an;
}
  
double anQ_u (double x,double u_0,double u_x,double a_ux,double L)
{
  const double pi = acos(-1); 
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L);
  return u_an;
} 
 
double anQ_rho (double x,double rho_0,double rho_x,double a_rhox,double L)
{ 
  const double pi = acos(-1);  
  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  return rho_an;
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
  Q_rho = rho_x * cos(a_rhox * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L)) * a_rhox * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L)) * a_ux * pi / L;
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
  Q_u = -sin(a_px * pi * x / L) * a_px * pi * p_x / L + rho_x * cos(a_rhox * pi * x / L) * pow(u_0 + u_x * sin(a_ux * pi * x / L), 0.2e1) * a_rhox * pi / L + 0.2e1 * u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L)) * (u_0 + u_x * sin(a_ux * pi * x / L)) * a_ux * pi / L;
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
  //const double pi = acos(-1);  
  const double pi = M_PI;
  double Q_e;
  Q_e = cos(a_rhox * pi * x / L) * rho_x * pow(u_0 + u_x * sin(a_ux * pi * x / L), 0.3e1) * a_rhox * pi / L / 0.2e1 + cos(a_ux * pi * x / L) * (p_0 + p_x * cos(a_px * pi * x / L)) * a_ux * pi * u_x * Gamma / L / (Gamma - 0.1e1) - Gamma * p_x * sin(a_px * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L)) * a_px * pi / L / (Gamma - 0.1e1) + 0.3e1 / 0.2e1 * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L)) * pow(u_0 + u_x * sin(a_ux * pi * x / L), 0.2e1) * a_ux * pi * u_x / L;
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
  int nx = 2000;  // number of points
  int lx=10;     // length
  double dx=(double)(lx)/(double)(nx);

  // solutions
  double ufield,ufield2,ufield3;
  double efield,efield2,efield3;
  double rho,rho2;

  double u_an,u_an2,u_an3;
  double p_an,p_an2;
  double rho_an,rho_an2;

  // initalize
  cmasa_init("euler-test","euler_1d");

  // initialize the default parameters
  cmasa_init_param();

  // get defaults for comparison to source terms
  cmasa_get_param("u_0",&u_0);
  cmasa_get_param("u_x",&u_x);

  cmasa_get_param("rho_0",&rho_0);
  cmasa_get_param("rho_x",&rho_x);

  cmasa_get_param("p_0",&p_0);
  cmasa_get_param("p_x",&p_x);

  cmasa_get_param("a_px",&a_px);
  cmasa_get_param("a_rhox",&a_rhox);
  cmasa_get_param("a_ux",&a_ux);

  cmasa_get_param("Gamma",&Gamma);
  cmasa_get_param("mu",&mu);
  cmasa_get_param("L",&L);

  // check that all terms have been initialized
  cmasa_sanity_check();

  // evaluate source terms (1D)
  for(i=0;i<nx;i++)
    {
      x=i*dx;
      
      // ask for masa answer
      cmasa_eval_1d_u_source  (x,&ufield);
      cmasa_eval_1d_e_source  (x,&efield);
      cmasa_eval_1d_rho_source(x,&rho);

      //evaluate analytical terms
      cmasa_eval_1d_u_an        (x,&u_an);
      cmasa_eval_1d_p_an        (x,&p_an);
      cmasa_eval_1d_rho_an      (x,&rho_an);

      // get fundamental source term solution
      ufield2   = SourceQ_u  (x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L);
      rho2      = SourceQ_rho(x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L);
      efield2   = SourceQ_e  (x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,Gamma,mu,L);
  
      u_an2   = anQ_u   (x,u_0,u_x,a_ux,L);
      rho_an2 = anQ_rho (x,rho_0,rho_x,a_rhox,L);
      p_an2   = anQ_p   (x,p_0,p_x,a_px,L);

      // test the result is roughly zero
      ufield3 = fabs(ufield-ufield2);
      efield3 = fabs(efield-efield2);
      rho     = fabs(rho-rho2);

      u_an3   = fabs(u_an-u_an2);
      rho_an = fabs(rho_an-rho_an2);
      p_an   = fabs(p_an-p_an2);
   
      if(ufield3 > threshold)
	{
	  printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	  printf("U Field Source Term\n");
	  printf("Threshold Exceeded: %g\n",ufield3);
	  printf("CMASA:              %5.16f\n",ufield);
	  printf("Maple:              %5.16f\n",ufield2);
	  printf("x:                  %g\n",x);
	  exit(1);
	}

      if(u_an3 > threshold)
	{
	  printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	  printf("U Field Analytical Term\n");
	  printf("Threshold Exceeded: %g\n",u_an3);
	  printf("CMASA:              %5.16f\n",u_an);
	  printf("Maple:              %5.16f\n",u_an2);
	  printf("x:                  %g\n",x);
	  exit(1);
	}

      if(efield3 > threshold)
	{
	  printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	  printf("Energy Source Term\n");
	  printf("Threshold Exceeded: %g\n",efield3);
	  printf("CMASA:              %5.16f\n",efield);
	  printf("Maple:              %5.16f\n",efield2);
	  printf("x:                  %g\n",x);
	  exit(1);
	}

      if(p_an > threshold)
	{
	  printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	  printf("P Field Analytical Term\n");
	  exit(1);
	}
      
      if(rho > threshold)
	{
	  printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	  printf("RHO Source Term\n");
	  exit(1);
	}
      
      if(rho_an > threshold)
	{
	  printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-1d\n");
	  printf("RHO Analytical Term\n");
	  exit(1);
	}

    } // done interating 

  // tests passed
  return 0;
}
