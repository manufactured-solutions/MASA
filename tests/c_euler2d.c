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
// c_euler2d.c :program that tests masa against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h> // for MASA_STRICT_REGRESSION
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double anQ_p (double x,double y,double p_0,double p_x,double p_y,double a_px,double a_py,double L)
{
  const double pi = acos(-1);

  double exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return exact_p;
}
  
double anQ_u (double x,double y,double u_0,double u_x,double u_y,double a_ux,double a_uy,double L)
{
  const double pi = acos(-1);

  double exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return exact_u;
} 
 
double anQ_v (double x,double y,double v_0,double v_x,double v_y,double a_vx,double a_vy,double L)
{
  const double pi = acos(-1);

  double exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return exact_v;
}

double anQ_rho (double x,double y,double rho_0,double rho_x,double rho_y,double a_rhox,double a_rhoy,double L)
{ 
  const double pi = acos(-1);

  double exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return exact_rho;
}

double SourceQ_e (
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double Gamma,
  double mu,
  double L)
{
  double pi = acos(-1);
  double Q_e;
  double RHO;
  double U;
  double V;
  double P;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  P = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);

  Q_e = -a_px * pi * p_x * Gamma * U * sin(a_px * pi * x / L) / (Gamma - 0.1e1) / L + a_py * pi * p_y * Gamma * V * cos(a_py * pi * y / L) / (Gamma - 0.1e1) / L + (U * U + V * V) * a_rhox * pi * rho_x * U * cos(a_rhox * pi * x / L) / L / 0.2e1 - (U * U + V * V) * a_rhoy * pi * rho_y * V * sin(a_rhoy * pi * y / L) / L / 0.2e1 + (0.3e1 * a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * U * U / L / 0.2e1 - (a_uy * u_y * sin(a_uy * pi * y / L) + a_vx * v_x * sin(a_vx * pi * x / L)) * pi * RHO * U * V / L + (a_ux * u_x * cos(a_ux * pi * x / L) + 0.3e1 * a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * V * V / L / 0.2e1 + (a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * Gamma * P / (Gamma - 0.1e1) / L;

  return(Q_e);
}

double SourceQ_u ( // 23 variables
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double L)
{
  double pi = acos(-1);
  double Q_u;
  double RHO;
  double U;
  double V;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);

  Q_u = a_rhox * pi * rho_x * U * U * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * U * V * sin(a_rhoy * pi * y / L) / L - a_uy * pi * u_y * RHO * V * sin(a_uy * pi * y / L) / L - a_px * pi * p_x * sin(a_px * pi * x / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * U / L;

  return(Q_u);
}

double SourceQ_v (
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double L)
{
  double pi = acos(-1);
  double Q_v;

  double RHO;
  double U;
  double V;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);

  Q_v = a_rhox * pi * rho_x * U * V * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * V * sin(a_rhoy * pi * y / L) / L - a_vx * pi * v_x * RHO * U * sin(a_vx * pi * x / L) / L + a_py * pi * p_y * cos(a_py * pi * y / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * V / L;

  return(Q_v);
}

double SourceQ_rho(
  double x,
  double y,
  double u_0,
  double u_x,
  double u_y,
  double v_0,
  double v_x,
  double v_y,
  double rho_0,
  double rho_x,
  double rho_y,
  double p_0,
  double p_x,
  double p_y,
  double a_px,
  double a_py,
  double a_rhox,
  double a_rhoy,
  double a_ux,
  double a_uy,
  double a_vx,
  double a_vy,
  double L)
{
  double pi = acos(-1);
  double Q_rho;

  double RHO;
  double U;
  double V;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);

  Q_rho = a_rhox * pi * rho_x * U * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * sin(a_rhoy * pi * y / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO / L;


  return(Q_rho);
}

int main()

{  //variables
  double u_0;
  double u_x;
  double u_y;
  double v_0;
  double v_x;
  double v_y;
  double rho_0;
  double rho_x;
  double rho_y;
  double p_0;
  double p_x;
  double p_y;
  double a_px;
  double a_py;
  double a_rhox;
  double a_rhoy;
  double a_ux;
  double a_uy;
  double a_vx;
  double a_vy;
  double Gamma;
  double mu;
  double L;

  // parameters
  double x;
  double y;
  int i,j,k;

  // solutions -- efield is MASA term, efield2 is maple, efield3 is abs error between them
  double ufield,ufield2,ufield3;
  double vfield,vfield2,vfield3;
  double efield,efield2,efield3;
  double rho,rho2,rho3;

  double exact_u,exact_u2,exact_u3;
  double exact_v,exact_v2,exact_v3;
  double exact_p,exact_p2,exact_p3;
  double exact_rho,exact_rho2,exact_rho3;

  // initalize
  int nx = 125;  // number of points
  int ny = 68;  
  int lx=3;     // length
  int ly=1; 
  
  double dx=(double)(lx)/(double)(nx);
  double dy=(double)(ly)/(double)(ny);

  masa_init("euler-test","euler_2d");

  // set params
  masa_init_param();
  
  // get vars
  u_0 = masa_get_param("u_0");
  u_x = masa_get_param("u_x");
  u_y = masa_get_param("u_y");

  v_0 = masa_get_param("v_0");
  v_x = masa_get_param("v_x");
  v_y = masa_get_param("v_y");

  rho_0 = masa_get_param("rho_0");
  rho_x = masa_get_param("rho_x");
  rho_y = masa_get_param("rho_y");

  p_0 = masa_get_param("p_0");
  p_x = masa_get_param("p_x");
  p_y = masa_get_param("p_y");

  a_px = masa_get_param("a_px");
  a_py = masa_get_param("a_py");

  a_rhox = masa_get_param("a_rhox");
  a_rhoy = masa_get_param("a_rhoy");

  a_ux = masa_get_param("a_ux");
  a_uy = masa_get_param("a_uy");

  a_vx = masa_get_param("a_vx");
  a_vy = masa_get_param("a_vy");

  Gamma = masa_get_param("Gamma");
  mu    = masa_get_param("mu");
  L     = masa_get_param("L");

  // check that all terms have been initialized
  masa_sanity_check();

  // evaluate source terms (2D)
  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)    
      {
	x=i*dx;
	y=j*dy;
	
	//evalulate source terms
	ufield = masa_eval_2d_source_rho_u(x,y);
	vfield = masa_eval_2d_source_rho_v(x,y);
	efield = masa_eval_2d_source_rho_e(x,y);
	rho    = masa_eval_2d_source_rho  (x,y);
	
	//evaluate analytical terms
	exact_u   = masa_eval_2d_exact_u  (x,y);
	exact_v   = masa_eval_2d_exact_v  (x,y);
	exact_p   = masa_eval_2d_exact_p  (x,y);
	exact_rho = masa_eval_2d_exact_rho(x,y);
	
	// check against maple
	ufield2 = SourceQ_u   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);
	vfield2 = SourceQ_v   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);
	rho2    = SourceQ_rho (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);  
	efield2 = SourceQ_e   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,Gamma,mu,L);
	
	exact_u2   = anQ_u   (x,y,u_0,u_x,u_y,a_ux,a_uy,L);
	exact_v2   = anQ_v   (x,y,v_0,v_x,v_y,a_vx,a_vy,L);
	exact_rho2 = anQ_rho (x,y,rho_0,rho_x,rho_y,a_rhox,a_rhoy,L);
	exact_p2   = anQ_p   (x,y,p_0,p_x,p_y,a_px,a_py,L);

	// test the result is roughly zero
	// choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	ufield3 = fabs(ufield-ufield2);
	vfield3 = fabs(vfield-vfield2);
	efield3 = fabs(efield-efield2);
	rho3    = fabs(rho-rho2);

	exact_u3   = fabs(exact_u-exact_u2);
	exact_v3   = fabs(exact_v-exact_v2);
	exact_rho3 = fabs(exact_rho-exact_rho2);
	exact_p3   = fabs(exact_p-exact_p2);

#else

	ufield3 = fabs(ufield-ufield2)/fabs(ufield2);
	vfield3 = fabs(vfield-vfield2)/fabs(vfield2);
	efield3 = fabs(efield-efield2)/fabs(efield2);
	rho3    = fabs(rho-rho2)/fabs(rho2);

	exact_u3   = fabs(exact_u-exact_u2)/fabs(exact_u2);
	exact_v3   = fabs(exact_v-exact_v2)/fabs(exact_v2);
	exact_rho3 = fabs(exact_rho-exact_rho2)/fabs(exact_rho2);
	exact_p3   = fabs(exact_p-exact_p2)/fabs(exact_p2);

#endif

	if(ufield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("U Field Source Term\n");
	    printf("Threshold Exceeded: %g\n",ufield3);
	    printf("MASA:              %5.16f\n",ufield);
	    printf("Maple:              %5.16f\n",ufield2);
	    printf("x,y:                %g %g\n",x,y);	   
	    exit(1);
	  }

	if(exact_u3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("U Field Analytical Term\n");
	    printf("Threshold Exceeded: %g\n",exact_u3);
	    printf("MASA:              %5.16f\n",exact_u);
	    printf("Maple:              %5.16f\n",exact_u2);
	    printf("x,y:                %g %g\n",x,y);	   
	    exit(1);
	  }

	if(vfield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("V Field Source Term\n");
	    printf("Threshold Exceeded: %g\n",vfield3);
	    printf("MASA:              %5.16f\n",vfield);
	    printf("Maple:              %5.16f\n",vfield2);
	    printf("x,y:                %g %g\n",x,y);	   
	    exit(1);
	  }

	if(exact_v3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("V Field Analytical Term\n");
	    exit(1);
	  }

	if(efield3 > threshold)
	  {

	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("Energy Source Term\n");
	    printf("Threshold Exceeded: %g\n",efield3);
	    printf("MASA:              %5.16f\n",efield);
	    printf("Maple:              %5.16f\n",efield2);
	    printf("x,y:                %g %g\n",x,y);	   
	    exit(1);
	  }

	if(exact_p3 > threshold)
	  {
	    
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("P Field Analytical Term\n");
	    exit(1);
	  }

	if(rho3 > threshold)
	  {

	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("RHO Source Term\n");
	    exit(1);
	  }

	if(exact_rho3 > threshold)
	  {	    
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("RHO Analytical Term\n");
	    exit(1);
	  }

      } // done iterating

  // tests passed
  return 0;
}
