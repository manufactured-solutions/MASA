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
// c_euler2d.c :program that tests masa against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double anQ_p (double x,double y,double p_0,double p_x,double p_y,double a_px,double a_py,double L)
{
  const double pi = acos(-1);

  double p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return p_an;
}
  
double anQ_u (double x,double y,double u_0,double u_x,double u_y,double a_ux,double a_uy,double L)
{
  const double pi = acos(-1);

  double u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return u_an;
} 
 
double anQ_v (double x,double y,double v_0,double v_x,double v_y,double a_vx,double a_vy,double L)
{
  const double pi = acos(-1);

  double v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return v_an;
}

double anQ_rho (double x,double y,double rho_0,double rho_x,double rho_y,double a_rhox,double a_rhoy,double L)
{ 
  const double pi = acos(-1);

  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return rho_an;
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
  double PI = acos(-1);
  double Q_e;
  Q_e = -Gamma * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * a_px * p_x * PI * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + Gamma * (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * a_py * p_y * PI * cos(a_py * PI * y / L) / (Gamma - 0.1e1) / L + a_rhox * PI * rho_x * cos(a_rhox * PI * x / L) * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * (pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) / L / 0.2e1 - a_rhoy * PI * rho_y * sin(a_rhoy * PI * y / L) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * (pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) / L / 0.2e1 + (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * a_ux * PI * u_x * cos(x * a_ux * PI / L) * Gamma / (Gamma - 0.1e1) / L + (pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + 0.3e1 * pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_ux * PI * u_x * cos(x * a_ux * PI / L) / L / 0.2e1 + (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * a_vy * PI * v_y * cos(y * a_vy * PI / L) * Gamma / (Gamma - 0.1e1) / L + (0.3e1 * pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_vy * PI * v_y * cos(y * a_vy * PI / L) / L / 0.2e1 - (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * PI * a_uy * u_y * sin(a_uy * PI * y / L) / L - (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * PI * a_vx * v_x * sin(a_vx * PI * x / L) / L;
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
  double PI = acos(-1);
  double Q_u;
  Q_u = -p_x * sin(a_px * PI * x / L) * a_px * PI / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L), 0.2e1) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhoy * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_ux * PI / L - u_y * sin(a_uy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_uy * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vy * PI / L;
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
  double PI = acos(-1);
  double Q_v;
  Q_v = p_y * cos(a_py * PI * y / L) * a_py * PI / L + rho_x * cos(a_rhox * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhox * PI / L - sin(a_rhoy * PI * y / L) * rho_y * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L), 0.2e1) * a_rhoy * PI / L + cos(a_ux * PI * x / L) * u_x * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_ux * PI / L - sin(a_vx * PI * x / L) * v_x * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vx * PI / L + 0.2e1 * cos(a_vy * PI * y / L) * v_y * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_vy * PI / L;
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
  double PI = acos(-1);
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_rhoy * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * a_ux * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * a_vy * PI / L;
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
  double rho,rho2;

  double u_an,u_an2,u_an3;
  double v_an,v_an2,v_an3;
  double p_an,p_an2;
  double rho_an,rho_an2;

  // initalize
  int nx = 125;  // number of points
  int ny = 68;  
  int lx=3;     // length
  int ly=1; 
  
  double dx=(double)(lx)/(double)(nx);
  double dy=(double)(ly)/(double)(ny);

  cmasa_init("euler-test","euler_2d");

  // set params
  cmasa_init_param();
  
  // get vars
  cmasa_get_param("u_0",&u_0);
  cmasa_get_param("u_x",&u_x);
  cmasa_get_param("u_y",&u_y);
  cmasa_get_param("v_0",&v_0);
  cmasa_get_param("v_x",&v_x);
  cmasa_get_param("v_y",&v_y);

  cmasa_get_param("rho_0",&rho_0);
  cmasa_get_param("rho_x",&rho_x);
  cmasa_get_param("rho_y",&rho_y);

  cmasa_get_param("p_0",&p_0);
  cmasa_get_param("p_x",&p_x);
  cmasa_get_param("p_y",&p_y);

  cmasa_get_param("a_px",&a_px);
  cmasa_get_param("a_py",&a_py);

  cmasa_get_param("a_rhox",&a_rhox);
  cmasa_get_param("a_rhoy",&a_rhoy);

  cmasa_get_param("a_ux",&a_ux);
  cmasa_get_param("a_uy",&a_uy);

  cmasa_get_param("a_vx",&a_vx);
  cmasa_get_param("a_vy",&a_vy);

  cmasa_get_param("Gamma",&Gamma);
  cmasa_get_param("mu",&mu);
  cmasa_get_param("L",&L);

  // check that all terms have been initialized
  cmasa_sanity_check();

  // evaluate source terms (2D)
  for(i=0;i<nx;i++)
    for(j=0;j<ny;j++)    
      {
	x=i*dx;
	y=j*dy;
	
	//evalulate source terms
	cmasa_eval_2d_u_source  (x,y,&ufield);
	cmasa_eval_2d_v_source  (x,y,&vfield);
	cmasa_eval_2d_e_source  (x,y,&efield);
	cmasa_eval_2d_rho_source(x,y,&rho);

	//evaluate analytical terms
	cmasa_eval_2d_u_an        (x,y,&u_an);
	cmasa_eval_2d_v_an        (x,y,&v_an);
	cmasa_eval_2d_p_an        (x,y,&p_an);
	cmasa_eval_2d_rho_an      (x,y,&rho_an);
	  
	// check against maple
	ufield2 = SourceQ_u   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);
	vfield2 = SourceQ_v   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);
	rho2    = SourceQ_rho (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);  
	efield2 = SourceQ_e   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,Gamma,mu,L);
	
	u_an2   = anQ_u   (x,y,u_0,u_x,u_y,a_ux,a_uy,L);
	v_an2   = anQ_v   (x,y,v_0,v_x,v_y,a_vx,a_vy,L);
	rho_an2 = anQ_rho (x,y,rho_0,rho_x,rho_y,a_rhox,a_rhoy,L);
	p_an2   = anQ_p   (x,y,p_0,p_x,p_y,a_px,a_py,L);
	
	// test the result is roughly zero
	ufield3 = fabs(ufield-ufield2);
	vfield3 = fabs(vfield-vfield2);
	efield3 = fabs(efield-efield2);
	rho     = fabs(rho-rho2);
	
	u_an3  = fabs(u_an-u_an2);
	v_an3  = fabs(v_an-v_an2);
	rho_an = fabs(rho_an-rho_an2);
	p_an   = fabs(p_an-p_an2);

	if(ufield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("U Field Source Term\n");
	    printf("Threshold Exceeded: %g\n",ufield3);
	    printf("CMASA:              %5.16f\n",ufield);
	    printf("Maple:              %5.16f\n",ufield2);
	    printf("x,y:                %g %g\n",x,y);	   
	    exit(1);
	  }

	if(u_an3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("U Field Analytical Term\n");
	    printf("Threshold Exceeded: %g\n",u_an3);
	    printf("CMASA:              %5.16f\n",u_an);
	    printf("Maple:              %5.16f\n",u_an2);
	    printf("x,y:                %g %g\n",x,y);	   
	    exit(1);
	  }

	if(vfield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("V Field Source Term\n");
	    printf("Threshold Exceeded: %g\n",vfield3);
	    printf("CMASA:              %5.16f\n",vfield);
	    printf("Maple:              %5.16f\n",vfield2);
	    printf("x,y:                %g %g\n",x,y);	   
	    exit(1);
	  }

	if(v_an3 > threshold)
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
	    printf("CMASA:              %5.16f\n",efield);
	    printf("Maple:              %5.16f\n",efield2);
	    printf("x,y:                %g %g\n",x,y);	   
	    exit(1);
	  }

	if(p_an > threshold)
	  {
	    
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("P Field Analytical Term\n");
	    exit(1);
	  }

	if(rho > threshold)
	  {

	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("RHO Source Term\n");
	    exit(1);
	  }

	if(rho_an > threshold)
	  {	    
	    printf("\nMASA REGRESSION TEST FAILED: C-binding Euler-2d\n");
	    printf("RHO Analytical Term\n");
	    exit(1);
	  }

      } // done iterating

  // tests passed
  return 0;
}
