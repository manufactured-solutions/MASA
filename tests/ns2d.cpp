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
// ns2d.cpp: program that tests navier-stokes-2d against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h> 
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
using namespace MASA;

typedef double Scalar;

const Scalar pi = acos(-1);
const Scalar threshold = 1.0e-15; // should be small enough to catch any obvious problems

Scalar nancheck(Scalar x)
{
  if(isnan(x))
    {
      cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
  return 1;
}

Scalar anQ_p (Scalar x,Scalar y,Scalar p_0,Scalar p_x,Scalar p_y,Scalar a_px,Scalar a_py,Scalar L)
{
  Scalar p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return p_an;
}
  
Scalar anQ_u (Scalar x,Scalar y,Scalar u_0,Scalar u_x,Scalar u_y,Scalar a_ux,Scalar a_uy,Scalar L)
{
  Scalar u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return u_an;
} 
 
Scalar anQ_v (Scalar x,Scalar y,Scalar v_0,Scalar v_x,Scalar v_y,Scalar a_vx,Scalar a_vy,Scalar L)
{
  Scalar v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return v_an;
}

Scalar anQ_rho (Scalar x,Scalar y,Scalar rho_0,Scalar rho_x,Scalar rho_y,Scalar a_rhox,Scalar a_rhoy,Scalar L)
{ 
  Scalar rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return rho_an;
}

Scalar SourceQ_e (
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar);

Scalar SourceQ_u (
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar);

Scalar SourceQ_v (
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar);

Scalar SourceQ_rho( 
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar);

int main()
{
  //variables
  Scalar u_0;
  Scalar u_x;
  Scalar u_y;
  Scalar v_0;
  Scalar v_x;
  Scalar v_y;
  Scalar rho_0;
  Scalar rho_x;
  Scalar rho_y;
  Scalar p_0;
  Scalar p_x;
  Scalar p_y;
  Scalar a_px;
  Scalar a_py;
  Scalar a_rhox;
  Scalar a_rhoy;
  Scalar a_ux;
  Scalar a_uy;
  Scalar a_vx;
  Scalar a_vy;
  Scalar Gamma;
  Scalar mu;
  Scalar L;
  Scalar R;
  Scalar k;

  // parameters
  Scalar x;
  Scalar y;
  Scalar z;

  // solutions
  Scalar u_an,u_an2,u_an3;
  Scalar v_an,v_an2,v_an3;
  Scalar w_an,w_an2,w_an3;
  Scalar p_an,p_an2,p_an3;
  Scalar rho_an,rho_an2,rho_an3;
  Scalar gradx,grady,gradz,gradp,gradrho;

  Scalar ufield,ufield2,ufield3;
  Scalar vfield,vfield2,vfield3;
  Scalar efield,efield2,efield3;
  Scalar rho,rho2,rho3;

  // initalize
  int err;
  int nx = 10;  // number of points
  int ny = 8;  
  int lx=2;     // length
  int ly=1; 
  
  Scalar dx=Scalar(lx)/Scalar(nx);
  Scalar dy=Scalar(ly)/Scalar(ny);

  masa_init<Scalar>("navier-stokes-test","navierstokes_2d_compressible");

  // set params
  masa_init_param<Scalar>();
  
  // get vars for comparison
  u_0 = masa_get_param<Scalar>("u_0");
  u_x = masa_get_param<Scalar>("u_x");
  u_y = masa_get_param<Scalar>("u_y");

  v_0 = masa_get_param<Scalar>("v_0");
  v_x = masa_get_param<Scalar>("v_x");
  v_y = masa_get_param<Scalar>("v_y");

  rho_0 = masa_get_param<Scalar>("rho_0");
  rho_x = masa_get_param<Scalar>("rho_x");
  rho_y = masa_get_param<Scalar>("rho_y");

  p_0 = masa_get_param<Scalar>("p_0");
  p_x = masa_get_param<Scalar>("p_x");
  p_y = masa_get_param<Scalar>("p_y");

  a_px = masa_get_param<Scalar>("a_px");
  a_py = masa_get_param<Scalar>("a_py");

  a_rhox = masa_get_param<Scalar>("a_rhox");
  a_rhoy = masa_get_param<Scalar>("a_rhoy");

  a_ux = masa_get_param<Scalar>("a_ux");
  a_uy = masa_get_param<Scalar>("a_uy");

  a_vx = masa_get_param<Scalar>("a_vx");
  a_vy = masa_get_param<Scalar>("a_vy");

  Gamma = masa_get_param<Scalar>("Gamma");
  mu    = masa_get_param<Scalar>("mu");
  L     = masa_get_param<Scalar>("L");

  R = masa_get_param<Scalar>("R");
  k = masa_get_param<Scalar>("k");

  // check that all terms have been initialized
  err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }  

  // evaluate source terms (2D)
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)    
      {
	x=i*dx;
	y=j*dy;
	
	//evalulate source terms
	ufield = masa_eval_u_source  <Scalar>(x,y);
	vfield = masa_eval_v_source  <Scalar>(x,y);
	efield = masa_eval_e_source  <Scalar>(x,y);
	rho    = masa_eval_rho_source<Scalar>(x,y);
	
	//evaluate analytical terms
	u_an = masa_eval_u_an        <Scalar>(x,y);
	v_an = masa_eval_v_an        <Scalar>(x,y);
	p_an = masa_eval_p_an        <Scalar>(x,y);
	rho_an = masa_eval_rho_an    <Scalar>(x,y);
	
	// eval gradient terms
	gradx = masa_eval_2d_grad_u<Scalar>(x,y,1);
	grady = masa_eval_2d_grad_u<Scalar>(x,y,2);		

	gradx = masa_eval_2d_grad_v<Scalar>(x,y,1);
	grady = masa_eval_2d_grad_v<Scalar>(x,y,2);		

	gradx = masa_eval_2d_grad_p<Scalar>(x,y,1);
	grady = masa_eval_2d_grad_p<Scalar>(x,y,2);		
  
	gradx = masa_eval_2d_grad_rho<Scalar>(x,y,1);
	grady = masa_eval_2d_grad_rho<Scalar>(x,y,2);		

	// check against maple
	ufield2 = SourceQ_u   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);
	vfield2 = SourceQ_v   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);
	rho2    = SourceQ_rho (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);  
	efield2 = SourceQ_e   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,Gamma,mu,L,R,k);
	
	u_an2   = anQ_u   (x,y,u_0,u_x,u_y,a_ux,a_uy,L);
	v_an2   = anQ_v   (x,y,v_0,v_x,v_y,a_vx,a_vy,L);
	rho_an2 = anQ_rho (x,y,rho_0,rho_x,rho_y,a_rhox,a_rhoy,L);
	p_an2   = anQ_p   (x,y,p_0,p_x,p_y,a_px,a_py,L);

	// test the result is roughly zero
	// choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	ufield3 = fabs(ufield-ufield2);
	vfield3 = fabs(vfield-vfield2);
	efield3 = fabs(efield-efield2);
	rho3    = fabs(rho-rho2);

	u_an3   = fabs(u_an-u_an2);
	v_an3   = fabs(v_an-v_an2);
	rho_an3 = fabs(rho_an-rho_an2);
	p_an3   = fabs(p_an-p_an2);

#else

	ufield3 = fabs(ufield-ufield2)/fabs(ufield2);
	vfield3 = fabs(vfield-vfield2)/fabs(vfield2);
	efield3 = fabs(efield-efield2)/fabs(efield2);
	rho3    = fabs(rho-rho2)/fabs(rho2);

	u_an3   = fabs(u_an-u_an2)/fabs(u_an2);
	v_an3   = fabs(v_an-v_an2)/fabs(v_an2);
	rho_an3 = fabs(rho_an-rho_an2)/fabs(rho_an2);
	p_an3   = fabs(p_an-p_an2)/fabs(p_an2);

#endif
	
	nancheck(ufield3);
	nancheck(vfield3);
	nancheck(efield3);
	nancheck(rho3);
	
	nancheck(u_an3);
	nancheck(v_an3);
	nancheck(rho_an3);
	nancheck(p_an3);

	if(ufield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n");
	    printf("U Field Source Term\n");
	    exit(1);
	  }

	if(u_an3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n");
	    printf("U Field Analytical Term\n");
	    exit(1);
	  }

	if(vfield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n");
	    printf("V Field Source Term\n");
	    printf("(Relative) Threshold Exceeded: %g\n",vfield3);
	    printf("CMASA:              %5.16f\n",vfield);
	    printf("Maple:              %5.16f\n",vfield2);
	    printf("x,y:                %g %g\n",x,y);

	    exit(1);
	  }

	if(v_an3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n");
	    printf("V Field Analytical Term\n");
	    exit(1);
	  }

	if(efield3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n");
	    printf("E Field Source Term\n");
	    printf("(Relative) Threshold Exceeded: %g\n",efield3);
	    printf("CMASA:              %5.16f\n",efield);
	    printf("Maple:              %5.16f\n",efield2);
	    printf("x,y:                %g %g\n",x,y);
	    exit(1);
	  }

	if(p_an3 > threshold)
	  {
	    
	    printf("\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n");
	    printf("P Field Analytical Term\n");
	    exit(1);
	  }

	if(rho3 > threshold)
	  {
	    printf("\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n");
	    printf("RHO Field Source Term\n");
	    exit(1);
	  }

	if(rho_an3 > threshold)
	  {	    
	    printf("\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n");
	    printf("RHO Analytical Term\n");
	    exit(1);
	  }

      } // done iterating

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.
  freopen("/dev/null","w",stdout);  
  
  // test gradient error terms
  Scalar derr = masa_eval_2d_grad_u<Scalar>(0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_2d_grad_u<Scalar>(0,0,3);
  if(derr != -1)
    {
      cout << "MASA :: gradient (3) error condition failed!\n";
      exit(1);
    }

  // v
  derr = masa_eval_2d_grad_v<Scalar>(0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: v gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_2d_grad_v<Scalar>(0,0,3);
  if(derr != -1)
    {
      cout << "MASA :: v gradient (3) error condition failed!\n";
      exit(1);
    }

  // p
  derr = masa_eval_2d_grad_p<Scalar>(0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: p gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_2d_grad_p<Scalar>(0,0,3);
  if(derr != -1)
    {
      cout << "MASA :: p gradient (3) error condition failed!\n";
      exit(1);
    }

  // rho
  derr = masa_eval_2d_grad_rho<Scalar>(0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: rho gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_2d_grad_rho<Scalar>(0,0,3);
  if(derr != -1)
    {
      cout << "MASA :: rho gradient (3) error condition failed!\n";
      exit(1);
    }

  // tests passed
  return 0;
}
