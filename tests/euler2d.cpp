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
// euler2d.cpp :program that tests euler2d from masa against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <limits>

using namespace MASA;
using namespace std;

template<typename Scalar>
Scalar nancheck(Scalar x)
{
  if(isnan(x))
    {
      cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
  return 1;
}

template<typename Scalar>
Scalar anQ_p (Scalar x,Scalar y,Scalar p_0,Scalar p_x,Scalar p_y,Scalar a_px,Scalar a_py,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return p_an;
}
  
template<typename Scalar>
Scalar anQ_u (Scalar x,Scalar y,Scalar u_0,Scalar u_x,Scalar u_y,Scalar a_ux,Scalar a_uy,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return u_an;
} 

template<typename Scalar> 
Scalar anQ_v (Scalar x,Scalar y,Scalar v_0,Scalar v_x,Scalar v_y,Scalar a_vx,Scalar a_vy,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return v_an;
}

template<typename Scalar>
Scalar anQ_rho (Scalar x,Scalar y,Scalar rho_0,Scalar rho_x,Scalar rho_y,Scalar a_rhox,Scalar a_rhoy,Scalar L)
{ 
  Scalar pi = acos(-1);
  Scalar rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return rho_an;
}

template<typename Scalar>
Scalar SourceQ_e (
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar Gamma,
  Scalar mu,
  Scalar L)
{
  Scalar pi = acos(-1);
  Scalar Q_e;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar P;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  P = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);

  Q_e = -a_px * pi * p_x * Gamma * U * sin(a_px * pi * x / L) / (Gamma - 0.1e1) / L + a_py * pi * p_y * Gamma * V * cos(a_py * pi * y / L) / (Gamma - 0.1e1) / L + (U * U + V * V) * a_rhox * pi * rho_x * U * cos(a_rhox * pi * x / L) / L / 0.2e1 - (U * U + V * V) * a_rhoy * pi * rho_y * V * sin(a_rhoy * pi * y / L) / L / 0.2e1 + (0.3e1 * a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * U * U / L / 0.2e1 - (a_uy * u_y * sin(a_uy * pi * y / L) + a_vx * v_x * sin(a_vx * pi * x / L)) * pi * RHO * U * V / L + (a_ux * u_x * cos(a_ux * pi * x / L) + 0.3e1 * a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * V * V / L / 0.2e1 + (a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * Gamma * P / (Gamma - 0.1e1) / L;

  return(Q_e);
}

template<typename Scalar>
Scalar SourceQ_u ( // 23 variables
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar L)
{
  Scalar pi = acos(-1);
  Scalar Q_u;
  Scalar RHO;
  Scalar U;
  Scalar V;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);

  Q_u = a_rhox * pi * rho_x * U * U * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * U * V * sin(a_rhoy * pi * y / L) / L - a_uy * pi * u_y * RHO * V * sin(a_uy * pi * y / L) / L - a_px * pi * p_x * sin(a_px * pi * x / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * U / L;

  return(Q_u);
}

template<typename Scalar>
Scalar SourceQ_v (
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar L)
{

  Scalar pi = acos(-1);
  Scalar Q_v;
  Scalar RHO;
  Scalar U;
  Scalar V;
  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);

  Q_v = a_rhox * pi * rho_x * U * V * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * V * sin(a_rhoy * pi * y / L) / L - a_vx * pi * v_x * RHO * U * sin(a_vx * pi * x / L) / L + a_py * pi * p_y * cos(a_py * pi * y / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO * V / L;

  return(Q_v);
}

template<typename Scalar>
Scalar SourceQ_rho(
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar L)
{
  Scalar pi = acos(-1);
  Scalar Q_rho;
  Scalar RHO;
  Scalar U;
  Scalar V;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);

  Q_rho = a_rhox * pi * rho_x * U * cos(a_rhox * pi * x / L) / L - a_rhoy * pi * rho_y * V * sin(a_rhoy * pi * y / L) / L + (a_ux * u_x * cos(a_ux * pi * x / L) + a_vy * v_y * cos(a_vy * pi * y / L)) * pi * RHO / L;

  return(Q_rho);
}

template<typename Scalar>
int run_regression()
{  
  
  // need to add variable based on precision
  Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();
  
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

  // parameters
  Scalar x;
  Scalar y;

  // solutions -- efield is MASA term, efield2 is maple, efield3 is abs error between them
  Scalar ufield,ufield2,ufield3;
  Scalar vfield,vfield2,vfield3;
  Scalar efield,efield2,efield3;
  Scalar rho,rho2,rho3;
  Scalar gradx,grady,gradz,gradp,gradprho;

  Scalar u_an,u_an2,u_an3;
  Scalar v_an,v_an2,v_an3;
  Scalar p_an,p_an2,p_an3;
  Scalar rho_an,rho_an2,rho_an3;

  // initalize
  int nx = 115;  // number of points
  int ny = 68;  
  int lx=3;     // length
  int ly=1; 
  
  Scalar dx=Scalar(lx)/Scalar(nx);
  Scalar dy=Scalar(ly)/Scalar(ny);

  masa_init<Scalar>("euler-test","euler_2d");

  // set params
  masa_init_param<Scalar>();
  
  // get vars
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

  // check that all terms have been initialized
  int err = masa_sanity_check<Scalar>();
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
	ufield2 = SourceQ_u<Scalar> (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);
	vfield2 = SourceQ_v<Scalar>   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);
	rho2    = SourceQ_rho<Scalar> (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);  
	efield2 = SourceQ_e<Scalar>   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,Gamma,mu,L);
	
	u_an2   = anQ_u<Scalar>   (x,y,u_0,u_x,u_y,a_ux,a_uy,L);
	v_an2   = anQ_v<Scalar>   (x,y,v_0,v_x,v_y,a_vx,a_vy,L);
	rho_an2 = anQ_rho<Scalar> (x,y,rho_0,rho_x,rho_y,a_rhox,a_rhoy,L);
	p_an2   = anQ_p<Scalar>   (x,y,p_0,p_x,p_y,a_px,a_py,L);

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
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "U Field Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << ufield3 << endl;
	    cout << "Source term is:                   " << ufield2 << endl;
	    cout << "MASA term is:                     " << ufield << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(u_an3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "U Field Analytical Term\n";
	    cout << "Exceeded Threshold by: " << u_an << endl;
	    cout.precision(16);
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(vfield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "V Field Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << vfield3 << endl;
	    cout << "Source term is:                   " << vfield2 << endl;
	    cout << "MASA term is:                     " << vfield << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(v_an3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "V Field Analytical Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << v_an3 << endl;
	    cout << "Source term is:        " << v_an2 << endl;
	    cout << "MASA term is:          " << v_an << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(efield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "Energy Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << efield3 << endl;
	    cout << "Threshold was: " << threshold << endl;
	    cout << "Source term is:        " << efield2 << endl;
	    cout << "MASA term is:          " << efield << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(p_an3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "P Field Analytical Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << p_an << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(rho3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout.precision(16);
	    cout << "RHO Source Term\n";
	    cout << "Exceeded Threshold by: " << rho << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(rho_an3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout.precision(16);
	    cout << "RHO Analytical Term\n";
	    cout << "Exceeded Threshold by: " << rho_an << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	// adding a new error check: ensure physical results are coming out!
	if(0 > rho)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "Initial Variables are returning non-physical results!\n";
	    cout << "RHO\n";
	    cout << rho << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(0 > rho_an)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "Initial Variables are returning non-physical results!\n";
	    cout << "RHO analytical\n";
	    exit(1);
	  }

	if(0 > p_an)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "Initial Variables are returning non-physical results!\n";
	    cout << "Pressure is negative!\n";
	    exit(1);
	  }

	/*
	if(0 > efield)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Euler-2d\n";
	    cout << "Initial Variables are returning non-physical results!\n";
	    cout << "Energy is negative!\n";
	    cout << efield << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }
	*/
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

  // all tests passed
  return 0;
}

// queue
int main()
{
  int err=0;

  err += run_regression<double>();
  //err += run_regression<long double>();

  return err;
}
