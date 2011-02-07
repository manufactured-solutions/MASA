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
// $Author: roystgnr $
// $Id: ns2d.cpp 15166 2010-10-20 03:26:27Z roystgnr $
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
#include <limits>

using namespace std;
using namespace MASA;

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
  Scalar exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return exact_p;
}
  
template<typename Scalar>
Scalar anQ_u (Scalar x,Scalar y,Scalar u_0,Scalar u_x,Scalar u_y,Scalar a_ux,Scalar a_uy,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return exact_u;
} 
 
template<typename Scalar>
Scalar anQ_v (Scalar x,Scalar y,Scalar v_0,Scalar v_x,Scalar v_y,Scalar a_vx,Scalar a_vy,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return exact_v;
}

template<typename Scalar>
Scalar anQ_rho (Scalar x,Scalar y,Scalar rho_0,Scalar rho_x,Scalar rho_y,Scalar a_rhox,Scalar a_rhoy,Scalar L)
{ 
  Scalar pi = acos(-1);
  Scalar exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return exact_rho;
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
  Scalar L,
  Scalar R,
  Scalar k)
{
  Scalar pi = acos(-1);
  Scalar Q_e;
  Q_e = -(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, 0.2e1)) * rho_y * sin(a_rhoy * pi * y / L) * a_rhoy * pi / L / 0.2e1 + (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, 0.2e1)) * rho_x * cos(a_rhox * pi * x / L) * a_rhox * pi / L / 0.2e1 + 0.4e1 / 0.3e1 * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * mu * v_y * sin(a_vy * pi * y / L) * a_vy * a_vy * pi * pi * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * v_y * v_y * pow(cos(a_vy * pi * y / L), 0.2e1) * a_vy * a_vy * pi * pi * pow(L, -0.2e1) - mu * v_x * v_x * pow(sin(a_vx * pi * x / L), 0.2e1) * a_vx * a_vx * pi * pi * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * u_x * u_x * pow(cos(a_ux * pi * x / L), 0.2e1) * a_ux * a_ux * pi * pi * pow(L, -0.2e1) - mu * u_y * u_y * pow(sin(a_uy * pi * y / L), 0.2e1) * a_uy * a_uy * pi * pi * pow(L, -0.2e1) + (Gamma * (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) / (Gamma - 0.1e1) + (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1) / 0.2e1 + 0.3e1 / 0.2e1 * pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, 0.2e1)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0)) * v_y * cos(a_vy * pi * y / L) * a_vy * pi / L + (Gamma * (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) / (Gamma - 0.1e1) + (0.3e1 / 0.2e1 * pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, 0.2e1) / 0.2e1) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0)) * u_x * cos(a_ux * pi * x / L) * a_ux * pi / L + (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * mu * v_x * cos(a_vx * pi * x / L) * a_vx * a_vx * pi * pi * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * mu * u_x * sin(a_ux * pi * x / L) * a_ux * a_ux * pi * pi * pow(L, -0.2e1) + (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * mu * u_y * cos(a_uy * pi * y / L) * a_uy * a_uy * pi * pi * pow(L, -0.2e1) - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * u_y * sin(a_uy * pi * y / L) * a_uy * pi / L - (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) * rho_x * k * sin(a_rhox * pi * x / L) * a_rhox * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (0.2e1 * p_x * cos(a_px * pi * x / L) + 0.2e1 * p_y * sin(a_py * pi * y / L) + 0.2e1 * p_0) * rho_x * rho_x * k * pow(cos(a_rhox * pi * x / L), 0.2e1) * a_rhox * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.3e1) * pow(L, -0.2e1) / R - (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) * rho_y * k * cos(a_rhoy * pi * y / L) * a_rhoy * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (0.2e1 * p_x * cos(a_px * pi * x / L) + 0.2e1 * p_y * sin(a_py * pi * y / L) + 0.2e1 * p_0) * rho_y * rho_y * k * pow(sin(a_rhoy * pi * y / L), 0.2e1) * a_rhoy * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.3e1) * pow(L, -0.2e1) / R + 0.4e1 / 0.3e1 * mu * u_x * v_y * cos(a_ux * pi * x / L) * cos(a_vy * pi * y / L) * a_ux * a_vy * pi * pi * pow(L, -0.2e1) - 0.2e1 * mu * u_y * v_x * sin(a_uy * pi * y / L) * sin(a_vx * pi * x / L) * a_uy * a_vx * pi * pi * pow(L, -0.2e1) - 0.2e1 * k * p_x * rho_x * cos(a_rhox * pi * x / L) * sin(a_px * pi * x / L) * a_px * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - 0.2e1 * k * p_y * rho_y * cos(a_py * pi * y / L) * sin(a_rhoy * pi * y / L) * a_py * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * v_x * sin(a_vx * pi * x / L) * a_vx * pi / L - Gamma * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * p_x * sin(a_px * pi * x / L) * a_px * pi / (Gamma - 0.1e1) / L + Gamma * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * p_y * cos(a_py * pi * y / L) * a_py * pi / (Gamma - 0.1e1) / L + k * p_x * cos(a_px * pi * x / L) * a_px * a_px * pi * pi / (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * pow(L, -0.2e1) / R + k * p_y * sin(a_py * pi * y / L) * a_py * a_py * pi * pi / (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * pow(L, -0.2e1) / R;
  return(Q_e);
}

template<typename Scalar>
Scalar SourceQ_u (
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
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k)
{
  Scalar pi = acos(-1);
  Scalar Q_u;
  Q_u = 0.4e1 / 0.3e1 * mu * u_x * sin(a_ux * pi * x / L) * a_ux * a_ux * pi * pi * pow(L, -0.2e1) + mu * u_y * cos(a_uy * pi * y / L) * a_uy * a_uy * pi * pi * pow(L, -0.2e1) - p_x * sin(a_px * pi * x / L) * a_px * pi / L + rho_x * cos(a_rhox * pi * x / L) * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L), 0.2e1) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_rhoy * pi / L + 0.2e1 * u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_ux * pi / L - u_y * sin(a_uy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_uy * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_vy * pi / L;
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
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k)
{
  Scalar pi = acos(-1);
  Scalar Q_v;
  Q_v = mu * v_x * cos(a_vx * pi * x / L) * a_vx * a_vx * pi * pi * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * mu * v_y * sin(a_vy * pi * y / L) * a_vy * a_vy * pi * pi * pow(L, -0.2e1) + p_y * cos(a_py * pi * y / L) * a_py * pi / L + rho_x * cos(a_rhox * pi * x / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L), 0.2e1) * a_rhoy * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_ux * pi / L - v_x * sin(a_vx * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_vx * pi / L + 0.2e1 * v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_vy * pi / L;
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
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k)

{
  Scalar pi = acos(-1);
  Scalar Q_rho;
  Q_rho = (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * a_rhox * pi * rho_x * cos(a_rhox * pi * x / L) / L - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * a_rhoy * pi * rho_y * sin(a_rhoy * pi * y / L) / L + (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * a_ux * pi * u_x * cos(a_ux * pi * x / L) / L + (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * a_vy * pi * v_y * cos(a_vy * pi * y / L) / L;
  return(Q_rho);
}

template<typename Scalar>
int run_regression()
{

  const Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();

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
  Scalar exact_u,exact_u2,exact_u3;
  Scalar exact_v,exact_v2,exact_v3;
  Scalar exact_w,exact_w2,exact_w3;
  Scalar exact_p,exact_p2,exact_p3;
  Scalar exact_rho,exact_rho2,exact_rho3;
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
	ufield = masa_eval_source_u  <Scalar>(x,y);
	vfield = masa_eval_source_v  <Scalar>(x,y);
	efield = masa_eval_source_e  <Scalar>(x,y);
	rho    = masa_eval_source_rho<Scalar>(x,y);
	
	//evaluate analytical terms
	exact_u = masa_eval_exact_u        <Scalar>(x,y);
	exact_v = masa_eval_exact_v        <Scalar>(x,y);
	exact_p = masa_eval_exact_p        <Scalar>(x,y);
	exact_rho = masa_eval_exact_rho    <Scalar>(x,y);
	
	// eval gradient terms
	gradx = masa_eval_grad_2d_u<Scalar>(x,y,1);
	grady = masa_eval_grad_2d_u<Scalar>(x,y,2);		

	gradx = masa_eval_grad_2d_v<Scalar>(x,y,1);
	grady = masa_eval_grad_2d_v<Scalar>(x,y,2);		

	gradx = masa_eval_grad_2d_p<Scalar>(x,y,1);
	grady = masa_eval_grad_2d_p<Scalar>(x,y,2);		
  
	gradx = masa_eval_grad_2d_rho<Scalar>(x,y,1);
	grady = masa_eval_grad_2d_rho<Scalar>(x,y,2);		

	// check against maple
	ufield2 = SourceQ_u   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);
	vfield2 = SourceQ_v   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);
	rho2    = SourceQ_rho (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);  
	efield2 = SourceQ_e   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,Gamma,mu,L,R,k);
	
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
	
	nancheck(ufield3);
	nancheck(vfield3);
	nancheck(efield3);
	nancheck(rho3);
	
	nancheck(exact_u3);
	nancheck(exact_v3);
	nancheck(exact_rho3);
	nancheck(exact_p3);


	if(ufield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n";
	    cout << "U Field Source Term\n";
	    cout << "Exceeded Threshold by: " << ufield << endl;
	    exit(1);
	  }

	if(exact_u3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n";
	    cout << "U Field Analytical Term\n";
	    cout << "Exceeded Threshold by: " << exact_u3 << endl;
	    exit(1);
	  }

	if(vfield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n";
	    cout << "V Field Source Term\n";
	    cout << "Exceeded Threshold of: " << threshold <<  " by: " << vfield3 << endl;
	    exit(1);
	  }

	if(exact_v3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n";
	    cout << "V Field Analytical Term\n";
	    cout << "Exceeded Threshold of: " << threshold <<  " by: " << exact_v3 << endl;
	    exit(1);
	  }

	if(efield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n";
	    cout << "Energy Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << efield3 << endl;
	    cout << "Source term is:                   " << efield2 << endl;
	    cout << "MASA term is:                     " << efield << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(exact_p3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n";
	    cout << "P Field Analytical Term\n";
	    cout << "Exceeded Threshold by: " << exact_p3 << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(rho3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n";
	    cout << "RHO Source Term\n";
	    cout << "Exceeded Threshold by: " << rho3 << endl;
	    cout << x << " " << y <<  endl;
	    exit(1);
	  }

	if(exact_rho3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 2d\n";
	    cout << "RHO Analytical Term\n";
	    cout << "Exceeded Threshold by: " << exact_rho3 << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

      } // done iterating

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.
  freopen("/dev/null","w",stdout);  
  
  // test gradient error terms
  Scalar derr = masa_eval_grad_2d_u<Scalar>(0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_grad_2d_u<Scalar>(0,0,3);
  if(derr != -1)
    {
      cout << "MASA :: gradient (3) error condition failed!\n";
      exit(1);
    }

  // v
  derr = masa_eval_grad_2d_v<Scalar>(0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: v gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_grad_2d_v<Scalar>(0,0,3);
  if(derr != -1)
    {
      cout << "MASA :: v gradient (3) error condition failed!\n";
      exit(1);
    }

  // p
  derr = masa_eval_grad_2d_p<Scalar>(0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: p gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_grad_2d_p<Scalar>(0,0,3);
  if(derr != -1)
    {
      cout << "MASA :: p gradient (3) error condition failed!\n";
      exit(1);
    }

  // rho
  derr = masa_eval_grad_2d_rho<Scalar>(0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: rho gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_grad_2d_rho<Scalar>(0,0,3);
  if(derr != -1)
    {
      cout << "MASA :: rho gradient (3) error condition failed!\n";
      exit(1);
    }

  // tests passed
  return 0;
}

int main()
{
  int err=0;

  err += run_regression<double>();
  //err += run_regression<long double>();

  return err;
}
