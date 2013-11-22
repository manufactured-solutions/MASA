// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
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
// $Id: ns2d.cpp 15166 2010-10-20 03:26:27Z roystgnr $
//
// ns2d.cpp: program that tests navier-stokes-2d against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <tests.h>
#include <cmath>

using namespace std;
using namespace MASA;

template<typename Scalar>
Scalar anQ_p (Scalar x,Scalar y,Scalar p_0,Scalar p_x,Scalar p_y,Scalar a_px,Scalar a_py,Scalar L)
{
  Scalar pi = std::acos(Scalar(-1));
  Scalar exact_p = p_0 + p_x * std::cos(a_px * pi * x / L) + p_y * std::sin(a_py * pi * y / L);
  return exact_p;
}
  
template<typename Scalar>
Scalar anQ_u (Scalar x,Scalar y,Scalar u_0,Scalar u_x,Scalar u_y,Scalar a_ux,Scalar a_uy,Scalar L)
{
  Scalar pi = std::acos(Scalar(-1));
  Scalar exact_u = u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L);
  return exact_u;
} 
 
template<typename Scalar>
Scalar anQ_v (Scalar x,Scalar y,Scalar v_0,Scalar v_x,Scalar v_y,Scalar a_vx,Scalar a_vy,Scalar L)
{
  Scalar pi = std::acos(Scalar(-1));
  Scalar exact_v = v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L);
  return exact_v;
}

template<typename Scalar>
Scalar anQ_rho (Scalar x,Scalar y,Scalar rho_0,Scalar rho_x,Scalar rho_y,Scalar a_rhox,Scalar a_rhoy,Scalar L)
{ 
  Scalar pi = std::acos(Scalar(-1));
  Scalar exact_rho = rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L);
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
  Scalar pi = std::acos(Scalar(-1));
  Scalar Q_e;
  Q_e = -(v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_0) * (std::pow(u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_0, Scalar(0.2e1)) + std::pow(v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_0, Scalar(0.2e1))) * rho_y * std::sin(a_rhoy * pi * y / L) * a_rhoy * pi / L / Scalar(0.2e1) + (u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_0) * (std::pow(u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_0, Scalar(0.2e1)) + std::pow(v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_0, Scalar(0.2e1))) * rho_x * std::cos(a_rhox * pi * x / L) * a_rhox * pi / L / Scalar(0.2e1) + Scalar(0.4e1) / Scalar(0.3e1) * (v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_0) * mu * v_y * std::sin(a_vy * pi * y / L) * a_vy * a_vy * pi * pi * std::pow(L, Scalar(-0.2e1)) - Scalar(0.4e1) / Scalar(0.3e1) * mu * v_y * v_y * std::pow(cos(a_vy * pi * y / L), Scalar(0.2e1)) * a_vy * a_vy * pi * pi * std::pow(L, Scalar(-0.2e1)) - mu * v_x * v_x * std::pow(sin(a_vx * pi * x / L), Scalar(0.2e1)) * a_vx * a_vx * pi * pi * std::pow(L, Scalar(-0.2e1)) - Scalar(0.4e1) / Scalar(0.3e1) * mu * u_x * u_x * std::pow(cos(a_ux * pi * x / L), Scalar(0.2e1)) * a_ux * a_ux * pi * pi * std::pow(L, Scalar(-0.2e1)) - mu * u_y * u_y * std::pow(sin(a_uy * pi * y / L), Scalar(0.2e1)) * a_uy * a_uy * pi * pi * std::pow(L, Scalar(-0.2e1)) + (Gamma * (p_x * std::cos(a_px * pi * x / L) + p_y * std::sin(a_py * pi * y / L) + p_0) / (Gamma - Scalar(0.1e1)) + (std::pow(u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_0, Scalar(0.2e1)) / Scalar(0.2e1) + Scalar(0.3e1) / Scalar(0.2e1) * std::pow(v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_0, Scalar(0.2e1))) * (rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0)) * v_y * std::cos(a_vy * pi * y / L) * a_vy * pi / L + (Gamma * (p_x * std::cos(a_px * pi * x / L) + p_y * std::sin(a_py * pi * y / L) + p_0) / (Gamma - Scalar(0.1e1)) + (Scalar(0.3e1) / Scalar(0.2e1) * std::pow(u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_0, Scalar(0.2e1)) + std::pow(v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_0, Scalar(0.2e1)) / Scalar(0.2e1)) * (rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0)) * u_x * std::cos(a_ux * pi * x / L) * a_ux * pi / L + (v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_0) * mu * v_x * std::cos(a_vx * pi * x / L) * a_vx * a_vx * pi * pi * std::pow(L, Scalar(-0.2e1)) + Scalar(0.4e1) / Scalar(0.3e1) * (u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_0) * mu * u_x * std::sin(a_ux * pi * x / L) * a_ux * a_ux * pi * pi * std::pow(L, Scalar(-0.2e1)) + (u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_0) * mu * u_y * std::cos(a_uy * pi * y / L) * a_uy * a_uy * pi * pi * std::pow(L, Scalar(-0.2e1)) - (v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_0) * (rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0) * (u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_0) * u_y * std::sin(a_uy * pi * y / L) * a_uy * pi / L - (p_x * std::cos(a_px * pi * x / L) + p_y * std::sin(a_py * pi * y / L) + p_0) * rho_x * k * std::sin(a_rhox * pi * x / L) * a_rhox * a_rhox * pi * pi * std::pow(rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0, Scalar(-0.2e1)) * std::pow(L, Scalar(-0.2e1)) / R - (Scalar(0.2e1) * p_x * std::cos(a_px * pi * x / L) + Scalar(0.2e1) * p_y * std::sin(a_py * pi * y / L) + Scalar(0.2e1) * p_0) * rho_x * rho_x * k * std::pow(cos(a_rhox * pi * x / L), Scalar(0.2e1)) * a_rhox * a_rhox * pi * pi * std::pow(rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0, Scalar(-0.3e1)) * std::pow(L, Scalar(-0.2e1)) / R - (p_x * std::cos(a_px * pi * x / L) + p_y * std::sin(a_py * pi * y / L) + p_0) * rho_y * k * std::cos(a_rhoy * pi * y / L) * a_rhoy * a_rhoy * pi * pi * std::pow(rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0, Scalar(-0.2e1)) * std::pow(L, Scalar(-0.2e1)) / R - (Scalar(0.2e1) * p_x * std::cos(a_px * pi * x / L) + Scalar(0.2e1) * p_y * std::sin(a_py * pi * y / L) + Scalar(0.2e1) * p_0) * rho_y * rho_y * k * std::pow(sin(a_rhoy * pi * y / L), Scalar(0.2e1)) * a_rhoy * a_rhoy * pi * pi * std::pow(rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0, Scalar(-0.3e1)) * std::pow(L, Scalar(-0.2e1)) / R + Scalar(0.4e1) / Scalar(0.3e1) * mu * u_x * v_y * std::cos(a_ux * pi * x / L) * std::cos(a_vy * pi * y / L) * a_ux * a_vy * pi * pi * std::pow(L, Scalar(-0.2e1)) - Scalar(0.2e1) * mu * u_y * v_x * std::sin(a_uy * pi * y / L) * std::sin(a_vx * pi * x / L) * a_uy * a_vx * pi * pi * std::pow(L, Scalar(-0.2e1)) - Scalar(0.2e1) * k * p_x * rho_x * std::cos(a_rhox * pi * x / L) * std::sin(a_px * pi * x / L) * a_px * a_rhox * pi * pi * std::pow(rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0, Scalar(-0.2e1)) * std::pow(L, Scalar(-0.2e1)) / R - Scalar(0.2e1) * k * p_y * rho_y * std::cos(a_py * pi * y / L) * std::sin(a_rhoy * pi * y / L) * a_py * a_rhoy * pi * pi * std::pow(rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0, Scalar(-0.2e1)) * std::pow(L, Scalar(-0.2e1)) / R - (v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_0) * (rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0) * (u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_0) * v_x * std::sin(a_vx * pi * x / L) * a_vx * pi / L - Gamma * (u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_0) * p_x * std::sin(a_px * pi * x / L) * a_px * pi / (Gamma - Scalar(0.1e1)) / L + Gamma * (v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_0) * p_y * std::cos(a_py * pi * y / L) * a_py * pi / (Gamma - Scalar(0.1e1)) / L + k * p_x * std::cos(a_px * pi * x / L) * a_px * a_px * pi * pi / (rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0) * std::pow(L, Scalar(-0.2e1)) / R + k * p_y * std::sin(a_py * pi * y / L) * a_py * a_py * pi * pi / (rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0) * std::pow(L, Scalar(-0.2e1)) / R;
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
  Scalar pi = std::acos(Scalar(-1));
  Scalar Q_u;
  Q_u = Scalar(0.4e1) / Scalar(0.3e1) * mu * u_x * std::sin(a_ux * pi * x / L) * a_ux * a_ux * pi * pi * std::pow(L, Scalar(-0.2e1)) + mu * u_y * std::cos(a_uy * pi * y / L) * a_uy * a_uy * pi * pi * std::pow(L, Scalar(-0.2e1)) - p_x * std::sin(a_px * pi * x / L) * a_px * pi / L + rho_x * std::cos(a_rhox * pi * x / L) * std::pow(u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L), Scalar(0.2e1)) * a_rhox * pi / L - rho_y * std::sin(a_rhoy * pi * y / L) * (v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L)) * (u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L)) * a_rhoy * pi / L + Scalar(0.2e1) * u_x * std::cos(a_ux * pi * x / L) * (rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L)) * (u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L)) * a_ux * pi / L - u_y * std::sin(a_uy * pi * y / L) * (rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L)) * (v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L)) * a_uy * pi / L + v_y * std::cos(a_vy * pi * y / L) * (rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L)) * (u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L)) * a_vy * pi / L;
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
  Scalar pi = std::acos(Scalar(-1));
  Scalar Q_v;
  Q_v = mu * v_x * std::cos(a_vx * pi * x / L) * a_vx * a_vx * pi * pi * std::pow(L, Scalar(-0.2e1)) + Scalar(0.4e1) / Scalar(0.3e1) * mu * v_y * std::sin(a_vy * pi * y / L) * a_vy * a_vy * pi * pi * std::pow(L, Scalar(-0.2e1)) + p_y * std::cos(a_py * pi * y / L) * a_py * pi / L + rho_x * std::cos(a_rhox * pi * x / L) * (v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L)) * (u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L)) * a_rhox * pi / L - rho_y * std::sin(a_rhoy * pi * y / L) * std::pow(v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L), Scalar(0.2e1)) * a_rhoy * pi / L + u_x * std::cos(a_ux * pi * x / L) * (rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L)) * (v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L)) * a_ux * pi / L - v_x * std::sin(a_vx * pi * x / L) * (rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L)) * (u_0 + u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L)) * a_vx * pi / L + Scalar(0.2e1) * v_y * std::cos(a_vy * pi * y / L) * (rho_0 + rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L)) * (v_0 + v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L)) * a_vy * pi / L;
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
  Scalar pi = std::acos(Scalar(-1));
  Scalar Q_rho;
  Q_rho = (u_x * std::sin(a_ux * pi * x / L) + u_y * std::cos(a_uy * pi * y / L) + u_0) * a_rhox * pi * rho_x * std::cos(a_rhox * pi * x / L) / L - (v_x * std::cos(a_vx * pi * x / L) + v_y * std::sin(a_vy * pi * y / L) + v_0) * a_rhoy * pi * rho_y * std::sin(a_rhoy * pi * y / L) / L + (rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0) * a_ux * pi * u_x * std::cos(a_ux * pi * x / L) / L + (rho_x * std::sin(a_rhox * pi * x / L) + rho_y * std::cos(a_rhoy * pi * y / L) + rho_0) * a_vy * pi * v_y * std::cos(a_vy * pi * y / L) / L;
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

  // solutions
  Scalar exact_u,exact_u2,exact_u3;
  Scalar exact_v,exact_v2,exact_v3;
  Scalar exact_p,exact_p2,exact_p3;
  Scalar exact_rho,exact_rho2,exact_rho3;
  Scalar gradx,grady;
  //Scalar gradp,gradrho;
  
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
	ufield = masa_eval_source_rho_u  <Scalar>(x,y);
	vfield = masa_eval_source_rho_v  <Scalar>(x,y);
	efield = masa_eval_source_rho_e  <Scalar>(x,y);
	rho    = masa_eval_source_rho<Scalar>(x,y);
	
	//evaluate analytical terms
	exact_u = masa_eval_exact_u        <Scalar>(x,y);
	exact_v = masa_eval_exact_v        <Scalar>(x,y);
	exact_p = masa_eval_exact_p        <Scalar>(x,y);
	exact_rho = masa_eval_exact_rho    <Scalar>(x,y);
	
	// eval gradient terms
	gradx = masa_eval_grad_u<Scalar>(x,y,1);
	nancheck(gradx);
	grady = masa_eval_grad_u<Scalar>(x,y,2);		
	nancheck(grady);

	gradx = masa_eval_grad_v<Scalar>(x,y,1);
	nancheck(gradx);
	grady = masa_eval_grad_v<Scalar>(x,y,2);		
	nancheck(grady);

	gradx = masa_eval_grad_p<Scalar>(x,y,1);
	nancheck(gradx);
	grady = masa_eval_grad_p<Scalar>(x,y,2);		
	nancheck(grady);
  
	gradx = masa_eval_grad_rho<Scalar>(x,y,1);
	nancheck(gradx);
	grady = masa_eval_grad_rho<Scalar>(x,y,2);		
	nancheck(grady);

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

	ufield3 = std::abs(ufield-ufield2);
	vfield3 = std::abs(vfield-vfield2);
	efield3 = std::abs(efield-efield2);
	rho3    = std::abs(rho-rho2);

	exact_u3   = std::abs(exact_u-exact_u2);
	exact_v3   = std::abs(exact_v-exact_v2);
	exact_rho3 = std::abs(exact_rho-exact_rho2);
	exact_p3   = std::abs(exact_p-exact_p2);

#else

	ufield3 = std::abs(ufield-ufield2)/std::abs(ufield2);
	vfield3 = std::abs(vfield-vfield2)/std::abs(vfield2);
	efield3 = std::abs(efield-efield2)/std::abs(efield2);
	rho3    = std::abs(rho-rho2)/std::abs(rho2);

	exact_u3   = std::abs(exact_u-exact_u2)/std::abs(exact_u2);
	exact_v3   = std::abs(exact_v-exact_v2)/std::abs(exact_v2);
	exact_rho3 = std::abs(exact_rho-exact_rho2)/std::abs(exact_rho2);
	exact_p3   = std::abs(exact_p-exact_p2)/std::abs(exact_p2);

#endif	

	threshcheck(ufield3,threshold);
	threshcheck(vfield3,threshold);
	threshcheck(efield3,threshold);
	threshcheck(rho3,threshold);
	
	threshcheck(exact_u3,threshold);
	threshcheck(exact_v3,threshold);
	threshcheck(exact_rho3,threshold);
	threshcheck(exact_p3,threshold);

      } // done iterating

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.
  freopen("/dev/null","w",stdout);  
  
  // test gradient error terms
  Scalar derr = masa_eval_grad_u<Scalar>(0,0,0);
  test_grad(derr);
  
  derr = masa_eval_grad_u<Scalar>(0,0,3);
  test_grad(derr);

  // v
  derr = masa_eval_grad_v<Scalar>(0,0,0);
  test_grad(derr);
  
  derr = masa_eval_grad_v<Scalar>(0,0,3);
  test_grad(derr);

  // p
  derr = masa_eval_grad_p<Scalar>(0,0,0);
  test_grad(derr);
  
  derr = masa_eval_grad_p<Scalar>(0,0,3);
  test_grad(derr);

  // rho
  derr = masa_eval_grad_rho<Scalar>(0,0,0);
  test_grad(derr);
  
  derr = masa_eval_grad_rho<Scalar>(0,0,3);
  test_grad(derr);

  // tests passed
  return 0;
}

int main()
{
  int err=0;

  err += run_regression<double>();
  err += run_regression<long double>();

  return err;
}
