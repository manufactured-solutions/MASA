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
// cns.cpp: These are the MASA class member functions and constructors
//          For the Compressible Navier-Stokes
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         Compressible Navier Stokes Equations
 *
 *         2D
 *
 * -----------------------------------------------
 */ 

MASA::navierstokes_2d_compressible::navierstokes_2d_compressible()
{
  mmsname = "navierstokes_2d_compressible";
  dimension=2;

  register_var("R",&R);
  register_var("k",&k);

  register_var("u_0",&u_0);
  register_var("u_x",&u_x);
  register_var("u_y",&u_y);
  register_var("v_0",&v_0);
  register_var("v_x",&v_x);
  register_var("v_y",&v_y);
  register_var("rho_0",&rho_0);
  register_var("rho_x",&rho_x);
  register_var("rho_y",&rho_y);
  register_var("p_0",&p_0);
  register_var("p_x",&p_x);
  register_var("p_y",&p_y);
  register_var("a_px",&a_px);
  register_var("a_py",&a_py);
  register_var("a_rhox",&a_rhox);
  register_var("a_rhoy",&a_rhoy);
  register_var("a_ux",&a_ux);
  register_var("a_uy",&a_uy);
  register_var("a_vx",&a_vx);
  register_var("a_vy",&a_vy);
  register_var("L",&L);
  register_var("Gamma",&Gamma);
  register_var("mu",&mu);

}//done with constructor

int MASA::navierstokes_2d_compressible::init_var()
{
  int err = 0;

  // currently randomly generated
  err += set_var("R",1.01);
  err += set_var("k",1.38);

  err += set_var("u_0",1.23);
  err += set_var("u_x",1.1);
  err += set_var("u_y",.08);
  err += set_var("v_0",12);
  err += set_var("v_x",1.6);
  err += set_var("v_y",.67);
  err += set_var("rho_0",1.02);
  err += set_var("rho_x",7.2);
  err += set_var("rho_y",9.8);
  err += set_var("p_0",1.2);
  err += set_var("p_x",.91);
  err += set_var("p_y",.623);
  err += set_var("a_px",.165);
  err += set_var("a_py",.612);
  err += set_var("a_rhox",.627);
  err += set_var("a_rhoy",.828);
  err += set_var("a_ux",.1987);
  err += set_var("a_uy",1.189);
  err += set_var("a_vx",1.91);
  err += set_var("a_vy",2.901);
  err += set_var("Gamma",1.01);
  err += set_var("mu",.918);
  err += set_var("L",3.02);

  return err;

} // done with variable initializer

// ----------------------------------------
//   Source Terms
// ----------------------------------------

double MASA::navierstokes_2d_compressible::eval_q_u(double x,double y)
{
  double Q_u;
  Q_u = 0.4e1 / 0.3e1 * mu * u_x * sin(a_ux * PI * x / L) * a_ux * a_ux * PI * PI * pow(L, -0.2e1) + mu * u_y * cos(a_uy * PI * y / L) * a_uy * a_uy * PI * PI * pow(L, -0.2e1) - p_x * sin(a_px * PI * x / L) * a_px * PI / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L), 0.2e1) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhoy * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_ux * PI / L - u_y * sin(a_uy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_uy * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vy * PI / L;
  return(Q_u);
}

double MASA::navierstokes_2d_compressible::eval_q_v(double x,double y)
{
  double Q_v;
  Q_v = mu * v_x * cos(a_vx * PI * x / L) * a_vx * a_vx * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * mu * v_y * sin(a_vy * PI * y / L) * a_vy * a_vy * PI * PI * pow(L, -0.2e1) + p_y * cos(a_py * PI * y / L) * a_py * PI / L + rho_x * cos(a_rhox * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L), 0.2e1) * a_rhoy * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_ux * PI / L - v_x * sin(a_vx * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vx * PI / L + 0.2e1 * v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_vy * PI / L;
  return(Q_v);
}

double MASA::navierstokes_2d_compressible::eval_q_rho(double x,double y)
{
  double Q_rho;
  Q_rho = (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * a_rhox * PI * rho_x * cos(a_rhox * PI * x / L) / L - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * a_rhoy * PI * rho_y * sin(a_rhoy * PI * y / L) / L + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_ux * PI * u_x * cos(a_ux * PI * x / L) / L + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_vy * PI * v_y * cos(a_vy * PI * y / L) / L;
  return(Q_rho);
}

double MASA::navierstokes_2d_compressible::eval_q_e(double x,double y)
{
  double Q_e;
  Q_e = -(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * (pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1)) * rho_y * sin(a_rhoy * PI * y / L) * a_rhoy * PI / L / 0.2e1 + (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * (pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1)) * rho_x * cos(a_rhox * PI * x / L) * a_rhox * PI / L / 0.2e1 + 0.4e1 / 0.3e1 * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * mu * v_y * sin(a_vy * PI * y / L) * a_vy * a_vy * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * v_y * v_y * pow(cos(a_vy * PI * y / L), 0.2e1) * a_vy * a_vy * PI * PI * pow(L, -0.2e1) - mu * v_x * v_x * pow(sin(a_vx * PI * x / L), 0.2e1) * a_vx * a_vx * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * u_x * u_x * pow(cos(a_ux * PI * x / L), 0.2e1) * a_ux * a_ux * PI * PI * pow(L, -0.2e1) - mu * u_y * u_y * pow(sin(a_uy * PI * y / L), 0.2e1) * a_uy * a_uy * PI * PI * pow(L, -0.2e1) + (Gamma * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) / (Gamma - 0.1e1) + (pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) / 0.2e1 + 0.3e1 / 0.2e1 * pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0)) * v_y * cos(a_vy * PI * y / L) * a_vy * PI / L + (Gamma * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) / (Gamma - 0.1e1) + (0.3e1 / 0.2e1 * pow(u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0, 0.2e1) / 0.2e1) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0)) * u_x * cos(a_ux * PI * x / L) * a_ux * PI / L + (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * mu * v_x * cos(a_vx * PI * x / L) * a_vx * a_vx * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * mu * u_x * sin(a_ux * PI * x / L) * a_ux * a_ux * PI * PI * pow(L, -0.2e1) + (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * mu * u_y * cos(a_uy * PI * y / L) * a_uy * a_uy * PI * PI * pow(L, -0.2e1) - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * u_y * sin(a_uy * PI * y / L) * a_uy * PI / L - (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * rho_x * k * sin(a_rhox * PI * x / L) * a_rhox * a_rhox * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (0.2e1 * p_x * cos(a_px * PI * x / L) + 0.2e1 * p_y * sin(a_py * PI * y / L) + 0.2e1 * p_0) * rho_x * rho_x * k * pow(cos(a_rhox * PI * x / L), 0.2e1) * a_rhox * a_rhox * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.3e1) * pow(L, -0.2e1) / R - (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * rho_y * k * cos(a_rhoy * PI * y / L) * a_rhoy * a_rhoy * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (0.2e1 * p_x * cos(a_px * PI * x / L) + 0.2e1 * p_y * sin(a_py * PI * y / L) + 0.2e1 * p_0) * rho_y * rho_y * k * pow(sin(a_rhoy * PI * y / L), 0.2e1) * a_rhoy * a_rhoy * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.3e1) * pow(L, -0.2e1) / R + 0.4e1 / 0.3e1 * mu * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) * a_ux * a_vy * PI * PI * pow(L, -0.2e1) - 0.2e1 * mu * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) * a_uy * a_vx * PI * PI * pow(L, -0.2e1) - 0.2e1 * k * p_x * rho_x * cos(a_rhox * PI * x / L) * sin(a_px * PI * x / L) * a_px * a_rhox * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - 0.2e1 * k * p_y * rho_y * cos(a_py * PI * y / L) * sin(a_rhoy * PI * y / L) * a_py * a_rhoy * PI * PI * pow(rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0, -0.2e1) * pow(L, -0.2e1) / R - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * v_x * sin(a_vx * PI * x / L) * a_vx * PI / L - Gamma * (u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_0) * p_x * sin(a_px * PI * x / L) * a_px * PI / (Gamma - 0.1e1) / L + Gamma * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_0) * p_y * cos(a_py * PI * y / L) * a_py * PI / (Gamma - 0.1e1) / L + k * p_x * cos(a_px * PI * x / L) * a_px * a_px * PI * PI / (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * pow(L, -0.2e1) / R + k * p_y * sin(a_py * PI * y / L) * a_py * a_py * PI * PI / (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * pow(L, -0.2e1) / R;
  return(Q_e);
}

// ----------------------------------------
//   Gradient of Source Terms
// ----------------------------------------

double MASA::navierstokes_2d_compressible::eval_2d_g_u(double x,double y, int i)
{

  double grad = -1;

  switch(i)
    {
    case 0:
      cout << "MASA error:: eval_2d_g_u has no 0th component\n";
      cout << "Try 1 or 2\n";
      break;
      
    case 1:
      grad =  u_x * cos(a_ux * pi * x / L) * a_ux * pi / L;      
      break;

    case 2:
      grad = -u_y * sin(a_uy * pi * y / L) * a_uy * pi / L;
      break;

    default:
      cout << "MASA error:: eval_2d_g_u has no " << i << "th component\n";
      cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

double MASA::navierstokes_2d_compressible::eval_2d_g_v(double x,double y, int i)
{

  double grad = -1;

  switch(i)
    {
    case 0:
      cout << "MASA error:: eval_2d_g_v has no 0th component\n";
      cout << "Try 1 or 2\n";
      break;
      
    case 1:
      grad = -v_x * sin(a_vx * pi * x / L) * a_vx * pi / L;
      break;

    case 2:
      grad = v_y * cos(a_vy * pi * y / L) * a_vy * pi / L;
      break;

    default:
      cout << "MASA error:: eval_2d_g_v has no " << i << "th component\n";
      cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

double MASA::navierstokes_2d_compressible::eval_2d_g_p(double x,double y, int i)
{

  double grad = -1;

  switch(i)
    {
    case 0:
      cout << "MASA error:: eval_2d_g_p has no 0th component\n";
      cout << "Try 1 or 2\n";
      break;
      
    case 1:
      grad = -p_x * sin(a_px * pi * x / L) * a_px * pi / L;
      break;

    case 2:
      grad =  p_y * cos(a_py * pi * y / L) * a_py * pi / L;
      break;

    default:
      cout << "MASA error:: eval_2d_g_p has no " << i << "th component\n";
      cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

double MASA::navierstokes_2d_compressible::eval_2d_g_rho(double x,double y, int i)
{

  double grad = -1;

  switch(i)
    {
    case 0:
      cout << "MASA error:: eval_2d_g_rho has no 0th component\n";
      cout << "Try 1 or 2\n";
      break;
      
    case 1:
      grad =  rho_x * cos(a_rhox * pi * x / L) * a_rhox * pi / L;
      break;

    case 2:
      grad = -rho_y * sin(a_rhoy * pi * y / L) * a_rhoy * pi / L;
      break;

    default:
      cout << "MASA error:: eval_2d_g_rho has no " << i << "th component\n";
      cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

double MASA::navierstokes_2d_compressible::eval_an_u(double x,double y)
{
  double u_an;
  u_an = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  return u_an;
}

double MASA::navierstokes_2d_compressible::eval_an_v(double x,double y)
{
  double v_an;
  v_an = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);
  return v_an;
}

double MASA::navierstokes_2d_compressible::eval_an_p(double x,double y)
{
  double p_an;
  p_an = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L);
  return p_an;
}

double MASA::navierstokes_2d_compressible::eval_an_rho(double x,double y)
{
  double rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L);
  return rho_an;
}

/* ------------------------------------------------
 *
 *         Compressible Navier Stokes Equations
 *
 *         3D
 *
 * -----------------------------------------------
 */ 

MASA::navierstokes_3d_compressible::navierstokes_3d_compressible()
{
  mmsname = "navierstokes_3d_compressible";
  dimension=3;

  register_var("R",&R);
  register_var("k",&k);

  register_var("u_0",&u_0);
  register_var("u_x",&u_x);
  register_var("u_y",&u_y);
  register_var("u_z",&u_z);
  register_var("v_0",&v_0);
  register_var("v_x",&v_x);
  register_var("v_y",&v_y);
  register_var("v_z",&v_z);
  register_var("rho_0",&rho_0);
  register_var("rho_x",&rho_x);
  register_var("rho_y",&rho_y);
  register_var("rho_z",&rho_z);
  register_var("p_0",&p_0);
  register_var("p_x",&p_x);
  register_var("p_y",&p_y);
  register_var("p_z",&p_z);
  register_var("a_px",&a_px);
  register_var("a_py",&a_py);
  register_var("a_rhox",&a_rhox);
  register_var("a_rhoy",&a_rhoy);
  register_var("a_ux",&a_ux);
  register_var("a_uy",&a_uy);
  register_var("a_uz",&a_uz);
  register_var("a_vx",&a_vx);
  register_var("a_vy",&a_vy);
  register_var("L",&L);
  register_var("Gamma",&Gamma);
  register_var("mu",&mu);
  register_var("a_pz",&a_pz);
  register_var("a_rhoz",&a_rhoz);
  register_var("a_vz",&a_vz);
  register_var("a_wz",&a_wz);
  register_var("a_wx",&a_wx);
  register_var("a_wy",&a_wy);
  register_var("w_0",&w_0);
  register_var("w_x",&w_x);
  register_var("w_y",&w_y);
  register_var("w_z",&w_z);

}//done with constructor

int MASA::navierstokes_3d_compressible::init_var()
{
  int err = 0;

  // set using values from table A.5 
  // in a paper in AIAA 2002 by C.J. Roy et. al.
  err += set_var("Gamma",1.4);
  err += set_var("L",1.0);
  err += set_var("R",287);

  err += set_var("a_px",1);
  err += set_var("a_py",1.25);

  err += set_var("a_rhox",.75);
  err += set_var("a_rhoy",1);

  err += set_var("a_ux",5/3);
  err += set_var("a_uy",1.5);

  err += set_var("a_vx",1.5);
  err += set_var("a_vy",1);

  err += set_var("a_wx",0);
  err += set_var("a_wy",0);

  err += set_var("k",0.0256833);
  err += set_var("mu",1.84e-5);

  err += set_var("p_0",100000);
  err += set_var("p_x",-30000);
  err += set_var("p_y",20000);

  err += set_var("rho_0",1);
  err += set_var("rho_x",0.1);
  err += set_var("rho_y",0.15);

  err += set_var("u_0",70);
  err += set_var("u_x",4);
  err += set_var("u_y",-12);

  err += set_var("v_0",90);
  err += set_var("v_x",-20);
  err += set_var("v_y",4);

  // all Z components are 0
  err += set_var("w_0",0);
  err += set_var("w_x",0);
  err += set_var("w_y",0);

  err += set_var("w_z",0);
  err += set_var("v_z",0);
  err += set_var("u_z",0);
  err += set_var("rho_z",0);
  err += set_var("a_vz",0);
  err += set_var("a_uz",0);
  err += set_var("a_rhoz",0);  
  err += set_var("a_pz",0);
  err += set_var("a_wz",0.0);
  err += set_var("p_z",0);

  /*
  // FIN-S defaults
  err += set_var("Gamma",1.4);
  err += set_var("L",1.0);
  err += set_var("R",286.9);

  //double k1 = R*gamma/(gamma-1)*mu/Pr
  double k1 = 286.9*1.4/(1.4-1)*(1.8e-5/0.71);
  
  err += set_var("k",k1);

  err += set_var("a_px",2);
  err += set_var("a_py",-1.5);

  err += set_var("a_rhox",2);
  err += set_var("a_rhoy",-1.5);

  err += set_var("a_ux",2);
  err += set_var("a_uy",-1.5);

  err += set_var("a_vx",2);
  err += set_var("a_vy",-1.5);

  err += set_var("a_wx",1);
  err += set_var("a_wy",-1);

  err += set_var("mu",1.8e-5);

  err += set_var("p_0",10000);
  err += set_var("p_x",1000);
  err += set_var("p_y",1000);

  err += set_var("rho_0",1);
  err += set_var("rho_x",0.1);
  err += set_var("rho_y",0.1);

  err += set_var("u_0",100);
  err += set_var("u_x",20);
  err += set_var("u_y",10);

  err += set_var("v_0",0);
  err += set_var("v_x",10);
  err += set_var("v_y",10);

  err += set_var("w_0",0);
  err += set_var("w_x",2);
  err += set_var("w_y",-2);

  err += set_var("w_z",4);
  err += set_var("v_z",5);
  err += set_var("u_z",5);
  err += set_var("rho_z",.05);
  err += set_var("a_vz",2);
  err += set_var("a_uz",-2.5);
  err += set_var("a_rhoz",-2);  
  err += set_var("a_pz",-2.5);
  err += set_var("a_wz",2);
  err += set_var("p_z",500);*/

  return err;

} // done with variable initializer

// ----------------------------------------
//   Gradient of Source Terms
// ----------------------------------------

double MASA::navierstokes_3d_compressible::eval_3d_g_u(double x,double y,double z,int i)
{

  double grad = -1;

  switch(i)
    {
    case 0:
      cout << "MASA error:: masa_eval_3d_grad_u has no 0th component\n";
      cout << "Try 1,2 or 3\n";
      break;
      
    case 1:
      grad =  u_x * cos(a_ux * pi * x / L) * a_ux * pi / L;
      
      break;

    case 2:
      grad = -u_y * sin(a_uy * pi * y / L) * a_uy * pi / L;
      break;

    case 3:
      grad = -u_z * sin(a_uz * pi * z / L) * a_uz * pi / L;
      break;

    default:
      cout << "MASA error:: masa_eval_3d_grad_u has no " << i << "th component\n";
      cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

double MASA::navierstokes_3d_compressible::eval_3d_g_v(double x,double y,double z,int i)
{

  double grad = -1;

  switch(i)
    {
    case 0:
      cout << "MASA error:: masa_eval_3d_grad_v has no 0th component\n";
      cout << "Try 1,2 or 3\n";
      break;
      
    case 1:
      grad = -v_x * sin(a_vx * pi * x / L) * a_vx * pi / L;
      break;

    case 2:
      grad =  v_y * cos(a_vy * pi * y / L) * a_vy * pi / L;
      break;

    case 3:
      grad =  v_z * cos(a_vz * pi * z / L) * a_vz * pi / L;      
      break;

    default:
      cout << "MASA error:: masa_eval_3d_grad_v has no " << i << "th component\n";
      cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

double MASA::navierstokes_3d_compressible::eval_3d_g_w(double x,double y,double z,int i)
{

  double grad = -1;

  switch(i)
    {
    case 0:
      cout << "MASA error:: masa_eval_3d_grad_w has no 0th component\n";
      cout << "Try 1,2 or 3\n";
      break;
      
    case 1:
      grad =  w_x * cos(a_wx * pi * x / L) * a_wx * pi / L;
      break;

    case 2:
      grad =  w_y * cos(a_wy * pi * y / L) * a_wy * pi / L;
      break;

    case 3:
      grad = -w_z * sin(a_wz * pi * z / L) * a_wz * pi / L;
      break;

    default:
      cout << "MASA error:: masa_eval_3d_grad_w has no " << i << "th component\n";
      cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}
double MASA::navierstokes_3d_compressible::eval_3d_g_p(double x,double y,double z,int i)
{

  double grad = -1;

  switch(i)
    {
    case 0:
      cout << "MASA error:: masa_eval_3d_grad_p has no 0th component\n";
      cout << "Try 1,2 or 3\n";
      break;
      
    case 1:
      grad = -p_x * sin(a_px * pi * x / L) * a_px * pi / L;
      break;

    case 2:
      grad = p_y * cos(a_py * pi * y / L) * a_py * pi / L;
      break;

    case 3:
      grad = -p_z * sin(a_pz * pi * z / L) * a_pz * pi / L;
      break;

    default:
      cout << "MASA error:: masa_eval_3d_grad_p has no " << i << "th component\n";
      cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

double MASA::navierstokes_3d_compressible::eval_3d_g_rho(double x,double y,double z,int i)
{

  double grad = -1;

  /*
    grad_p_an[0] = -p_x * sin(a_px * pi * x / L) * a_px * pi / L;
    grad_p_an[1] = p_y * cos(a_py * pi * y / L) * a_py * pi / L;
    grad_p_an[2] = -p_z * sin(a_pz * pi * z / L) * a_pz * pi / L;
    grad_v_an[0] = -v_x * sin(a_vx * pi * x / L) * a_vx * pi / L;
    grad_v_an[1] = v_y * cos(a_vy * pi * y / L) * a_vy * pi / L;
    grad_v_an[2] = v_z * cos(a_vz * pi * z / L) * a_vz * pi / L;
    grad_w_an[0] = w_x * cos(a_wx * pi * x / L) * a_wx * pi / L;
    grad_w_an[1] = w_y * cos(a_wy * pi * y / L) * a_wy * pi / L;
    grad_w_an[2] = -w_z * sin(a_wz * pi * z / L) * a_wz * pi / L;
  */

  switch(i)
    {
    case 0:
      cout << "MASA error:: masa_eval_3d_grad_rho has no 0th component\n";
      cout << "Try 1,2 or 3\n";
      break;
      
    case 1:
      grad = rho_x * cos(a_rhox * pi * x / L) * a_rhox * pi / L;
      break;

    case 2:
      grad = -rho_y * sin(a_rhoy * pi * y / L) * a_rhoy * pi / L;
      break;

    case 3:
      grad = rho_z * cos(a_rhoz * pi * z / L) * a_rhoz * pi / L;
      break;

    default:
      cout << "MASA error:: masa_eval_3d_grad_rho has no " << i << "th component\n";
      cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

// ----------------------------------------
//   Source Terms
// ----------------------------------------

double MASA::navierstokes_3d_compressible::eval_q_u(double x,double y,double z)
{
  double Q_u;
  Q_u = 0.4e1 / 0.3e1 * mu * u_x * sin(a_ux * PI * x / L) * a_ux * a_ux * PI * PI * pow(L, -0.2e1) + mu * u_y * cos(a_uy * PI * y / L) * a_uy * a_uy * PI * PI * pow(L, -0.2e1) + mu * u_z * cos(a_uz * PI * z / L) * a_uz * a_uz * PI * PI * pow(L, -0.2e1) - p_x * sin(a_px * PI * x / L) * a_px * PI / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhoz * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_ux * PI / L - u_y * sin(a_uy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_uy * PI / L - u_z * sin(a_uz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_uz * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_vy * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_wz * PI / L;
  return(Q_u);
}

double MASA::navierstokes_3d_compressible::eval_q_v(double x,double y,double z)
{
  double Q_v;
  Q_v = mu * v_x * cos(a_vx * PI * x / L) * a_vx * a_vx * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * mu * v_y * sin(a_vy * PI * y / L) * a_vy * a_vy * PI * PI * pow(L, -0.2e1) + mu * v_z * sin(a_vz * PI * z / L) * a_vz * a_vz * PI * PI * pow(L, -0.2e1) + p_y * cos(a_py * PI * y / L) * a_py * PI / L + rho_x * cos(a_rhox * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_ux * PI / L - v_x * sin(a_vx * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_vx * PI / L + 0.2e1 * v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_vy * PI / L + v_z * cos(a_vz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_vz * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_wz * PI / L;
  return(Q_v);
}

double MASA::navierstokes_3d_compressible::eval_q_w(double x,double y,double z)
{
  double Q_w;
  Q_w = mu * w_x * sin(a_wx * PI * x / L) * a_wx * a_wx * PI * PI * pow(L, -0.2e1) + mu * w_y * sin(a_wy * PI * y / L) * a_wy * a_wy * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * mu * w_z * cos(a_wz * PI * z / L) * a_wz * a_wz * PI * PI * pow(L, -0.2e1) - p_z * sin(a_pz * PI * z / L) * a_pz * PI / L + rho_x * cos(a_rhox * PI * x / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_ux * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_vy * PI / L + w_x * cos(a_wx * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_wx * PI / L + w_y * cos(a_wy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_wy * PI / L - 0.2e1 * w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_wz * PI / L;
  return(Q_w);
}

double MASA::navierstokes_3d_compressible::eval_q_rho(double x,double y,double z)
{
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_ux * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_vy * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_wz * PI / L;
  return(Q_rho);
}

double MASA::navierstokes_3d_compressible::eval_q_e(double x,double y,double z)
{
  double Q_e = cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1)) * rho_x * a_rhox * PI / L / 0.2e1 - sin(a_rhoy * PI * y / L) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1)) * rho_y * a_rhoy * PI / L / 0.2e1 + cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1)) * rho_z * a_rhoz * PI / L / 0.2e1 + ((pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1) + 0.3e1 * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) / L / 0.2e1 + Gamma * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) / L / (Gamma - 0.1e1)) * u_x * cos(a_ux * PI * x / L) * a_ux * PI + ((pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + 0.3e1 * pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) / L / 0.2e1 + Gamma * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) / L / (Gamma - 0.1e1)) * v_y * cos(a_vy * PI * y / L) * a_vy * PI + (-(0.3e1 * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) + pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0, 0.2e1) + pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) / L / 0.2e1 - Gamma * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) / L / (Gamma - 0.1e1)) * w_z * sin(a_wz * PI * z / L) * a_wz * PI + 0.4e1 / 0.3e1 * (-pow(cos(a_ux * PI * x / L), 0.2e1) * u_x + sin(a_ux * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L))) * mu * u_x * a_ux * a_ux * PI * PI * pow(L, -0.2e1) + (-pow(sin(a_uy * PI * y / L), 0.2e1) * u_y + cos(a_uy * PI * y / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L))) * mu * u_y * a_uy * a_uy * PI * PI * pow(L, -0.2e1) + (-pow(sin(a_uz * PI * z / L), 0.2e1) * u_z + cos(a_uz * PI * z / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L))) * mu * u_z * a_uz * a_uz * PI * PI * pow(L, -0.2e1) - (pow(sin(a_vx * PI * x / L), 0.2e1) * v_x - cos(a_vx * PI * x / L) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0)) * mu * v_x * a_vx * a_vx * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * (pow(cos(a_vy * PI * y / L), 0.2e1) * v_y - sin(a_vy * PI * y / L) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0)) * mu * v_y * a_vy * a_vy * PI * PI * pow(L, -0.2e1) - (pow(cos(a_vz * PI * z / L), 0.2e1) * v_z - sin(a_vz * PI * z / L) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0)) * mu * v_z * a_vz * a_vz * PI * PI * pow(L, -0.2e1) + (-pow(cos(a_wx * PI * x / L), 0.2e1) * w_x + sin(a_wx * PI * x / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L))) * mu * w_x * a_wx * a_wx * PI * PI * pow(L, -0.2e1) + (-pow(cos(a_wy * PI * y / L), 0.2e1) * w_y + sin(a_wy * PI * y / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L))) * mu * w_y * a_wy * a_wy * PI * PI * pow(L, -0.2e1) + 0.4e1 / 0.3e1 * (-pow(sin(a_wz * PI * z / L), 0.2e1) * w_z + cos(a_wz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L))) * mu * w_z * a_wz * a_wz * PI * PI * pow(L, -0.2e1) + sin(a_py * PI * y / L) * k * p_y * a_py * a_py * PI * PI * pow(L, -0.2e1) / R / (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) - 0.2e1 * cos(a_rhox * PI * x / L) * rho_x * sin(a_px * PI * x / L) * k * p_x * a_px * a_rhox * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.2e1) - 0.2e1 * sin(a_rhoy * PI * y / L) * rho_y * cos(a_py * PI * y / L) * k * p_y * a_py * a_rhoy * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.2e1) - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * u_y * sin(a_uy * PI * y / L) * a_uy * PI / L + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * v_z * cos(a_vz * PI * z / L) * a_vz * PI / L + cos(a_px * PI * x / L) * k * p_x * a_px * a_px * PI * PI * pow(L, -0.2e1) / R / (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) + cos(a_pz * PI * z / L) * k * p_z * a_pz * a_pz * PI * PI * pow(L, -0.2e1) / R / (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * w_x * cos(a_wx * PI * x / L) * a_wx * PI / L - (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * v_x * sin(a_vx * PI * x / L) * a_vx * PI / L - (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * u_z * sin(a_uz * PI * z / L) * a_uz * PI / L + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * w_y * cos(a_wy * PI * y / L) * a_wy * PI / L - 0.2e1 * cos(a_rhoz * PI * z / L) * rho_z * sin(a_pz * PI * z / L) * k * p_z * a_pz * a_rhoz * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.2e1) - (0.2e1 * pow(cos(a_rhox * PI * x / L), 0.2e1) * rho_x + sin(a_rhox * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L))) * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) * k * rho_x * a_rhox * a_rhox * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.3e1) - (0.2e1 * pow(sin(a_rhoy * PI * y / L), 0.2e1) * rho_y + cos(a_rhoy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L))) * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) * k * rho_y * a_rhoy * a_rhoy * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.3e1) - (0.2e1 * pow(cos(a_rhoz * PI * z / L), 0.2e1) * rho_z + sin(a_rhoz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L))) * (p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L)) * k * rho_z * a_rhoz * a_rhoz * PI * PI * pow(L, -0.2e1) / R * pow(rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L), -0.3e1) + 0.4e1 / 0.3e1 * mu * u_x * v_y * cos(a_ux * PI * x / L) * cos(a_vy * PI * y / L) * a_ux * a_vy * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * u_x * w_z * cos(a_ux * PI * x / L) * sin(a_wz * PI * z / L) * a_ux * a_wz * PI * PI * pow(L, -0.2e1) - 0.2e1 * mu * u_y * v_x * sin(a_uy * PI * y / L) * sin(a_vx * PI * x / L) * a_uy * a_vx * PI * PI * pow(L, -0.2e1) + 0.2e1 * mu * u_z * w_x * cos(a_wx * PI * x / L) * sin(a_uz * PI * z / L) * a_uz * a_wx * PI * PI * pow(L, -0.2e1) - 0.4e1 / 0.3e1 * mu * v_y * w_z * cos(a_vy * PI * y / L) * sin(a_wz * PI * z / L) * a_vy * a_wz * PI * PI * pow(L, -0.2e1) - 0.2e1 * mu * v_z * w_y * cos(a_vz * PI * z / L) * cos(a_wy * PI * y / L) * a_vz * a_wy * PI * PI * pow(L, -0.2e1) - Gamma * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * sin(a_px * PI * x / L) * p_x * a_px * PI / L / (Gamma - 0.1e1) + Gamma * (v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L) + v_0) * cos(a_py * PI * y / L) * p_y * a_py * PI / L / (Gamma - 0.1e1) - Gamma * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * sin(a_pz * PI * z / L) * p_z * a_pz * PI / L / (Gamma - 0.1e1);  
  return(Q_e);
}


// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

double MASA::navierstokes_3d_compressible::eval_an_u(double x,double y,double z)
{
  double u_an;
  u_an = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L);  
  return u_an;
}

double MASA::navierstokes_3d_compressible::eval_an_v(double x,double y,double z)
{
  double v_an;
  v_an = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L);
  return v_an;
}

double MASA::navierstokes_3d_compressible::eval_an_w(double x,double y,double z)
{
  double w_an;
  w_an = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L);
  return w_an;
}

double MASA::navierstokes_3d_compressible::eval_an_p(double x,double y,double z)
{
  double p_an;
  p_an = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L);
  return p_an;
}

double MASA::navierstokes_3d_compressible::eval_an_rho(double x,double y,double z)
{
  double rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L);
  return rho_an;
}
