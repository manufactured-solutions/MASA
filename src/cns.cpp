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

template <typename Scalar>
MASA::navierstokes_2d_compressible<Scalar>::navierstokes_2d_compressible()
{
  this->mmsname = "navierstokes_2d_compressible";
  this->dimension=2;

  this->register_var("R",&R);
  this->register_var("k",&k);

  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("u_y",&u_y);
  this->register_var("v_0",&v_0);
  this->register_var("v_x",&v_x);
  this->register_var("v_y",&v_y);
  this->register_var("rho_0",&rho_0);
  this->register_var("rho_x",&rho_x);
  this->register_var("rho_y",&rho_y);
  this->register_var("p_0",&p_0);
  this->register_var("p_x",&p_x);
  this->register_var("p_y",&p_y);
  this->register_var("a_px",&a_px);
  this->register_var("a_py",&a_py);
  this->register_var("a_rhox",&a_rhox);
  this->register_var("a_rhoy",&a_rhoy);
  this->register_var("a_ux",&a_ux);
  this->register_var("a_uy",&a_uy);
  this->register_var("a_vx",&a_vx);
  this->register_var("a_vy",&a_vy);
  this->register_var("L",&L);
  this->register_var("Gamma",&Gamma);
  this->register_var("mu",&mu);

}//done with constructor

template <typename Scalar>
int MASA::navierstokes_2d_compressible<Scalar>::init_var()
{
  int err = 0;

  // currently randomly generated
  err += this->set_var("R",1.01);
  err += this->set_var("k",1.38);

  err += this->set_var("u_0",1.23);
  err += this->set_var("u_x",1.1);
  err += this->set_var("u_y",.08);
  err += this->set_var("v_0",12);
  err += this->set_var("v_x",1.6);
  err += this->set_var("v_y",.67);
  err += this->set_var("rho_0",1.02);
  err += this->set_var("rho_x",7.2);
  err += this->set_var("rho_y",9.8);
  err += this->set_var("p_0",1.2);
  err += this->set_var("p_x",.91);
  err += this->set_var("p_y",.623);
  err += this->set_var("a_px",.165);
  err += this->set_var("a_py",.612);
  err += this->set_var("a_rhox",.627);
  err += this->set_var("a_rhoy",.828);
  err += this->set_var("a_ux",.1987);
  err += this->set_var("a_uy",1.189);
  err += this->set_var("a_vx",1.91);
  err += this->set_var("a_vy",2.901);
  err += this->set_var("Gamma",1.01);
  err += this->set_var("mu",.918);
  err += this->set_var("L",3.02);

  return err;

} // done with variable initializer

// ----------------------------------------
//   Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_q_u(Scalar x,Scalar y)
{
  Scalar Q_u;
  Q_u = Scalar(4) / Scalar(3) * mu * u_x * sin(a_ux * this->PI * x / L) * a_ux * a_ux * this->PI * this->PI * pow(L, -Scalar(2)) + mu * u_y * cos(a_uy * this->PI * y / L) * a_uy * a_uy * this->PI * this->PI * pow(L, -Scalar(2)) - p_x * sin(a_px * this->PI * x / L) * a_px * this->PI / L + rho_x * cos(a_rhox * this->PI * x / L) * pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L), Scalar(2)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L)) * a_rhoy * this->PI / L + Scalar(2) * u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L)) * a_ux * this->PI / L - u_y * sin(a_uy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L)) * a_uy * this->PI / L + v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L)) * a_vy * this->PI / L;
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_q_v(Scalar x,Scalar y)
{
  Scalar Q_v;
  Q_v = mu * v_x * cos(a_vx * this->PI * x / L) * a_vx * a_vx * this->PI * this->PI * pow(L, -Scalar(2)) + Scalar(4) / Scalar(3) * mu * v_y * sin(a_vy * this->PI * y / L) * a_vy * a_vy * this->PI * this->PI * pow(L, -Scalar(2)) + p_y * cos(a_py * this->PI * y / L) * a_py * this->PI / L + rho_x * cos(a_rhox * this->PI * x / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * pow(v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L), Scalar(2)) * a_rhoy * this->PI / L + u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L)) * a_ux * this->PI / L - v_x * sin(a_vx * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L)) * a_vx * this->PI / L + Scalar(2) * v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L)) * a_vy * this->PI / L;
  return(Q_v);
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_q_rho(Scalar x,Scalar y)
{
  Scalar Q_rho;
  Q_rho = (u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_0) * a_rhox * this->PI * rho_x * cos(a_rhox * this->PI * x / L) / L - (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_0) * a_rhoy * this->PI * rho_y * sin(a_rhoy * this->PI * y / L) / L + (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0) * a_ux * this->PI * u_x * cos(a_ux * this->PI * x / L) / L + (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0) * a_vy * this->PI * v_y * cos(a_vy * this->PI * y / L) / L;
  return(Q_rho);
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_q_e(Scalar x,Scalar y)
{
  Scalar Q_e;
  Q_e = -(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_0) * (pow(u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_0, Scalar(2)) + pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_0, Scalar(2))) * rho_y * sin(a_rhoy * this->PI * y / L) * a_rhoy * this->PI / L / Scalar(2) + (u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_0) * (pow(u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_0, Scalar(2)) + pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_0, Scalar(2))) * rho_x * cos(a_rhox * this->PI * x / L) * a_rhox * this->PI / L / Scalar(2) + Scalar(4) / Scalar(3) * (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_0) * mu * v_y * sin(a_vy * this->PI * y / L) * a_vy * a_vy * this->PI * this->PI * pow(L, -Scalar(2)) - Scalar(4) / Scalar(3) * mu * v_y * v_y * pow(cos(a_vy * this->PI * y / L), Scalar(2)) * a_vy * a_vy * this->PI * this->PI * pow(L, -Scalar(2)) - mu * v_x * v_x * pow(sin(a_vx * this->PI * x / L), Scalar(2)) * a_vx * a_vx * this->PI * this->PI * pow(L, -Scalar(2)) - Scalar(4) / Scalar(3) * mu * u_x * u_x * pow(cos(a_ux * this->PI * x / L), Scalar(2)) * a_ux * a_ux * this->PI * this->PI * pow(L, -Scalar(2)) - mu * u_y * u_y * pow(sin(a_uy * this->PI * y / L), Scalar(2)) * a_uy * a_uy * this->PI * this->PI * pow(L, -Scalar(2)) + (Gamma * (p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_0) / (Gamma - Scalar(1)) + (pow(u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_0, Scalar(2)) / Scalar(2) + Scalar(3) / Scalar(2) * pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_0, Scalar(2))) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0)) * v_y * cos(a_vy * this->PI * y / L) * a_vy * this->PI / L + (Gamma * (p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_0) / (Gamma - Scalar(1)) + (Scalar(3) / Scalar(2) * pow(u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_0, Scalar(2)) + pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_0, Scalar(2)) / Scalar(2)) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0)) * u_x * cos(a_ux * this->PI * x / L) * a_ux * this->PI / L + (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_0) * mu * v_x * cos(a_vx * this->PI * x / L) * a_vx * a_vx * this->PI * this->PI * pow(L, -Scalar(2)) + Scalar(4) / Scalar(3) * (u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_0) * mu * u_x * sin(a_ux * this->PI * x / L) * a_ux * a_ux * this->PI * this->PI * pow(L, -Scalar(2)) + (u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_0) * mu * u_y * cos(a_uy * this->PI * y / L) * a_uy * a_uy * this->PI * this->PI * pow(L, -Scalar(2)) - (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_0) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0) * (u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_0) * u_y * sin(a_uy * this->PI * y / L) * a_uy * this->PI / L - (p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_0) * rho_x * k * sin(a_rhox * this->PI * x / L) * a_rhox * a_rhox * this->PI * this->PI * pow(rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0, -Scalar(2)) * pow(L, -Scalar(2)) / R - (Scalar(2) * p_x * cos(a_px * this->PI * x / L) + Scalar(2) * p_y * sin(a_py * this->PI * y / L) + Scalar(2) * p_0) * rho_x * rho_x * k * pow(cos(a_rhox * this->PI * x / L), Scalar(2)) * a_rhox * a_rhox * this->PI * this->PI * pow(rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0, -Scalar(3)) * pow(L, -Scalar(2)) / R - (p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_0) * rho_y * k * cos(a_rhoy * this->PI * y / L) * a_rhoy * a_rhoy * this->PI * this->PI * pow(rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0, -Scalar(2)) * pow(L, -Scalar(2)) / R - (Scalar(2) * p_x * cos(a_px * this->PI * x / L) + Scalar(2) * p_y * sin(a_py * this->PI * y / L) + Scalar(2) * p_0) * rho_y * rho_y * k * pow(sin(a_rhoy * this->PI * y / L), Scalar(2)) * a_rhoy * a_rhoy * this->PI * this->PI * pow(rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0, -Scalar(3)) * pow(L, -Scalar(2)) / R + Scalar(4) / Scalar(3) * mu * u_x * v_y * cos(a_ux * this->PI * x / L) * cos(a_vy * this->PI * y / L) * a_ux * a_vy * this->PI * this->PI * pow(L, -Scalar(2)) - Scalar(2) * mu * u_y * v_x * sin(a_uy * this->PI * y / L) * sin(a_vx * this->PI * x / L) * a_uy * a_vx * this->PI * this->PI * pow(L, -Scalar(2)) - Scalar(2) * k * p_x * rho_x * cos(a_rhox * this->PI * x / L) * sin(a_px * this->PI * x / L) * a_px * a_rhox * this->PI * this->PI * pow(rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0, -Scalar(2)) * pow(L, -Scalar(2)) / R - Scalar(2) * k * p_y * rho_y * cos(a_py * this->PI * y / L) * sin(a_rhoy * this->PI * y / L) * a_py * a_rhoy * this->PI * this->PI * pow(rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0, -Scalar(2)) * pow(L, -Scalar(2)) / R - (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_0) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0) * (u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_0) * v_x * sin(a_vx * this->PI * x / L) * a_vx * this->PI / L - Gamma * (u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_0) * p_x * sin(a_px * this->PI * x / L) * a_px * this->PI / (Gamma - Scalar(1)) / L + Gamma * (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_0) * p_y * cos(a_py * this->PI * y / L) * a_py * this->PI / (Gamma - Scalar(1)) / L + k * p_x * cos(a_px * this->PI * x / L) * a_px * a_px * this->PI * this->PI / (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0) * pow(L, -Scalar(2)) / R + k * p_y * sin(a_py * this->PI * y / L) * a_py * a_py * this->PI * this->PI / (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0) * pow(L, -Scalar(2)) / R;
  return(Q_e);
}

// ----------------------------------------
//   Gradient of Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_2d_g_u(Scalar x,Scalar y, int i)
{

  Scalar grad = -1;

  switch(i)
    {
    case 0:
      std::cout << "MASA error:: eval_2d_g_u has no 0th component\n";
      std::cout << "Try 1 or 2\n";
      break;
      
    case 1:
      grad =  u_x * cos(a_ux * this->pi * x / L) * a_ux * this->pi / L;      
      break;

    case 2:
      grad = -u_y * sin(a_uy * this->pi * y / L) * a_uy * this->pi / L;
      break;

    default:
      std::cout << "MASA error:: eval_2d_g_u has no " << i << "th component\n";
      std::cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_2d_g_v(Scalar x,Scalar y, int i)
{

  Scalar grad = -1;

  switch(i)
    {
    case 0:
      std::cout << "MASA error:: eval_2d_g_v has no 0th component\n";
      std::cout << "Try 1 or 2\n";
      break;
      
    case 1:
      grad = -v_x * sin(a_vx * this->pi * x / L) * a_vx * this->pi / L;
      break;

    case 2:
      grad = v_y * cos(a_vy * this->pi * y / L) * a_vy * this->pi / L;
      break;

    default:
      std::cout << "MASA error:: eval_2d_g_v has no " << i << "th component\n";
      std::cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_2d_g_p(Scalar x,Scalar y, int i)
{

  Scalar grad = -1;

  switch(i)
    {
    case 0:
      std::cout << "MASA error:: eval_2d_g_p has no 0th component\n";
      std::cout << "Try 1 or 2\n";
      break;
      
    case 1:
      grad = -p_x * sin(a_px * this->pi * x / L) * a_px * this->pi / L;
      break;

    case 2:
      grad =  p_y * cos(a_py * this->pi * y / L) * a_py * this->pi / L;
      break;

    default:
      std::cout << "MASA error:: eval_2d_g_p has no " << i << "th component\n";
      std::cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_2d_g_rho(Scalar x,Scalar y, int i)
{

  Scalar grad = -1;

  switch(i)
    {
    case 0:
      std::cout << "MASA error:: eval_2d_g_rho has no 0th component\n";
      std::cout << "Try 1 or 2\n";
      break;
      
    case 1:
      grad =  rho_x * cos(a_rhox * this->pi * x / L) * a_rhox * this->pi / L;
      break;

    case 2:
      grad = -rho_y * sin(a_rhoy * this->pi * y / L) * a_rhoy * this->pi / L;
      break;

    default:
      std::cout << "MASA error:: eval_2d_g_rho has no " << i << "th component\n";
      std::cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_an_u(Scalar x,Scalar y)
{
  Scalar u_an;
  u_an = u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L);
  return u_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_an_v(Scalar x,Scalar y)
{
  Scalar v_an;
  v_an = v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L);
  return v_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_an_p(Scalar x,Scalar y)
{
  Scalar p_an;
  p_an = p_0 + p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L);
  return p_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_an_rho(Scalar x,Scalar y)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L);
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

template <typename Scalar>
MASA::navierstokes_3d_compressible<Scalar>::navierstokes_3d_compressible()
{
  this->mmsname = "navierstokes_3d_compressible";
  this->dimension=3;

  this->register_var("R",&R);
  this->register_var("k",&k);

  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("u_y",&u_y);
  this->register_var("u_z",&u_z);
  this->register_var("v_0",&v_0);
  this->register_var("v_x",&v_x);
  this->register_var("v_y",&v_y);
  this->register_var("v_z",&v_z);
  this->register_var("rho_0",&rho_0);
  this->register_var("rho_x",&rho_x);
  this->register_var("rho_y",&rho_y);
  this->register_var("rho_z",&rho_z);
  this->register_var("p_0",&p_0);
  this->register_var("p_x",&p_x);
  this->register_var("p_y",&p_y);
  this->register_var("p_z",&p_z);
  this->register_var("a_px",&a_px);
  this->register_var("a_py",&a_py);
  this->register_var("a_rhox",&a_rhox);
  this->register_var("a_rhoy",&a_rhoy);
  this->register_var("a_ux",&a_ux);
  this->register_var("a_uy",&a_uy);
  this->register_var("a_uz",&a_uz);
  this->register_var("a_vx",&a_vx);
  this->register_var("a_vy",&a_vy);
  this->register_var("L",&L);
  this->register_var("Gamma",&Gamma);
  this->register_var("mu",&mu);
  this->register_var("a_pz",&a_pz);
  this->register_var("a_rhoz",&a_rhoz);
  this->register_var("a_vz",&a_vz);
  this->register_var("a_wz",&a_wz);
  this->register_var("a_wx",&a_wx);
  this->register_var("a_wy",&a_wy);
  this->register_var("w_0",&w_0);
  this->register_var("w_x",&w_x);
  this->register_var("w_y",&w_y);
  this->register_var("w_z",&w_z);

}//done with constructor

template <typename Scalar>
int MASA::navierstokes_3d_compressible<Scalar>::init_var()
{
  int err = 0;

  // set using values from table A.5 
  // in a paper in AIAA 2002 by C.J. Roy et. al.
  err += this->set_var("Gamma",1.4);
  err += this->set_var("L",1.0);
  err += this->set_var("R",287);

  err += this->set_var("a_px",1);
  err += this->set_var("a_py",1.25);

  err += this->set_var("a_rhox",.75);
  err += this->set_var("a_rhoy",1);

  err += this->set_var("a_ux",5/3);
  err += this->set_var("a_uy",1.5);

  err += this->set_var("a_vx",1.5);
  err += this->set_var("a_vy",1);

  err += this->set_var("a_wx",0);
  err += this->set_var("a_wy",0);

  err += this->set_var("k",0.0256833);
  err += this->set_var("mu",1.84e-5);

  err += this->set_var("p_0",100000);
  err += this->set_var("p_x",-30000);
  err += this->set_var("p_y",20000);

  err += this->set_var("rho_0",1);
  err += this->set_var("rho_x",0.1);
  err += this->set_var("rho_y",0.15);

  err += this->set_var("u_0",70);
  err += this->set_var("u_x",4);
  err += this->set_var("u_y",-12);

  err += this->set_var("v_0",90);
  err += this->set_var("v_x",-20);
  err += this->set_var("v_y",4);

  // all Z components are 0
  err += this->set_var("w_0",0);
  err += this->set_var("w_x",0);
  err += this->set_var("w_y",0);

  err += this->set_var("w_z",0);
  err += this->set_var("v_z",0);
  err += this->set_var("u_z",0);
  err += this->set_var("rho_z",0);
  err += this->set_var("a_vz",0);
  err += this->set_var("a_uz",0);
  err += this->set_var("a_rhoz",0);  
  err += this->set_var("a_pz",0);
  err += this->set_var("a_wz",0.0);
  err += this->set_var("p_z",0);

  /*
  // FIN-S defaults
  err += this->set_var("Gamma",1.4);
  err += this->set_var("L",1.0);
  err += this->set_var("R",286.9);

  //Scalar k1 = R*gamma/(gamma-1)*mu/Pr
  Scalar k1 = 286.9*1.4/(1.4-1)*(1.8e-5/0.71);
  
  err += this->set_var("k",k1);

  err += this->set_var("a_px",2);
  err += this->set_var("a_py",-1.5);

  err += this->set_var("a_rhox",2);
  err += this->set_var("a_rhoy",-1.5);

  err += this->set_var("a_ux",2);
  err += this->set_var("a_uy",-1.5);

  err += this->set_var("a_vx",2);
  err += this->set_var("a_vy",-1.5);

  err += this->set_var("a_wx",1);
  err += this->set_var("a_wy",-1);

  err += this->set_var("mu",1.8e-5);

  err += this->set_var("p_0",10000);
  err += this->set_var("p_x",1000);
  err += this->set_var("p_y",1000);

  err += this->set_var("rho_0",1);
  err += this->set_var("rho_x",0.1);
  err += this->set_var("rho_y",0.1);

  err += this->set_var("u_0",100);
  err += this->set_var("u_x",20);
  err += this->set_var("u_y",10);

  err += this->set_var("v_0",0);
  err += this->set_var("v_x",10);
  err += this->set_var("v_y",10);

  err += this->set_var("w_0",0);
  err += this->set_var("w_x",2);
  err += this->set_var("w_y",-2);

  err += this->set_var("w_z",4);
  err += this->set_var("v_z",5);
  err += this->set_var("u_z",5);
  err += this->set_var("rho_z",.05);
  err += this->set_var("a_vz",2);
  err += this->set_var("a_uz",-2.5);
  err += this->set_var("a_rhoz",-2);  
  err += this->set_var("a_pz",-2.5);
  err += this->set_var("a_wz",2);
  err += this->set_var("p_z",500);*/

  return err;

} // done with variable initializer

// ----------------------------------------
//   Gradient of Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_3d_g_u(Scalar x,Scalar y,Scalar z,int i)
{

  Scalar grad = -1;

  switch(i)
    {
    case 0:
      std::cout << "MASA error:: masa_eval_3d_grad_u has no 0th component\n";
      std::cout << "Try 1,2 or 3\n";
      break;
      
    case 1:
      grad =  u_x * cos(a_ux * this->pi * x / L) * a_ux * this->pi / L;
      
      break;

    case 2:
      grad = -u_y * sin(a_uy * this->pi * y / L) * a_uy * this->pi / L;
      break;

    case 3:
      grad = -u_z * sin(a_uz * this->pi * z / L) * a_uz * this->pi / L;
      break;

    default:
      std::cout << "MASA error:: masa_eval_3d_grad_u has no " << i << "th component\n";
      std::cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_3d_g_v(Scalar x,Scalar y,Scalar z,int i)
{

  Scalar grad = -1;

  switch(i)
    {
    case 0:
      std::cout << "MASA error:: masa_eval_3d_grad_v has no 0th component\n";
      std::cout << "Try 1,2 or 3\n";
      break;
      
    case 1:
      grad = -v_x * sin(a_vx * this->pi * x / L) * a_vx * this->pi / L;
      break;

    case 2:
      grad =  v_y * cos(a_vy * this->pi * y / L) * a_vy * this->pi / L;
      break;

    case 3:
      grad =  v_z * cos(a_vz * this->pi * z / L) * a_vz * this->pi / L;      
      break;

    default:
      std::cout << "MASA error:: masa_eval_3d_grad_v has no " << i << "th component\n";
      std::cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_3d_g_w(Scalar x,Scalar y,Scalar z,int i)
{

  Scalar grad = -1;

  switch(i)
    {
    case 0:
      std::cout << "MASA error:: masa_eval_3d_grad_w has no 0th component\n";
      std::cout << "Try 1,2 or 3\n";
      break;
      
    case 1:
      grad =  w_x * cos(a_wx * this->pi * x / L) * a_wx * this->pi / L;
      break;

    case 2:
      grad =  w_y * cos(a_wy * this->pi * y / L) * a_wy * this->pi / L;
      break;

    case 3:
      grad = -w_z * sin(a_wz * this->pi * z / L) * a_wz * this->pi / L;
      break;

    default:
      std::cout << "MASA error:: masa_eval_3d_grad_w has no " << i << "th component\n";
      std::cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_3d_g_p(Scalar x,Scalar y,Scalar z,int i)
{

  Scalar grad = -1;

  switch(i)
    {
    case 0:
      std::cout << "MASA error:: masa_eval_3d_grad_p has no 0th component\n";
      std::cout << "Try 1,2 or 3\n";
      break;
      
    case 1:
      grad = -p_x * sin(a_px * this->pi * x / L) * a_px * this->pi / L;
      break;

    case 2:
      grad = p_y * cos(a_py * this->pi * y / L) * a_py * this->pi / L;
      break;

    case 3:
      grad = -p_z * sin(a_pz * this->pi * z / L) * a_pz * this->pi / L;
      break;

    default:
      std::cout << "MASA error:: masa_eval_3d_grad_p has no " << i << "th component\n";
      std::cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_3d_g_rho(Scalar x,Scalar y,Scalar z,int i)
{

  Scalar grad = -1;

  /*
    grad_p_an[0] = -p_x * sin(a_px * this->pi * x / L) * a_px * this->pi / L;
    grad_p_an[1] = p_y * cos(a_py * this->pi * y / L) * a_py * this->pi / L;
    grad_p_an[2] = -p_z * sin(a_pz * this->pi * z / L) * a_pz * this->pi / L;
    grad_v_an[0] = -v_x * sin(a_vx * this->pi * x / L) * a_vx * this->pi / L;
    grad_v_an[1] = v_y * cos(a_vy * this->pi * y / L) * a_vy * this->pi / L;
    grad_v_an[2] = v_z * cos(a_vz * this->pi * z / L) * a_vz * this->pi / L;
    grad_w_an[0] = w_x * cos(a_wx * this->pi * x / L) * a_wx * this->pi / L;
    grad_w_an[1] = w_y * cos(a_wy * this->pi * y / L) * a_wy * this->pi / L;
    grad_w_an[2] = -w_z * sin(a_wz * this->pi * z / L) * a_wz * this->pi / L;
  */

  switch(i)
    {
    case 0:
      std::cout << "MASA error:: masa_eval_3d_grad_rho has no 0th component\n";
      std::cout << "Try 1,2 or 3\n";
      break;
      
    case 1:
      grad = rho_x * cos(a_rhox * this->pi * x / L) * a_rhox * this->pi / L;
      break;

    case 2:
      grad = -rho_y * sin(a_rhoy * this->pi * y / L) * a_rhoy * this->pi / L;
      break;

    case 3:
      grad = rho_z * cos(a_rhoz * this->pi * z / L) * a_rhoz * this->pi / L;
      break;

    default:
      std::cout << "MASA error:: masa_eval_3d_grad_rho has no " << i << "th component\n";
      std::cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

// ----------------------------------------
//   Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_q_u(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_u;
  Q_u = Scalar(4) / Scalar(3) * mu * u_x * sin(a_ux * this->PI * x / L) * a_ux * a_ux * this->PI * this->PI * pow(L, -Scalar(2)) + mu * u_y * cos(a_uy * this->PI * y / L) * a_uy * a_uy * this->PI * this->PI * pow(L, -Scalar(2)) + mu * u_z * cos(a_uz * this->PI * z / L) * a_uz * a_uz * this->PI * this->PI * pow(L, -Scalar(2)) - p_x * sin(a_px * this->PI * x / L) * a_px * this->PI / L + rho_x * cos(a_rhox * this->PI * x / L) * pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_rhoy * this->PI / L + rho_z * cos(a_rhoz * this->PI * z / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_rhoz * this->PI / L + Scalar(2) * u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_ux * this->PI / L - u_y * sin(a_uy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_uy * this->PI / L - u_z * sin(a_uz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_uz * this->PI / L + v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_vy * this->PI / L - w_z * sin(a_wz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_wz * this->PI / L;
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_q_v(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_v;
  Q_v = mu * v_x * cos(a_vx * this->PI * x / L) * a_vx * a_vx * this->PI * this->PI * pow(L, -Scalar(2)) + Scalar(4) / Scalar(3) * mu * v_y * sin(a_vy * this->PI * y / L) * a_vy * a_vy * this->PI * this->PI * pow(L, -Scalar(2)) + mu * v_z * sin(a_vz * this->PI * z / L) * a_vz * a_vz * this->PI * this->PI * pow(L, -Scalar(2)) + p_y * cos(a_py * this->PI * y / L) * a_py * this->PI / L + rho_x * cos(a_rhox * this->PI * x / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * pow(v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L), Scalar(2)) * a_rhoy * this->PI / L + rho_z * cos(a_rhoz * this->PI * z / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_rhoz * this->PI / L + u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_ux * this->PI / L - v_x * sin(a_vx * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_vx * this->PI / L + Scalar(2) * v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_vy * this->PI / L + v_z * cos(a_vz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_vz * this->PI / L - w_z * sin(a_wz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_wz * this->PI / L;
  return(Q_v);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_q_w(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_w;
  Q_w = mu * w_x * sin(a_wx * this->PI * x / L) * a_wx * a_wx * this->PI * this->PI * pow(L, -Scalar(2)) + mu * w_y * sin(a_wy * this->PI * y / L) * a_wy * a_wy * this->PI * this->PI * pow(L, -Scalar(2)) + Scalar(4) / Scalar(3) * mu * w_z * cos(a_wz * this->PI * z / L) * a_wz * a_wz * this->PI * this->PI * pow(L, -Scalar(2)) - p_z * sin(a_pz * this->PI * z / L) * a_pz * this->PI / L + rho_x * cos(a_rhox * this->PI * x / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_rhoy * this->PI / L + rho_z * cos(a_rhoz * this->PI * z / L) * pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2)) * a_rhoz * this->PI / L + u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_ux * this->PI / L + v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_vy * this->PI / L + w_x * cos(a_wx * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_wx * this->PI / L + w_y * cos(a_wy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_wy * this->PI / L - Scalar(2) * w_z * sin(a_wz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_wz * this->PI / L;
  return(Q_w);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_q_rho(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_rho;
  Q_rho = rho_x * cos(a_rhox * this->PI * x / L) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_rhoy * this->PI / L + rho_z * cos(a_rhoz * this->PI * z / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_rhoz * this->PI / L + u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * a_ux * this->PI / L + v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * a_vy * this->PI / L - w_z * sin(a_wz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * a_wz * this->PI / L;
  return(Q_rho);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_q_e(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_e = cos(a_rhox * this->PI * x / L) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * (pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) + pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2)) + pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0, Scalar(2))) * rho_x * a_rhox * this->PI / L / Scalar(2) - sin(a_rhoy * this->PI * y / L) * (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0) * (pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) + pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2)) + pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0, Scalar(2))) * rho_y * a_rhoy * this->PI / L / Scalar(2) + cos(a_rhoz * this->PI * z / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) + pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2)) + pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0, Scalar(2))) * rho_z * a_rhoz * this->PI / L / Scalar(2) + ((pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2)) + pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0, Scalar(2)) + Scalar(3) * pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2))) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) / L / Scalar(2) + Gamma * (p_0 + p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_z * cos(a_pz * this->PI * z / L)) / L / (Gamma - Scalar(1))) * u_x * cos(a_ux * this->PI * x / L) * a_ux * this->PI + ((pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) + pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2)) + Scalar(3) * pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0, Scalar(2))) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) / L / Scalar(2) + Gamma * (p_0 + p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_z * cos(a_pz * this->PI * z / L)) / L / (Gamma - Scalar(1))) * v_y * cos(a_vy * this->PI * y / L) * a_vy * this->PI + (-(Scalar(3) * pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2)) + pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0, Scalar(2)) + pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2))) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) / L / Scalar(2) - Gamma * (p_0 + p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_z * cos(a_pz * this->PI * z / L)) / L / (Gamma - Scalar(1))) * w_z * sin(a_wz * this->PI * z / L) * a_wz * this->PI + Scalar(4) / Scalar(3) * (-pow(cos(a_ux * this->PI * x / L), Scalar(2)) * u_x + sin(a_ux * this->PI * x / L) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L))) * mu * u_x * a_ux * a_ux * this->PI * this->PI * pow(L, -Scalar(2)) + (-pow(sin(a_uy * this->PI * y / L), Scalar(2)) * u_y + cos(a_uy * this->PI * y / L) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L))) * mu * u_y * a_uy * a_uy * this->PI * this->PI * pow(L, -Scalar(2)) + (-pow(sin(a_uz * this->PI * z / L), Scalar(2)) * u_z + cos(a_uz * this->PI * z / L) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L))) * mu * u_z * a_uz * a_uz * this->PI * this->PI * pow(L, -Scalar(2)) - (pow(sin(a_vx * this->PI * x / L), Scalar(2)) * v_x - cos(a_vx * this->PI * x / L) * (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0)) * mu * v_x * a_vx * a_vx * this->PI * this->PI * pow(L, -Scalar(2)) - Scalar(4) / Scalar(3) * (pow(cos(a_vy * this->PI * y / L), Scalar(2)) * v_y - sin(a_vy * this->PI * y / L) * (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0)) * mu * v_y * a_vy * a_vy * this->PI * this->PI * pow(L, -Scalar(2)) - (pow(cos(a_vz * this->PI * z / L), Scalar(2)) * v_z - sin(a_vz * this->PI * z / L) * (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0)) * mu * v_z * a_vz * a_vz * this->PI * this->PI * pow(L, -Scalar(2)) + (-pow(cos(a_wx * this->PI * x / L), Scalar(2)) * w_x + sin(a_wx * this->PI * x / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L))) * mu * w_x * a_wx * a_wx * this->PI * this->PI * pow(L, -Scalar(2)) + (-pow(cos(a_wy * this->PI * y / L), Scalar(2)) * w_y + sin(a_wy * this->PI * y / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L))) * mu * w_y * a_wy * a_wy * this->PI * this->PI * pow(L, -Scalar(2)) + Scalar(4) / Scalar(3) * (-pow(sin(a_wz * this->PI * z / L), Scalar(2)) * w_z + cos(a_wz * this->PI * z / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L))) * mu * w_z * a_wz * a_wz * this->PI * this->PI * pow(L, -Scalar(2)) + sin(a_py * this->PI * y / L) * k * p_y * a_py * a_py * this->PI * this->PI * pow(L, -Scalar(2)) / R / (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) - Scalar(2) * cos(a_rhox * this->PI * x / L) * rho_x * sin(a_px * this->PI * x / L) * k * p_x * a_px * a_rhox * this->PI * this->PI * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L), -Scalar(2)) - Scalar(2) * sin(a_rhoy * this->PI * y / L) * rho_y * cos(a_py * this->PI * y / L) * k * p_y * a_py * a_rhoy * this->PI * this->PI * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L), -Scalar(2)) - (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * u_y * sin(a_uy * this->PI * y / L) * a_uy * this->PI / L + (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0) * v_z * cos(a_vz * this->PI * z / L) * a_vz * this->PI / L + cos(a_px * this->PI * x / L) * k * p_x * a_px * a_px * this->PI * this->PI * pow(L, -Scalar(2)) / R / (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) + cos(a_pz * this->PI * z / L) * k * p_z * a_pz * a_pz * this->PI * this->PI * pow(L, -Scalar(2)) / R / (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) + (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * w_x * cos(a_wx * this->PI * x / L) * a_wx * this->PI / L - (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * v_x * sin(a_vx * this->PI * x / L) * a_vx * this->PI / L - (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * u_z * sin(a_uz * this->PI * z / L) * a_uz * this->PI / L + (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0) * w_y * cos(a_wy * this->PI * y / L) * a_wy * this->PI / L - Scalar(2) * cos(a_rhoz * this->PI * z / L) * rho_z * sin(a_pz * this->PI * z / L) * k * p_z * a_pz * a_rhoz * this->PI * this->PI * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L), -Scalar(2)) - (Scalar(2) * pow(cos(a_rhox * this->PI * x / L), Scalar(2)) * rho_x + sin(a_rhox * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L))) * (p_0 + p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_z * cos(a_pz * this->PI * z / L)) * k * rho_x * a_rhox * a_rhox * this->PI * this->PI * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L), -Scalar(3)) - (Scalar(2) * pow(sin(a_rhoy * this->PI * y / L), Scalar(2)) * rho_y + cos(a_rhoy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L))) * (p_0 + p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_z * cos(a_pz * this->PI * z / L)) * k * rho_y * a_rhoy * a_rhoy * this->PI * this->PI * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L), -Scalar(3)) - (Scalar(2) * pow(cos(a_rhoz * this->PI * z / L), Scalar(2)) * rho_z + sin(a_rhoz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L))) * (p_0 + p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_z * cos(a_pz * this->PI * z / L)) * k * rho_z * a_rhoz * a_rhoz * this->PI * this->PI * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L), -Scalar(3)) + Scalar(4) / Scalar(3) * mu * u_x * v_y * cos(a_ux * this->PI * x / L) * cos(a_vy * this->PI * y / L) * a_ux * a_vy * this->PI * this->PI * pow(L, -Scalar(2)) - Scalar(4) / Scalar(3) * mu * u_x * w_z * cos(a_ux * this->PI * x / L) * sin(a_wz * this->PI * z / L) * a_ux * a_wz * this->PI * this->PI * pow(L, -Scalar(2)) - Scalar(2) * mu * u_y * v_x * sin(a_uy * this->PI * y / L) * sin(a_vx * this->PI * x / L) * a_uy * a_vx * this->PI * this->PI * pow(L, -Scalar(2)) + Scalar(2) * mu * u_z * w_x * cos(a_wx * this->PI * x / L) * sin(a_uz * this->PI * z / L) * a_uz * a_wx * this->PI * this->PI * pow(L, -Scalar(2)) - Scalar(4) / Scalar(3) * mu * v_y * w_z * cos(a_vy * this->PI * y / L) * sin(a_wz * this->PI * z / L) * a_vy * a_wz * this->PI * this->PI * pow(L, -Scalar(2)) - Scalar(2) * mu * v_z * w_y * cos(a_vz * this->PI * z / L) * cos(a_wy * this->PI * y / L) * a_vz * a_wy * this->PI * this->PI * pow(L, -Scalar(2)) - Gamma * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * sin(a_px * this->PI * x / L) * p_x * a_px * this->PI / L / (Gamma - Scalar(1)) + Gamma * (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L) + v_0) * cos(a_py * this->PI * y / L) * p_y * a_py * this->PI / L / (Gamma - Scalar(1)) - Gamma * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * sin(a_pz * this->PI * z / L) * p_z * a_pz * this->PI / L / (Gamma - Scalar(1));  
  return(Q_e);
}


// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_an_u(Scalar x,Scalar y,Scalar z)
{
  Scalar u_an;
  u_an = u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L);  
  return u_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_an_v(Scalar x,Scalar y,Scalar z)
{
  Scalar v_an;
  v_an = v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L);
  return v_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_an_w(Scalar x,Scalar y,Scalar z)
{
  Scalar w_an;
  w_an = w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L);
  return w_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_an_p(Scalar x,Scalar y,Scalar z)
{
  Scalar p_an;
  p_an = p_0 + p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_z * cos(a_pz * this->PI * z / L);
  return p_an;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_an_rho(Scalar x,Scalar y,Scalar z)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L);
  return rho_an;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::navierstokes_2d_compressible);
MASA_INSTANTIATE_ALL(MASA::navierstokes_3d_compressible);
