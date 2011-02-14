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

  this->init_var();

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
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_q_rho_u(Scalar x,Scalar y)
{
  Scalar Q_u;
  Q_u = Scalar(4) / Scalar(3) * mu * u_x * sin(a_ux * pi * x / L) * a_ux * a_ux * pi * pi * pow(L, -Scalar(2)) + mu * u_y * cos(a_uy * pi * y / L) * a_uy * a_uy * pi * pi * pow(L, -Scalar(2)) - p_x * sin(a_px * pi * x / L) * a_px * pi / L + rho_x * cos(a_rhox * pi * x / L) * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L), Scalar(2)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_rhoy * pi / L + Scalar(2) * u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_ux * pi / L - u_y * sin(a_uy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_uy * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_vy * pi / L;
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_q_rho_v(Scalar x,Scalar y)
{
  Scalar Q_v;
  Q_v = mu * v_x * cos(a_vx * pi * x / L) * a_vx * a_vx * pi * pi * pow(L, -Scalar(2)) + Scalar(4) / Scalar(3) * mu * v_y * sin(a_vy * pi * y / L) * a_vy * a_vy * pi * pi * pow(L, -Scalar(2)) + p_y * cos(a_py * pi * y / L) * a_py * pi / L + rho_x * cos(a_rhox * pi * x / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L), Scalar(2)) * a_rhoy * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_ux * pi / L - v_x * sin(a_vx * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L)) * a_vx * pi / L + Scalar(2) * v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L)) * a_vy * pi / L;
  return(Q_v);
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_q_rho(Scalar x,Scalar y)
{
  Scalar Q_rho;
  Q_rho = (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * a_rhox * pi * rho_x * cos(a_rhox * pi * x / L) / L - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * a_rhoy * pi * rho_y * sin(a_rhoy * pi * y / L) / L + (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * a_ux * pi * u_x * cos(a_ux * pi * x / L) / L + (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * a_vy * pi * v_y * cos(a_vy * pi * y / L) / L;
  return(Q_rho);
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_q_e(Scalar x,Scalar y)
{
  Scalar Q_e;
  Q_e = -(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, Scalar(2)) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, Scalar(2))) * rho_y * sin(a_rhoy * pi * y / L) * a_rhoy * pi / L / Scalar(2) + (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, Scalar(2)) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, Scalar(2))) * rho_x * cos(a_rhox * pi * x / L) * a_rhox * pi / L / Scalar(2) + Scalar(4) / Scalar(3) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * mu * v_y * sin(a_vy * pi * y / L) * a_vy * a_vy * pi * pi * pow(L, -Scalar(2)) - Scalar(4) / Scalar(3) * mu * v_y * v_y * pow(cos(a_vy * pi * y / L), Scalar(2)) * a_vy * a_vy * pi * pi * pow(L, -Scalar(2)) - mu * v_x * v_x * pow(sin(a_vx * pi * x / L), Scalar(2)) * a_vx * a_vx * pi * pi * pow(L, -Scalar(2)) - Scalar(4) / Scalar(3) * mu * u_x * u_x * pow(cos(a_ux * pi * x / L), Scalar(2)) * a_ux * a_ux * pi * pi * pow(L, -Scalar(2)) - mu * u_y * u_y * pow(sin(a_uy * pi * y / L), Scalar(2)) * a_uy * a_uy * pi * pi * pow(L, -Scalar(2)) + (Gamma * (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) / (Gamma - Scalar(1)) + (pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, Scalar(2)) / Scalar(2) + Scalar(3) / Scalar(2) * pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, Scalar(2))) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0)) * v_y * cos(a_vy * pi * y / L) * a_vy * pi / L + (Gamma * (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) / (Gamma - Scalar(1)) + (Scalar(3) / Scalar(2) * pow(u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0, Scalar(2)) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0, Scalar(2)) / Scalar(2)) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0)) * u_x * cos(a_ux * pi * x / L) * a_ux * pi / L + (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * mu * v_x * cos(a_vx * pi * x / L) * a_vx * a_vx * pi * pi * pow(L, -Scalar(2)) + Scalar(4) / Scalar(3) * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * mu * u_x * sin(a_ux * pi * x / L) * a_ux * a_ux * pi * pi * pow(L, -Scalar(2)) + (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * mu * u_y * cos(a_uy * pi * y / L) * a_uy * a_uy * pi * pi * pow(L, -Scalar(2)) - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * u_y * sin(a_uy * pi * y / L) * a_uy * pi / L - (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) * rho_x * k * sin(a_rhox * pi * x / L) * a_rhox * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -Scalar(2)) * pow(L, -Scalar(2)) / R - (Scalar(2) * p_x * cos(a_px * pi * x / L) + Scalar(2) * p_y * sin(a_py * pi * y / L) + Scalar(2) * p_0) * rho_x * rho_x * k * pow(cos(a_rhox * pi * x / L), Scalar(2)) * a_rhox * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -Scalar(3)) * pow(L, -Scalar(2)) / R - (p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_0) * rho_y * k * cos(a_rhoy * pi * y / L) * a_rhoy * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -Scalar(2)) * pow(L, -Scalar(2)) / R - (Scalar(2) * p_x * cos(a_px * pi * x / L) + Scalar(2) * p_y * sin(a_py * pi * y / L) + Scalar(2) * p_0) * rho_y * rho_y * k * pow(sin(a_rhoy * pi * y / L), Scalar(2)) * a_rhoy * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -Scalar(3)) * pow(L, -Scalar(2)) / R + Scalar(4) / Scalar(3) * mu * u_x * v_y * cos(a_ux * pi * x / L) * cos(a_vy * pi * y / L) * a_ux * a_vy * pi * pi * pow(L, -Scalar(2)) - Scalar(2) * mu * u_y * v_x * sin(a_uy * pi * y / L) * sin(a_vx * pi * x / L) * a_uy * a_vx * pi * pi * pow(L, -Scalar(2)) - Scalar(2) * k * p_x * rho_x * cos(a_rhox * pi * x / L) * sin(a_px * pi * x / L) * a_px * a_rhox * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -Scalar(2)) * pow(L, -Scalar(2)) / R - Scalar(2) * k * p_y * rho_y * cos(a_py * pi * y / L) * sin(a_rhoy * pi * y / L) * a_py * a_rhoy * pi * pi * pow(rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0, -Scalar(2)) * pow(L, -Scalar(2)) / R - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * v_x * sin(a_vx * pi * x / L) * a_vx * pi / L - Gamma * (u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_0) * p_x * sin(a_px * pi * x / L) * a_px * pi / (Gamma - Scalar(1)) / L + Gamma * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_0) * p_y * cos(a_py * pi * y / L) * a_py * pi / (Gamma - Scalar(1)) / L + k * p_x * cos(a_px * pi * x / L) * a_px * a_px * pi * pi / (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * pow(L, -Scalar(2)) / R + k * p_y * sin(a_py * pi * y / L) * a_py * a_py * pi * pi / (rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_0) * pow(L, -Scalar(2)) / R;
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
      grad =  u_x * cos(a_ux * pi * x / L) * a_ux * pi / L;      
      break;

    case 2:
      grad = -u_y * sin(a_uy * pi * y / L) * a_uy * pi / L;
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
      grad = -v_x * sin(a_vx * pi * x / L) * a_vx * pi / L;
      break;

    case 2:
      grad = v_y * cos(a_vy * pi * y / L) * a_vy * pi / L;
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
      grad = -p_x * sin(a_px * pi * x / L) * a_px * pi / L;
      break;

    case 2:
      grad =  p_y * cos(a_py * pi * y / L) * a_py * pi / L;
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
      grad =  rho_x * cos(a_rhox * pi * x / L) * a_rhox * pi / L;
      break;

    case 2:
      grad = -rho_y * sin(a_rhoy * pi * y / L) * a_rhoy * pi / L;
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
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_exact_u(Scalar x,Scalar y)
{
  Scalar exact_u;
  exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return exact_u;
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_exact_v(Scalar x,Scalar y)
{
  Scalar exact_v;
  exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return exact_v;
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_exact_p(Scalar x,Scalar y)
{
  Scalar exact_p;
  exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return exact_p;
}

template <typename Scalar>
Scalar MASA::navierstokes_2d_compressible<Scalar>::eval_exact_rho(Scalar x,Scalar y)
{
  Scalar exact_rho;
  exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return exact_rho;
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

  this->init_var();

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
      std::cout << "MASA error:: masa_eval_grad_3d_u has no 0th component\n";
      std::cout << "Try 1,2 or 3\n";
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
      std::cout << "MASA error:: masa_eval_grad_3d_u has no " << i << "th component\n";
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
      std::cout << "MASA error:: masa_eval_grad_3d_v has no 0th component\n";
      std::cout << "Try 1,2 or 3\n";
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
      std::cout << "MASA error:: masa_eval_grad_3d_v has no " << i << "th component\n";
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
      std::cout << "MASA error:: masa_eval_grad_3d_w has no 0th component\n";
      std::cout << "Try 1,2 or 3\n";
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
      std::cout << "MASA error:: masa_eval_grad_3d_w has no " << i << "th component\n";
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
      std::cout << "MASA error:: masa_eval_grad_3d_p has no 0th component\n";
      std::cout << "Try 1,2 or 3\n";
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
      std::cout << "MASA error:: masa_eval_grad_3d_p has no " << i << "th component\n";
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
    grad_exact_p[0] = -p_x * sin(a_px * pi * x / L) * a_px * pi / L;
    grad_exact_p[1] = p_y * cos(a_py * pi * y / L) * a_py * pi / L;
    grad_exact_p[2] = -p_z * sin(a_pz * pi * z / L) * a_pz * pi / L;
    grad_exact_v[0] = -v_x * sin(a_vx * pi * x / L) * a_vx * pi / L;
    grad_exact_v[1] = v_y * cos(a_vy * pi * y / L) * a_vy * pi / L;
    grad_exact_v[2] = v_z * cos(a_vz * pi * z / L) * a_vz * pi / L;
    grad_exact_w[0] = w_x * cos(a_wx * pi * x / L) * a_wx * pi / L;
    grad_exact_w[1] = w_y * cos(a_wy * pi * y / L) * a_wy * pi / L;
    grad_exact_w[2] = -w_z * sin(a_wz * pi * z / L) * a_wz * pi / L;
  */

  switch(i)
    {
    case 0:
      std::cout << "MASA error:: masa_eval_grad_3d_rho has no 0th component\n";
      std::cout << "Try 1,2 or 3\n";
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
      std::cout << "MASA error:: masa_eval_grad_3d_rho has no " << i << "th component\n";
      std::cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

// ----------------------------------------
//   Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_q_rho_u(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_u;
  Q_u = Scalar(4) / Scalar(3) * mu * u_x * sin(a_ux * pi * x / L) * a_ux * a_ux * pi * pi * pow(L, -Scalar(2)) + mu * u_y * cos(a_uy * pi * y / L) * a_uy * a_uy * pi * pi * pow(L, -Scalar(2)) + mu * u_z * cos(a_uz * pi * z / L) * a_uz * a_uz * pi * pi * pow(L, -Scalar(2)) - p_x * sin(a_px * pi * x / L) * a_px * pi / L + rho_x * cos(a_rhox * pi * x / L) * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), Scalar(2)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhoz * pi / L + Scalar(2) * u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_ux * pi / L - u_y * sin(a_uy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_uy * pi / L - u_z * sin(a_uz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_uz * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_vy * pi / L - w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_wz * pi / L;
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_q_rho_v(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_v;
  Q_v = mu * v_x * cos(a_vx * pi * x / L) * a_vx * a_vx * pi * pi * pow(L, -Scalar(2)) + Scalar(4) / Scalar(3) * mu * v_y * sin(a_vy * pi * y / L) * a_vy * a_vy * pi * pi * pow(L, -Scalar(2)) + mu * v_z * sin(a_vz * pi * z / L) * a_vz * a_vz * pi * pi * pow(L, -Scalar(2)) + p_y * cos(a_py * pi * y / L) * a_py * pi / L + rho_x * cos(a_rhox * pi * x / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * pow(v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L), Scalar(2)) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_rhoz * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_ux * pi / L - v_x * sin(a_vx * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_vx * pi / L + Scalar(2) * v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_vy * pi / L + v_z * cos(a_vz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_vz * pi / L - w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_wz * pi / L;
  return(Q_v);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_q_rho_w(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_w;
  Q_w = mu * w_x * sin(a_wx * pi * x / L) * a_wx * a_wx * pi * pi * pow(L, -Scalar(2)) + mu * w_y * sin(a_wy * pi * y / L) * a_wy * a_wy * pi * pi * pow(L, -Scalar(2)) + Scalar(4) / Scalar(3) * mu * w_z * cos(a_wz * pi * z / L) * a_wz * a_wz * pi * pi * pow(L, -Scalar(2)) - p_z * sin(a_pz * pi * z / L) * a_pz * pi / L + rho_x * cos(a_rhox * pi * x / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), Scalar(2)) * a_rhoz * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_ux * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_vy * pi / L + w_x * cos(a_wx * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_wx * pi / L + w_y * cos(a_wy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_wy * pi / L - Scalar(2) * w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_wz * pi / L;
  return(Q_w);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_q_rho(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_rho;
  Q_rho = rho_x * cos(a_rhox * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * a_rhox * pi / L - rho_y * sin(a_rhoy * pi * y / L) * (v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L)) * a_rhoy * pi / L + rho_z * cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * a_rhoz * pi / L + u_x * cos(a_ux * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * a_ux * pi / L + v_y * cos(a_vy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * a_vy * pi / L - w_z * sin(a_wz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * a_wz * pi / L;
  return(Q_rho);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_q_e(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_e = cos(a_rhox * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * (pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), Scalar(2)) + pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), Scalar(2)) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, Scalar(2))) * rho_x * a_rhox * pi / L / Scalar(2) - sin(a_rhoy * pi * y / L) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * (pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), Scalar(2)) + pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), Scalar(2)) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, Scalar(2))) * rho_y * a_rhoy * pi / L / Scalar(2) + cos(a_rhoz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), Scalar(2)) + pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), Scalar(2)) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, Scalar(2))) * rho_z * a_rhoz * pi / L / Scalar(2) + ((pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), Scalar(2)) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, Scalar(2)) + Scalar(3) * pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), Scalar(2))) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) / L / Scalar(2) + Gamma * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) / L / (Gamma - Scalar(1))) * u_x * cos(a_ux * pi * x / L) * a_ux * pi + ((pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), Scalar(2)) + pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), Scalar(2)) + Scalar(3) * pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, Scalar(2))) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) / L / Scalar(2) + Gamma * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) / L / (Gamma - Scalar(1))) * v_y * cos(a_vy * pi * y / L) * a_vy * pi + (-(Scalar(3) * pow(w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L), Scalar(2)) + pow(v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0, Scalar(2)) + pow(u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L), Scalar(2))) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) / L / Scalar(2) - Gamma * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) / L / (Gamma - Scalar(1))) * w_z * sin(a_wz * pi * z / L) * a_wz * pi + Scalar(4) / Scalar(3) * (-pow(cos(a_ux * pi * x / L), Scalar(2)) * u_x + sin(a_ux * pi * x / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L))) * mu * u_x * a_ux * a_ux * pi * pi * pow(L, -Scalar(2)) + (-pow(sin(a_uy * pi * y / L), Scalar(2)) * u_y + cos(a_uy * pi * y / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L))) * mu * u_y * a_uy * a_uy * pi * pi * pow(L, -Scalar(2)) + (-pow(sin(a_uz * pi * z / L), Scalar(2)) * u_z + cos(a_uz * pi * z / L) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L))) * mu * u_z * a_uz * a_uz * pi * pi * pow(L, -Scalar(2)) - (pow(sin(a_vx * pi * x / L), Scalar(2)) * v_x - cos(a_vx * pi * x / L) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0)) * mu * v_x * a_vx * a_vx * pi * pi * pow(L, -Scalar(2)) - Scalar(4) / Scalar(3) * (pow(cos(a_vy * pi * y / L), Scalar(2)) * v_y - sin(a_vy * pi * y / L) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0)) * mu * v_y * a_vy * a_vy * pi * pi * pow(L, -Scalar(2)) - (pow(cos(a_vz * pi * z / L), Scalar(2)) * v_z - sin(a_vz * pi * z / L) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0)) * mu * v_z * a_vz * a_vz * pi * pi * pow(L, -Scalar(2)) + (-pow(cos(a_wx * pi * x / L), Scalar(2)) * w_x + sin(a_wx * pi * x / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L))) * mu * w_x * a_wx * a_wx * pi * pi * pow(L, -Scalar(2)) + (-pow(cos(a_wy * pi * y / L), Scalar(2)) * w_y + sin(a_wy * pi * y / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L))) * mu * w_y * a_wy * a_wy * pi * pi * pow(L, -Scalar(2)) + Scalar(4) / Scalar(3) * (-pow(sin(a_wz * pi * z / L), Scalar(2)) * w_z + cos(a_wz * pi * z / L) * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L))) * mu * w_z * a_wz * a_wz * pi * pi * pow(L, -Scalar(2)) + sin(a_py * pi * y / L) * k * p_y * a_py * a_py * pi * pi * pow(L, -Scalar(2)) / R / (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) - Scalar(2) * cos(a_rhox * pi * x / L) * rho_x * sin(a_px * pi * x / L) * k * p_x * a_px * a_rhox * pi * pi * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -Scalar(2)) - Scalar(2) * sin(a_rhoy * pi * y / L) * rho_y * cos(a_py * pi * y / L) * k * p_y * a_py * a_rhoy * pi * pi * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -Scalar(2)) - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * u_y * sin(a_uy * pi * y / L) * a_uy * pi / L + (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * v_z * cos(a_vz * pi * z / L) * a_vz * pi / L + cos(a_px * pi * x / L) * k * p_x * a_px * a_px * pi * pi * pow(L, -Scalar(2)) / R / (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) + cos(a_pz * pi * z / L) * k * p_z * a_pz * a_pz * pi * pi * pow(L, -Scalar(2)) / R / (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) + (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * w_x * cos(a_wx * pi * x / L) * a_wx * pi / L - (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * v_x * sin(a_vx * pi * x / L) * a_vx * pi / L - (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * u_z * sin(a_uz * pi * z / L) * a_uz * pi / L + (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L)) * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * w_y * cos(a_wy * pi * y / L) * a_wy * pi / L - Scalar(2) * cos(a_rhoz * pi * z / L) * rho_z * sin(a_pz * pi * z / L) * k * p_z * a_pz * a_rhoz * pi * pi * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -Scalar(2)) - (Scalar(2) * pow(cos(a_rhox * pi * x / L), Scalar(2)) * rho_x + sin(a_rhox * pi * x / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L))) * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) * k * rho_x * a_rhox * a_rhox * pi * pi * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -Scalar(3)) - (Scalar(2) * pow(sin(a_rhoy * pi * y / L), Scalar(2)) * rho_y + cos(a_rhoy * pi * y / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L))) * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) * k * rho_y * a_rhoy * a_rhoy * pi * pi * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -Scalar(3)) - (Scalar(2) * pow(cos(a_rhoz * pi * z / L), Scalar(2)) * rho_z + sin(a_rhoz * pi * z / L) * (rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L))) * (p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L)) * k * rho_z * a_rhoz * a_rhoz * pi * pi * pow(L, -Scalar(2)) / R * pow(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L), -Scalar(3)) + Scalar(4) / Scalar(3) * mu * u_x * v_y * cos(a_ux * pi * x / L) * cos(a_vy * pi * y / L) * a_ux * a_vy * pi * pi * pow(L, -Scalar(2)) - Scalar(4) / Scalar(3) * mu * u_x * w_z * cos(a_ux * pi * x / L) * sin(a_wz * pi * z / L) * a_ux * a_wz * pi * pi * pow(L, -Scalar(2)) - Scalar(2) * mu * u_y * v_x * sin(a_uy * pi * y / L) * sin(a_vx * pi * x / L) * a_uy * a_vx * pi * pi * pow(L, -Scalar(2)) + Scalar(2) * mu * u_z * w_x * cos(a_wx * pi * x / L) * sin(a_uz * pi * z / L) * a_uz * a_wx * pi * pi * pow(L, -Scalar(2)) - Scalar(4) / Scalar(3) * mu * v_y * w_z * cos(a_vy * pi * y / L) * sin(a_wz * pi * z / L) * a_vy * a_wz * pi * pi * pow(L, -Scalar(2)) - Scalar(2) * mu * v_z * w_y * cos(a_vz * pi * z / L) * cos(a_wy * pi * y / L) * a_vz * a_wy * pi * pi * pow(L, -Scalar(2)) - Gamma * (u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L)) * sin(a_px * pi * x / L) * p_x * a_px * pi / L / (Gamma - Scalar(1)) + Gamma * (v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_0) * cos(a_py * pi * y / L) * p_y * a_py * pi / L / (Gamma - Scalar(1)) - Gamma * (w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L)) * sin(a_pz * pi * z / L) * p_z * a_pz * pi / L / (Gamma - Scalar(1));  
  return(Q_e);
}


// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_exact_u(Scalar x,Scalar y,Scalar z)
{
  Scalar exact_u;
  exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);  
  return exact_u;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_exact_v(Scalar x,Scalar y,Scalar z)
{
  Scalar exact_v;
  exact_v = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  return exact_v;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_exact_w(Scalar x,Scalar y,Scalar z)
{
  Scalar exact_w;
  exact_w = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);
  return exact_w;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_exact_p(Scalar x,Scalar y,Scalar z)
{
  Scalar exact_p;
  exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L);
  return exact_p;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_compressible<Scalar>::eval_exact_rho(Scalar x,Scalar y,Scalar z)
{
  Scalar exact_rho;
  exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  return exact_rho;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::navierstokes_2d_compressible);
MASA_INSTANTIATE_ALL(MASA::navierstokes_3d_compressible);
