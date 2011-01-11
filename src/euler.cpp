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
// euler.cpp: These are the MASA class member functions and constructors
//          For the Euler Equations
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         EULER EQUATION 1D
 *
 *
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::euler_1d<Scalar>::euler_1d()
{
  this->mmsname = "euler_1d";
  this->dimension=1;

  this->register_var("R",&R);
  this->register_var("k",&k);

  this->register_var("u_0",&u_0);
  this->register_var("u_x",&u_x);
  this->register_var("rho_0",&rho_0);
  this->register_var("rho_x",&rho_x);
  this->register_var("p_0",&p_0);
  this->register_var("p_x",&p_x);
  this->register_var("a_px",&a_px);
  this->register_var("a_rhox",&a_rhox);
  this->register_var("a_ux",&a_ux);
  this->register_var("L",&L);
  this->register_var("Gamma",&Gamma);
  this->register_var("mu",&mu);

}//done with constructor

template <typename Scalar>
int MASA::euler_1d<Scalar>::init_var()
{
  int err = 0;

  // randomly generated
  err += this->set_var("R",1.01);
  err += this->set_var("k",1.38);

  err += this->set_var("u_0",.191);
  err += this->set_var("u_x",1.63);
  err += this->set_var("rho_0",91.5);
  err += this->set_var("rho_x",5.13);
  err += this->set_var("p_0",.1984);
  err += this->set_var("p_x",3.151);
  err += this->set_var("a_px",6.151);
  err += this->set_var("a_rhox",1.2);
  err += this->set_var("a_ux",.03);
  err += this->set_var("L",3.02);
  err += this->set_var("Gamma",16.1);
  err += this->set_var("mu",.091);

  return err;

} // done with variable initializer


// ----------------------------------------
//   Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_q_u(Scalar x)
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

  Scalar Q_u;
  Q_u = -sin(a_px * PI * x / L) * a_px * this->PI * p_x / L + rho_x * cos(a_rhox * this->PI * x / L) * pow(u_0 + u_x * sin(a_ux * this->PI * x / L), Scalar(2)) * a_rhox * this->PI / L + Scalar(2) * u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L)) * a_ux * this->PI / L;
  return Q_u;
}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_q_e(Scalar x)
{
  Scalar Q_e;
  Q_e = cos(a_rhox * this->PI * x / L) * rho_x * pow(u_0 + u_x * sin(a_ux * this->PI * x / L), Scalar(3)) * a_rhox * this->PI / L / Scalar(2) + cos(a_ux * this->PI * x / L) * (p_0 + p_x * cos(a_px * this->PI * x / L)) * a_ux * this->PI * u_x * Gamma / L / (Gamma - Scalar(1)) - Gamma * p_x * sin(a_px * this->PI * x / L) * (u_0 + u_x * sin(a_ux * this->PI * x / L)) * a_px * this->PI / L / (Gamma - Scalar(1)) + Scalar(3) / Scalar(2) * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L)) * pow(u_0 + u_x * sin(a_ux * this->PI * x / L), Scalar(2)) * a_ux * this->PI * u_x / L;
  return(Q_e);
}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_q_rho(Scalar x)
{
  Scalar Q_rho;
  Q_rho = rho_x * cos(a_rhox * this->PI * x / L) * (u_0 + u_x * sin(a_ux * this->PI * x / L)) * a_rhox * this->PI / L + u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L)) * a_ux * this->PI / L;
  return(Q_rho);
}

// ----------------------------------------
//   Gradient of Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_1d_g_u(Scalar x)
{

  Scalar grad_u = u_x * cos(a_ux * this->pi * x / L) * a_ux * this->pi / L;
  return grad_u;


}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_1d_g_p(Scalar x)
{

  Scalar grad_p = -p_x * sin(a_px * this->pi * x / L) * a_px * this->pi / L;
  return grad_p;

}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_1d_g_rho(Scalar x)
{

  Scalar grad_rho = rho_x * cos(a_rhox * this->pi * x / L) * a_rhox * this->pi / L;
  return grad_rho;

}


// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_an_u(Scalar x)
{
  Scalar u_an;
  u_an = u_0 + u_x * sin(a_ux * this->PI * x / L);
  return u_an;
}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_an_p(Scalar x)
{
  Scalar p_an;
  p_an = p_0 + p_x * cos(a_px * this->PI * x / L);
  return p_an;
}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_an_rho(Scalar x)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * this->PI * x / L);                                                                                                                                                                                           
  return rho_an;
}


/* ------------------------------------------------
 *
 *         EULER EQUATION 2D
 *
 *
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::euler_2d<Scalar>::euler_2d()
{
  this->mmsname = "euler_2d";
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
int MASA::euler_2d<Scalar>::init_var()
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
Scalar MASA::euler_2d<Scalar>::eval_q_u(Scalar x,Scalar y)
{
  Scalar Q_u;
  Q_u = -p_x * sin(a_px * this->PI * x / L) * a_px * this->PI / L + rho_x * cos(a_rhox * this->PI * x / L) * pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L), Scalar(2)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L)) * a_rhoy * this->PI / L + Scalar(2) * u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L)) * a_ux * this->PI / L - u_y * sin(a_uy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L)) * a_uy * this->PI / L + v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L)) * a_vy * this->PI / L;
  return Q_u;
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_q_v(Scalar x,Scalar y)
{
  Scalar Q_v;
  Q_v = p_y * cos(a_py * this->PI * y / L) * a_py * this->PI / L + rho_x * cos(a_rhox * this->PI * x / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L)) * a_rhox * this->PI / L - sin(a_rhoy * this->PI * y / L) * rho_y * pow(v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L), Scalar(2)) * a_rhoy * this->PI / L + cos(a_ux * this->PI * x / L) * u_x * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L)) * a_ux * this->PI / L - sin(a_vx * this->PI * x / L) * v_x * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L)) * a_vx * this->PI / L + Scalar(2) * cos(a_vy * this->PI * y / L) * v_y * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L)) * a_vy * this->PI / L;
  return Q_v;
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_q_e(Scalar x,Scalar y)
{
  Scalar Q_e;
  Q_e = -Gamma * (u_x * sin(x * a_ux * this->PI / L) + u_y * cos(a_uy * this->PI * y / L) + u_0) * a_px * p_x * this->PI * sin(a_px * this->PI * x / L) / (Gamma - Scalar(1)) / L + Gamma * (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(y * a_vy * this->PI / L) + v_0) * a_py * p_y * this->PI * cos(a_py * this->PI * y / L) / (Gamma - Scalar(1)) / L + a_rhox * this->PI * rho_x * cos(a_rhox * this->PI * x / L) * (u_x * sin(x * a_ux * this->PI / L) + u_y * cos(a_uy * this->PI * y / L) + u_0) * (pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(y * a_vy * this->PI / L) + v_0, Scalar(2)) + pow(u_x * sin(x * a_ux * this->PI / L) + u_y * cos(a_uy * this->PI * y / L) + u_0, Scalar(2))) / L / Scalar(2) - a_rhoy * this->PI * rho_y * sin(a_rhoy * this->PI * y / L) * (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(y * a_vy * this->PI / L) + v_0) * (pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(y * a_vy * this->PI / L) + v_0, Scalar(2)) + pow(u_x * sin(x * a_ux * this->PI / L) + u_y * cos(a_uy * this->PI * y / L) + u_0, Scalar(2))) / L / Scalar(2) + (p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_0) * a_ux * this->PI * u_x * cos(x * a_ux * this->PI / L) * Gamma / (Gamma - Scalar(1)) / L + (pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(y * a_vy * this->PI / L) + v_0, Scalar(2)) + Scalar(3) * pow(u_x * sin(x * a_ux * this->PI / L) + u_y * cos(a_uy * this->PI * y / L) + u_0, Scalar(2))) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0) * a_ux * this->PI * u_x * cos(x * a_ux * this->PI / L) / L / Scalar(2) + (p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_0) * a_vy * this->PI * v_y * cos(y * a_vy * this->PI / L) * Gamma / (Gamma - Scalar(1)) / L + (Scalar(3) * pow(v_x * cos(a_vx * this->PI * x / L) + v_y * sin(y * a_vy * this->PI / L) + v_0, Scalar(2)) + pow(u_x * sin(x * a_ux * this->PI / L) + u_y * cos(a_uy * this->PI * y / L) + u_0, Scalar(2))) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0) * a_vy * this->PI * v_y * cos(y * a_vy * this->PI / L) / L / Scalar(2) - (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(y * a_vy * this->PI / L) + v_0) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0) * (u_x * sin(x * a_ux * this->PI / L) + u_y * cos(a_uy * this->PI * y / L) + u_0) * this->PI * a_uy * u_y * sin(a_uy * this->PI * y / L) / L - (v_x * cos(a_vx * this->PI * x / L) + v_y * sin(y * a_vy * this->PI / L) + v_0) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_0) * (u_x * sin(x * a_ux * this->PI / L) + u_y * cos(a_uy * this->PI * y / L) + u_0) * this->PI * a_vx * v_x * sin(a_vx * this->PI * x / L) / L;
  return(Q_e);
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_q_rho(Scalar x,Scalar y)
{
  Scalar Q_rho;
  Q_rho = rho_x * cos(a_rhox * this->PI * x / L) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L)) * a_rhoy * this->PI / L + u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * a_ux * this->PI / L + v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L)) * a_vy * this->PI / L;
  return(Q_rho);
}

// ----------------------------------------
//   Gradient of Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_2d_g_u(Scalar x,Scalar y, int i)
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
Scalar MASA::euler_2d<Scalar>::eval_2d_g_v(Scalar x,Scalar y, int i)
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
      grad =  v_y * cos(a_vy * this->pi * y / L) * a_vy * this->pi / L;
      break;

    default:
      std::cout << "MASA error:: eval_2d_g_v has no " << i << "th component\n";
      std::cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_2d_g_p(Scalar x,Scalar y, int i)
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
Scalar MASA::euler_2d<Scalar>::eval_2d_g_rho(Scalar x,Scalar y, int i)
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
Scalar MASA::euler_2d<Scalar>::eval_an_u(Scalar x,Scalar y)
{
  Scalar u_an;
  u_an = u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L); 
  return u_an;
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_an_v(Scalar x,Scalar y)
{
  Scalar v_an;
  v_an = v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L);
  return v_an;
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_an_p(Scalar x,Scalar y)
{
  Scalar p_an;
  p_an = p_0 + p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L);
  return p_an;
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_an_rho(Scalar x,Scalar y)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L); 
  return rho_an;
}

/* ------------------------------------------------
 *
 *         EULER EQUATION 3d
 *
 *
 *
 * -----------------------------------------------
 */ 


template <typename Scalar>
MASA::euler_3d<Scalar>::euler_3d()
{
  this->mmsname = "euler_3d";
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
int MASA::euler_3d<Scalar>::init_var()
{
  int err = 0;

  // set params (random currenly)
  err += this->set_var("R",1.01);
  err += this->set_var("k",1.38);

  err += this->set_var("u_0",2.27);
  err += this->set_var("u_x",6.00);
  err += this->set_var("u_y",5.35);
  err += this->set_var("u_z",5.34);
  err += this->set_var("v_0",3.46);
  err += this->set_var("v_x",6.13);
  err += this->set_var("v_y",.54);
  err += this->set_var("v_z",.30);
  err += this->set_var("w_0",.411);
  err += this->set_var("w_x",3.14);
  err += this->set_var("w_y",5.68);
  err += this->set_var("w_z",6.51);
  err += this->set_var("rho_0",1.63);
  err += this->set_var("rho_x",4.7);
  err += this->set_var("rho_y",20.85);
  err += this->set_var("rho_z",12.15);
  err += this->set_var("p_0",50.135);
  err += this->set_var("p_x",.73);
  err += this->set_var("p_y",49);
  err += this->set_var("p_z",60.8);
  err += this->set_var("a_px",388.8);
  err += this->set_var("a_py",40.1);
  err += this->set_var("a_pz",38.5);
  err += this->set_var("a_rhox",.82);
  err += this->set_var("a_rhoy",.41);
  err += this->set_var("a_rhoz",.44);
  err += this->set_var("a_ux",.46);
  err += this->set_var("a_uy",.425);
  err += this->set_var("a_uz",.42);
  err += this->set_var("a_vx",.52);
  err += this->set_var("a_vy",.23);
  err += this->set_var("a_vz",16.2);
  err += this->set_var("a_wx",11.05);
  err += this->set_var("a_wy",21.8);
  err += this->set_var("a_wz",13.6);
  err += this->set_var("Gamma",27.5);
  err += this->set_var("mu",12.01);
  err += this->set_var("L",3.02);

  return err;

} // done with variable initializer

// ----------------------------------------
//   Gradient of Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_3d_g_u(Scalar x,Scalar y,Scalar z,int i)
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
Scalar MASA::euler_3d<Scalar>::eval_3d_g_v(Scalar x,Scalar y,Scalar z,int i)
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
Scalar MASA::euler_3d<Scalar>::eval_3d_g_w(Scalar x,Scalar y,Scalar z,int i)
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
Scalar MASA::euler_3d<Scalar>::eval_3d_g_p(Scalar x,Scalar y,Scalar z,int i)
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
      grad =  p_y * cos(a_py * this->pi * y / L) * a_py * this->pi / L;
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
Scalar MASA::euler_3d<Scalar>::eval_3d_g_rho(Scalar x,Scalar y,Scalar z,int i)
{

  Scalar grad = -1;

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
//   Source Term
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_q_u(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_u;
  Q_u = -p_x * sin(a_px * this->PI * x / L) * a_px * this->PI / L + rho_x * cos(a_rhox * this->PI * x / L) * pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_rhoy * this->PI / L + rho_z * cos(a_rhoz * this->PI * z / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_rhoz * this->PI / L + Scalar(2) * u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_ux * this->PI / L - u_y * sin(a_uy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_uy * this->PI / L - u_z * sin(a_uz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_uz * this->PI / L + v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_vy * this->PI / L - w_z * sin(a_wz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_wz * this->PI / L;
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_q_v(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_v;
  Q_v = p_y * cos(a_py * this->PI * y / L) * a_py * this->PI / L + rho_x * cos(a_rhox * this->PI * x / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * pow(v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L), Scalar(2)) * a_rhoy * this->PI / L + rho_z * cos(a_rhoz * this->PI * z / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_rhoz * this->PI / L + u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_ux * this->PI / L - v_x * sin(a_vx * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_vx * this->PI / L + Scalar(2) * v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_vy * this->PI / L + v_z * cos(a_vz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_vz * this->PI / L - w_z * sin(a_wz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_wz * this->PI / L;
  return(Q_v);
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_q_w(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_w;
  Q_w = -p_z * sin(a_pz * this->PI * z / L) * a_pz * this->PI / L + rho_x * cos(a_rhox * this->PI * x / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_rhoy * this->PI / L + rho_z * cos(a_rhoz * this->PI * z / L) * pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2)) * a_rhoz * this->PI / L + u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_ux * this->PI / L + v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_vy * this->PI / L + w_x * cos(a_wx * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_wx * this->PI / L + w_y * cos(a_wy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_wy * this->PI / L - Scalar(2) * w_z * sin(a_wz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_wz * this->PI / L;
  return(Q_w);
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_q_e(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_e;
  Q_e = -Gamma * p_x * sin(a_px * this->PI * x / L) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_px * this->PI / L / (Gamma - Scalar(1)) + Gamma * p_y * cos(a_py * this->PI * y / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_py * this->PI / L / (Gamma - Scalar(1)) - Gamma * p_z * sin(a_pz * this->PI * z / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_pz * this->PI / L / (Gamma - Scalar(1)) + (pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) + pow(v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L), Scalar(2)) + pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2))) * a_rhox * this->PI * rho_x * cos(a_rhox * this->PI * x / L) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) / L / Scalar(2) - (pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) + pow(v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L), Scalar(2)) + pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2))) * a_rhoy * this->PI * rho_y * sin(a_rhoy * this->PI * y / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) / L / Scalar(2) + (pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) + pow(v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L), Scalar(2)) + pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2))) * a_rhoz * this->PI * rho_z * cos(a_rhoz * this->PI * z / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) / L / Scalar(2) + a_ux * this->PI * u_x * cos(a_ux * this->PI * x / L) * ((Scalar(3) * rho_x * sin(a_rhox * this->PI * x / L) + Scalar(3) * rho_y * cos(a_rhoy * this->PI * y / L) + Scalar(3) * rho_z * sin(a_rhoz * this->PI * z / L) + Scalar(3) * rho_0) * pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) + (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * pow(v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L), Scalar(2)) + (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2))) / L / Scalar(2) + Gamma * a_ux * this->PI * u_x * cos(a_ux * this->PI * x / L) * (p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_z * cos(a_pz * this->PI * z / L) + p_0) / L / (Gamma - Scalar(1)) - (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * this->PI * a_uy * u_y * sin(a_uy * this->PI * y / L) / L - (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * this->PI * a_uz * u_z * sin(a_uz * this->PI * z / L) / L - (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * this->PI * a_vx * v_x * sin(a_vx * this->PI * x / L) / L + a_vy * this->PI * v_y * cos(a_vy * this->PI * y / L) * ((rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) + (Scalar(3) * rho_x * sin(a_rhox * this->PI * x / L) + Scalar(3) * rho_y * cos(a_rhoy * this->PI * y / L) + Scalar(3) * rho_z * sin(a_rhoz * this->PI * z / L) + Scalar(3) * rho_0) * pow(v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L), Scalar(2)) + (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2))) / L / Scalar(2) + Gamma * a_vy * this->PI * v_y * cos(a_vy * this->PI * y / L) * (p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_z * cos(a_pz * this->PI * z / L) + p_0) / L / (Gamma - Scalar(1)) + (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * this->PI * a_vz * v_z * cos(a_vz * this->PI * z / L) / L + (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * this->PI * a_wx * w_x * cos(a_wx * this->PI * x / L) / L + (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * this->PI * a_wy * w_y * cos(a_wy * this->PI * y / L) / L - a_wz * this->PI * w_z * sin(a_wz * this->PI * z / L) * ((rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * pow(u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L), Scalar(2)) + (rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L) + rho_0) * pow(v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L), Scalar(2)) + (Scalar(3) * rho_x * sin(a_rhox * this->PI * x / L) + Scalar(3) * rho_y * cos(a_rhoy * this->PI * y / L) + Scalar(3) * rho_z * sin(a_rhoz * this->PI * z / L) + Scalar(3) * rho_0) * pow(w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L), Scalar(2))) / L / Scalar(2) - Gamma * a_wz * this->PI * w_z * sin(a_wz * this->PI * z / L) * (p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_z * cos(a_pz * this->PI * z / L) + p_0) / L / (Gamma - Scalar(1));
  return(Q_e);
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_q_rho(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_rho;
  Q_rho = rho_x * cos(a_rhox * this->PI * x / L) * (u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L)) * a_rhox * this->PI / L - rho_y * sin(a_rhoy * this->PI * y / L) * (v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L)) * a_rhoy * this->PI / L + rho_z * cos(a_rhoz * this->PI * z / L) * (w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L)) * a_rhoz * this->PI / L + u_x * cos(a_ux * this->PI * x / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * a_ux * this->PI / L + v_y * cos(a_vy * this->PI * y / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * a_vy * this->PI / L - w_z * sin(a_wz * this->PI * z / L) * (rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L)) * a_wz * this->PI / L;
  return(Q_rho);
}

// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_an_u(Scalar x,Scalar y,Scalar z)
{
  Scalar u_an;
  u_an = u_0 + u_x * sin(a_ux * this->PI * x / L) + u_y * cos(a_uy * this->PI * y / L) + u_z * cos(a_uz * this->PI * z / L);
  return u_an;
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_an_v(Scalar x,Scalar y,Scalar z)
{
  Scalar v_an;
  v_an = v_0 + v_x * cos(a_vx * this->PI * x / L) + v_y * sin(a_vy * this->PI * y / L) + v_z * sin(a_vz * this->PI * z / L);
  return v_an;
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_an_w(Scalar x,Scalar y,Scalar z)
{
  Scalar w_an;
  w_an = w_0 + w_x * sin(a_wx * this->PI * x / L) + w_y * sin(a_wy * this->PI * y / L) + w_z * cos(a_wz * this->PI * z / L);  
  return w_an;
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_an_p(Scalar x,Scalar y,Scalar z)
{
  Scalar p_an;
  p_an = p_0 + p_x * cos(a_px * this->PI * x / L) + p_y * sin(a_py * this->PI * y / L) + p_z * cos(a_pz * this->PI * z / L);
  return p_an;
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_an_rho(Scalar x,Scalar y,Scalar z)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * this->PI * x / L) + rho_y * cos(a_rhoy * this->PI * y / L) + rho_z * sin(a_rhoz * this->PI * z / L);
  return rho_an;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::euler_1d);
MASA_INSTANTIATE_ALL(MASA::euler_2d);
MASA_INSTANTIATE_ALL(MASA::euler_3d);
