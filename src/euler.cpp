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

  err += this->set_var("u_0",14.191);
  err += this->set_var("u_x",1.63);
  err += this->set_var("rho_0",91.5);
  err += this->set_var("rho_x",.13);
  err += this->set_var("p_0",12.1984);
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
  Scalar Q_u;
  Scalar RHO;
  Scalar U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);

  Q_u = 0.2e1 * cos(a_ux * PI * x / L) * RHO * U * a_ux * PI * u_x / L + cos(a_rhox * PI * x / L) * U * U * a_rhox * PI * rho_x / L - sin(a_px * PI * x / L) * a_px * PI * p_x / L;
  
  return(Q_u);
}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_q_e(Scalar x)
{
  Scalar Q_e;

  Scalar RHO;
  Scalar U;
  Scalar P;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  P = p_0 + p_x * cos(a_px * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  Q_e = cos(a_rhox * PI * x / L) * pow(U, 0.3e1) * a_rhox * PI * rho_x / L / 0.2e1 + cos(a_ux * PI * x / L) * P * a_ux * PI * u_x * Gamma / L / (Gamma - 0.1e1) + 0.3e1 / 0.2e1 * cos(a_ux * PI * x / L) * RHO * U * U * a_ux * PI * u_x / L - sin(a_px * PI * x / L) * U * a_px * PI * p_x * Gamma / L / (Gamma - 0.1e1);

  return(Q_e);
}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_q_rho(Scalar x)
{
  Scalar Q_rho;
  Scalar RHO;
  Scalar U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);

  Q_rho = cos(a_ux * PI * x / L) * RHO * a_ux * PI * u_x / L + cos(a_rhox * PI * x / L) * U * a_rhox * PI * rho_x / L;
  return(Q_rho);


  return(Q_rho);
}

// ----------------------------------------
//   Gradient of Source Terms
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_1d_g_u(Scalar x)
{

  Scalar grad_u = u_x * cos(a_ux * pi * x / L) * a_ux * pi / L;
  return grad_u;


}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_1d_g_p(Scalar x)
{
  Scalar grad_p = -p_x * sin(a_px * pi * x / L) * a_px * pi / L;
  return grad_p;

}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_1d_g_rho(Scalar x)
{
  Scalar grad_rho = rho_x * cos(a_rhox * pi * x / L) * a_rhox * pi / L;
  return grad_rho;

}


// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_an_u(Scalar x)
{
  Scalar u_an;
  u_an = u_0 + u_x * sin(a_ux * pi * x / L);
  return u_an;
}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_an_p(Scalar x)
{
  Scalar p_an;
  p_an = p_0 + p_x * cos(a_px * pi * x / L);
  return p_an;
}

template <typename Scalar>
Scalar MASA::euler_1d<Scalar>::eval_an_rho(Scalar x)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L);
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

  err += this->set_var("u_0",200.23);
  err += this->set_var("u_x",1.1);
  err += this->set_var("u_y",1.08);
  err += this->set_var("v_0",1.2);
  err += this->set_var("v_x",1.6);
  err += this->set_var("v_y",.47);
  err += this->set_var("rho_0",100.02);
  err += this->set_var("rho_x",2.22);
  err += this->set_var("rho_y",0.8);
  err += this->set_var("p_0",150.2);
  err += this->set_var("p_x",.91);
  err += this->set_var("p_y",.623);
  err += this->set_var("a_px",.165);
  err += this->set_var("a_py",.612);
  err += this->set_var("a_rhox",1.0);
  err += this->set_var("a_rhoy",1.0);
  err += this->set_var("a_ux",.1987);
  err += this->set_var("a_uy",1.189);
  err += this->set_var("a_vx",1.91);
  err += this->set_var("a_vy",1.0);
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
  Scalar RHO;
  Scalar U;
  Scalar V;

  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);

  Q_u = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * sin(a_rhoy * PI * y / L) / L - a_uy * PI * u_y * RHO * V * sin(a_uy * PI * y / L) / L - a_px * PI * p_x * sin(a_px * PI * x / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U / L;

  return Q_u;
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_q_v(Scalar x,Scalar y)
{
  Scalar Q_v;
  Scalar RHO;
  Scalar U;
  Scalar V;
  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);

  Q_v = a_rhox * PI * rho_x * U * V * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * V * sin(a_rhoy * PI * y / L) / L - a_vx * PI * v_x * RHO * U * sin(a_vx * PI * x / L) / L + a_py * PI * p_y * cos(a_py * PI * y / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V / L;

  return Q_v;
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_q_e(Scalar x,Scalar y)
{
  Scalar Q_e;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar P;

  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L);

  Q_e = -a_px * PI * p_x * Gamma * U * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + a_py * PI * p_y * Gamma * V * cos(a_py * PI * y / L) / (Gamma - 0.1e1) / L + (U * U + V * V) * a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L / 0.2e1 - (U * U + V * V) * a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L / 0.2e1 + (0.3e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * U * U / L / 0.2e1 - (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * PI * RHO * U * V / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.3e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO * V * V / L / 0.2e1 + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * Gamma * P / (Gamma - 0.1e1) / L;

  return(Q_e);
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_q_rho(Scalar x,Scalar y)
{
  Scalar Q_rho;
  Scalar RHO;
  Scalar U;
  Scalar V;

  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);

  Q_rho = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * RHO / L;

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
      grad = u_x * cos(a_ux * pi * x / L) * a_ux * pi / L;      
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
      grad = -v_x * sin(a_vx * pi * x / L) * a_vx * pi / L;
      break;

    case 2:
      grad =  v_y * cos(a_vy * pi * y / L) * a_vy * pi / L;
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
Scalar MASA::euler_2d<Scalar>::eval_an_u(Scalar x,Scalar y)
{
  Scalar u_an;
  u_an = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L); 
  return u_an;
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_an_v(Scalar x,Scalar y)
{
  Scalar v_an;
  v_an = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);
  return v_an;
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_an_p(Scalar x,Scalar y)
{
  Scalar p_an;
  p_an = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L);
  return p_an;
}

template <typename Scalar>
Scalar MASA::euler_2d<Scalar>::eval_an_rho(Scalar x,Scalar y)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L); 
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

  err += this->set_var("u_0",20.27);
  err += this->set_var("u_x",6.00);
  err += this->set_var("u_y",5.35);
  err += this->set_var("u_z",5.34);
  err += this->set_var("v_0",30.46);
  err += this->set_var("v_x",6.13);
  err += this->set_var("v_y",.54);
  err += this->set_var("v_z",.30);
  err += this->set_var("w_0",20.411);
  err += this->set_var("w_x",3.14);
  err += this->set_var("w_y",5.68);
  err += this->set_var("w_z",6.51);
  err += this->set_var("rho_0",10.63);
  err += this->set_var("rho_x",0.3);
  err += this->set_var("rho_y",0.15);
  err += this->set_var("rho_z",1.15);
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
      grad =  u_x * cos(a_ux * pi * x / L) * a_ux * pi / L;
      break;

    case 2:
      grad = -u_y * sin(a_uy * pi * y / L) * a_uy * pi / L;
      break;

    case 3:
      grad = -u_z * sin(a_uz * pi * z / L) * a_uz * pi / L;
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
      grad = -v_x * sin(a_vx * pi * x / L) * a_vx * pi / L;
      break;

    case 2:
      grad =  v_y * cos(a_vy * pi * y / L) * a_vy * pi / L;
      break;

    case 3:
      grad =  v_z * cos(a_vz * pi * z / L) * a_vz * pi / L;
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
      grad =  w_x * cos(a_wx * pi * x / L) * a_wx * pi / L;
      break;

    case 2:
      grad =  w_y * cos(a_wy * pi * y / L) * a_wy * pi / L;
      break;

    case 3:
      grad = -w_z * sin(a_wz * pi * z / L) * a_wz * pi / L;
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
      grad = -p_x * sin(a_px * pi * x / L) * a_px * pi / L;
      break;

    case 2:
      grad =  p_y * cos(a_py * pi * y / L) * a_py * pi / L;
      break;

    case 3:
      grad = -p_z * sin(a_pz * pi * z / L) * a_pz * pi / L;
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
      grad = rho_x * cos(a_rhox * pi * x / L) * a_rhox * pi / L;
      break;
      
    case 2:
      grad = -rho_y * sin(a_rhoy * pi * y / L) * a_rhoy * pi / L;
      break;
      
    case 3:
      grad = rho_z * cos(a_rhoz * pi * z / L) * a_rhoz * pi / L;
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
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L);
  W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L);

  Q_u = a_rhox * PI * rho_x * U * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * U * V * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * U * W * cos(a_rhoz * PI * z / L) / L - a_uy * PI * u_y * RHO * V * sin(a_uy * PI * y / L) / L - a_uz * PI * u_z * RHO * W * sin(a_uz * PI * z / L) / L - a_px * PI * p_x * sin(a_px * PI * x / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * U / L;

  return(Q_u);
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_q_v(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_v;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L);
  W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L);

  Q_v = a_rhox * PI * rho_x * U * V * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * V * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * V * W * cos(a_rhoz * PI * z / L) / L - a_vx * PI * v_x * RHO * U * sin(a_vx * PI * x / L) / L + a_vz * PI * v_z * RHO * W * cos(a_vz * PI * z / L) / L + a_py * PI * p_y * cos(a_py * PI * y / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * V / L;

  return(Q_v);
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_q_w(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_w;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L);
  W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L);

  Q_w = a_rhox * PI * rho_x * U * W * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * W * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * W * W * cos(a_rhoz * PI * z / L) / L + a_wx * PI * w_x * RHO * U * cos(a_wx * PI * x / L) / L + a_wy * PI * w_y * RHO * V * cos(a_wy * PI * y / L) / L - a_pz * PI * p_z * sin(a_pz * PI * z / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - 0.2e1 * a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * W / L;

  return(Q_w);
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_q_e(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_e;
  Scalar RHO;
  Scalar P;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L);
  W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L);
  P = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L);

  Q_e = -a_px * PI * p_x * Gamma * U * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + a_py * PI * p_y * Gamma * V * cos(a_py * PI * y / L) / (Gamma - 0.1e1) / L - a_pz * PI * p_z * Gamma * W * sin(a_pz * PI * z / L) / (Gamma - 0.1e1) / L + (U * U + V * V + W * W) * a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L / 0.2e1 - (U * U + V * V + W * W) * a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L / 0.2e1 + (U * U + V * V + W * W) * a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L / 0.2e1 - (-0.3e1 * a_ux * u_x * cos(a_ux * PI * x / L) - a_vy * v_y * cos(a_vy * PI * y / L) + a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * U * U / L / 0.2e1 - (a_uy * u_y * sin(a_uy * PI * y / L) + a_vx * v_x * sin(a_vx * PI * x / L)) * PI * RHO * U * V / L - (a_uz * u_z * sin(a_uz * PI * z / L) - a_wx * w_x * cos(a_wx * PI * x / L)) * PI * RHO * U * W / L - (-a_ux * u_x * cos(a_ux * PI * x / L) - 0.3e1 * a_vy * v_y * cos(a_vy * PI * y / L) + a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * V * V / L / 0.2e1 + (a_vz * v_z * cos(a_vz * PI * z / L) + a_wy * w_y * cos(a_wy * PI * y / L)) * PI * RHO * V * W / L - (-a_ux * u_x * cos(a_ux * PI * x / L) - a_vy * v_y * cos(a_vy * PI * y / L) + 0.3e1 * a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO * W * W / L / 0.2e1 - (-a_ux * u_x * cos(a_ux * PI * x / L) - a_vy * v_y * cos(a_vy * PI * y / L) + a_wz * w_z * sin(a_wz * PI * z / L)) * PI * Gamma * P / (Gamma - 0.1e1) / L;

  return(Q_e);
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_q_rho(Scalar x,Scalar y,Scalar z)
{
  Scalar Q_rho;
  Scalar RHO;
  Scalar U;
  Scalar V;
  Scalar W;

  RHO = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L);
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L);
  W = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L);

  Q_rho = a_rhox * PI * rho_x * U * cos(a_rhox * PI * x / L) / L - a_rhoy * PI * rho_y * V * sin(a_rhoy * PI * y / L) / L + a_rhoz * PI * rho_z * W * cos(a_rhoz * PI * z / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L) - a_wz * w_z * sin(a_wz * PI * z / L)) * PI * RHO / L;

  return(Q_rho);
}

// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_an_u(Scalar x,Scalar y,Scalar z)
{
  Scalar u_an;
  u_an = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L);
  return u_an;
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_an_v(Scalar x,Scalar y,Scalar z)
{
  Scalar v_an;
  v_an = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L);
  return v_an;
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_an_w(Scalar x,Scalar y,Scalar z)
{
  Scalar w_an;
  w_an = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L);  
  return w_an;
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_an_p(Scalar x,Scalar y,Scalar z)
{
  Scalar p_an;
  p_an = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L);
  return p_an;
}

template <typename Scalar>
Scalar MASA::euler_3d<Scalar>::eval_an_rho(Scalar x,Scalar y,Scalar z)
{
  Scalar rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L);
  return rho_an;
}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::euler_1d);
MASA_INSTANTIATE_ALL(MASA::euler_2d);
MASA_INSTANTIATE_ALL(MASA::euler_3d);
