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

MASA::euler_1d::euler_1d()
{
  mmsname = "euler_1d";
  dimension=1;

  register_var("R",&R);
  register_var("k",&k);

  register_var("u_0",&u_0);
  register_var("u_x",&u_x);
  register_var("rho_0",&rho_0);
  register_var("rho_x",&rho_x);
  register_var("p_0",&p_0);
  register_var("p_x",&p_x);
  register_var("a_px",&a_px);
  register_var("a_rhox",&a_rhox);
  register_var("a_ux",&a_ux);
  register_var("L",&L);
  register_var("Gamma",&Gamma);
  register_var("mu",&mu);

}//done with constructor

int MASA::euler_1d::init_var()
{
  int err = 0;

  // randomly generated
  err += set_var("R",1.01);
  err += set_var("k",1.38);

  err += set_var("u_0",.191);
  err += set_var("u_x",1.63);
  err += set_var("rho_0",91.5);
  err += set_var("rho_x",5.13);
  err += set_var("p_0",.1984);
  err += set_var("p_x",3.151);
  err += set_var("a_px",6.151);
  err += set_var("a_rhox",1.2);
  err += set_var("a_ux",.03);
  err += set_var("L",3.02);
  err += set_var("Gamma",16.1);
  err += set_var("mu",.091);

  return err;

} // done with variable initializer


// ----------------------------------------
//   Source Terms
// ----------------------------------------

double MASA::euler_1d::eval_q_u(double x)
{
  double Q_u;
  Q_u = -sin(a_px * PI * x / L) * a_px * PI * p_x / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.2e1) * a_rhox * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L)) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_ux * PI / L;
  return Q_u;
}

double MASA::euler_1d::eval_q_e(double x)
{
  double Q_e;
  Q_e = cos(a_rhox * PI * x / L) * rho_x * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.3e1) * a_rhox * PI / L / 0.2e1 + cos(a_ux * PI * x / L) * (p_0 + p_x * cos(a_px * PI * x / L)) * a_ux * PI * u_x * Gamma / L / (Gamma - 0.1e1) - Gamma * p_x * sin(a_px * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_px * PI / L / (Gamma - 0.1e1) + 0.3e1 / 0.2e1 * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L)) * pow(u_0 + u_x * sin(a_ux * PI * x / L), 0.2e1) * a_ux * PI * u_x / L;
  return(Q_e);
}

double MASA::euler_1d::eval_q_rho(double x)
{
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L)) * a_rhox * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L)) * a_ux * PI / L;
  return(Q_rho);
}

// ----------------------------------------
//   Gradient of Source Terms
// ----------------------------------------

double MASA::euler_1d::eval_1d_g(double x)
{

  return 0;


}

// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

double MASA::euler_1d::eval_an_u(double x)
{
  double u_an;
  u_an = u_0 + u_x * sin(a_ux * PI * x / L);
  return u_an;
}

double MASA::euler_1d::eval_an_p(double x)
{
  double p_an;
  p_an = p_0 + p_x * cos(a_px * PI * x / L);
  return p_an;
}

double MASA::euler_1d::eval_an_rho(double x)
{
  double rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * PI * x / L);                                                                                                                                                                                           
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

MASA::euler_2d::euler_2d()
{
  mmsname = "euler_2d";
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

int MASA::euler_2d::init_var()
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

double MASA::euler_2d::eval_q_u(double x,double y)
{
  double Q_u;
  Q_u = -p_x * sin(a_px * PI * x / L) * a_px * PI / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L), 0.2e1) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhoy * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_ux * PI / L - u_y * sin(a_uy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_uy * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vy * PI / L;
  return Q_u;
}

double MASA::euler_2d::eval_q_v(double x,double y)
{
  double Q_v;
  Q_v = p_y * cos(a_py * PI * y / L) * a_py * PI / L + rho_x * cos(a_rhox * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhox * PI / L - sin(a_rhoy * PI * y / L) * rho_y * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L), 0.2e1) * a_rhoy * PI / L + cos(a_ux * PI * x / L) * u_x * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_ux * PI / L - sin(a_vx * PI * x / L) * v_x * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_vx * PI / L + 0.2e1 * cos(a_vy * PI * y / L) * v_y * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_vy * PI / L;
  return Q_v;
}

double MASA::euler_2d::eval_q_e(double x,double y)
{
  double Q_e;
  Q_e = -Gamma * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * a_px * p_x * PI * sin(a_px * PI * x / L) / (Gamma - 0.1e1) / L + Gamma * (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * a_py * p_y * PI * cos(a_py * PI * y / L) / (Gamma - 0.1e1) / L + a_rhox * PI * rho_x * cos(a_rhox * PI * x / L) * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * (pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) / L / 0.2e1 - a_rhoy * PI * rho_y * sin(a_rhoy * PI * y / L) * (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * (pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) / L / 0.2e1 + (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * a_ux * PI * u_x * cos(x * a_ux * PI / L) * Gamma / (Gamma - 0.1e1) / L + (pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + 0.3e1 * pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_ux * PI * u_x * cos(x * a_ux * PI / L) / L / 0.2e1 + (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_0) * a_vy * PI * v_y * cos(y * a_vy * PI / L) * Gamma / (Gamma - 0.1e1) / L + (0.3e1 * pow(v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0, 0.2e1) + pow(u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0, 0.2e1)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * a_vy * PI * v_y * cos(y * a_vy * PI / L) / L / 0.2e1 - (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * PI * a_uy * u_y * sin(a_uy * PI * y / L) / L - (v_x * cos(a_vx * PI * x / L) + v_y * sin(y * a_vy * PI / L) + v_0) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_0) * (u_x * sin(x * a_ux * PI / L) + u_y * cos(a_uy * PI * y / L) + u_0) * PI * a_vx * v_x * sin(a_vx * PI * x / L) / L;
  return(Q_e);
}

double MASA::euler_2d::eval_q_rho(double x,double y)
{
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L)) * a_rhoy * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * a_ux * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L)) * a_vy * PI / L;
  return(Q_rho);
}

// ----------------------------------------
//   Gradient of Source Terms
// ----------------------------------------

double MASA::euler_2d::eval_2d_g(double x,double y, int i)
{

  double grad = -1;

  switch(i)
    {
    case 0:
      cout << "MASA error:: eval_2d_g has no 0th component\n";
      cout << "Try 1 or 2\n";
      break;
      
    case 1:
      
      break;

    case 2:

      break;

    default:
      cout << "MASA error:: eval_2d_g has no " << i << "th component\n";
      cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

double MASA::euler_2d::eval_an_u(double x,double y)
{
  double u_an;
  u_an = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L); 
  return u_an;
}

double MASA::euler_2d::eval_an_v(double x,double y)
{
  double v_an;
  v_an = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);
  return v_an;
}

double MASA::euler_2d::eval_an_p(double x,double y)
{
  double p_an;
  p_an = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L);
  return p_an;
}

double MASA::euler_2d::eval_an_rho(double x,double y)
{
  double rho_an;
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


MASA::euler_3d::euler_3d()
{
  mmsname = "euler_3d";
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

int MASA::euler_3d::init_var()
{
  int err = 0;

  // set params (random currenly)
  err += set_var("R",1.01);
  err += set_var("k",1.38);

  err += set_var("u_0",2.27);
  err += set_var("u_x",6.00);
  err += set_var("u_y",5.35);
  err += set_var("u_z",5.34);
  err += set_var("v_0",3.46);
  err += set_var("v_x",6.13);
  err += set_var("v_y",.54);
  err += set_var("v_z",.30);
  err += set_var("w_0",.411);
  err += set_var("w_x",3.14);
  err += set_var("w_y",5.68);
  err += set_var("w_z",6.51);
  err += set_var("rho_0",1.63);
  err += set_var("rho_x",4.7);
  err += set_var("rho_y",20.85);
  err += set_var("rho_z",12.15);
  err += set_var("p_0",50.135);
  err += set_var("p_x",.73);
  err += set_var("p_y",49);
  err += set_var("p_z",60.8);
  err += set_var("a_px",388.8);
  err += set_var("a_py",40.1);
  err += set_var("a_pz",38.5);
  err += set_var("a_rhox",.82);
  err += set_var("a_rhoy",.41);
  err += set_var("a_rhoz",.44);
  err += set_var("a_ux",.46);
  err += set_var("a_uy",.425);
  err += set_var("a_uz",.42);
  err += set_var("a_vx",.52);
  err += set_var("a_vy",.23);
  err += set_var("a_vz",16.2);
  err += set_var("a_wx",11.05);
  err += set_var("a_wy",21.8);
  err += set_var("a_wz",13.6);
  err += set_var("Gamma",27.5);
  err += set_var("mu",12.01);
  err += set_var("L",3.02);

  return err;

} // done with variable initializer

// ----------------------------------------
//   Gradient of Source Terms
// ----------------------------------------

double MASA::euler_3d::eval_3d_g(double x,double y,double z,int i)
{

  double grad = -1;

  switch(i)
    {
    case 0:
      cout << "MASA error:: masa_eval_3d_grad has no 0th component\n";
      cout << "Try 1,2 or 3\n";
      break;
      
    case 1:
      
      break;

    case 2:

      break;

    case 3:

      break;

    default:
      cout << "MASA error:: masa_eval_3d_grad has no " << i << "th component\n";
      cout << "Try 1 or 2\n";
      break;

    }// done with switch
  
  return grad;

}

// ----------------------------------------
//   Source Term
// ----------------------------------------

double MASA::euler_3d::eval_q_u(double x,double y,double z)
{
  double Q_u;
  Q_u = -p_x * sin(a_px * PI * x / L) * a_px * PI / L + rho_x * cos(a_rhox * PI * x / L) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhoz * PI / L + 0.2e1 * u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_ux * PI / L - u_y * sin(a_uy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_uy * PI / L - u_z * sin(a_uz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_uz * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_vy * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_wz * PI / L;
  return(Q_u);
}

double MASA::euler_3d::eval_q_v(double x,double y,double z)
{
  double Q_v;
  Q_v = p_y * cos(a_py * PI * y / L) * a_py * PI / L + rho_x * cos(a_rhox * PI * x / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_ux * PI / L - v_x * sin(a_vx * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_vx * PI / L + 0.2e1 * v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_vy * PI / L + v_z * cos(a_vz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_vz * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_wz * PI / L;
  return(Q_v);
}

double MASA::euler_3d::eval_q_w(double x,double y,double z)
{
  double Q_w;
  Q_w = -p_z * sin(a_pz * PI * z / L) * a_pz * PI / L + rho_x * cos(a_rhox * PI * x / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_ux * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_vy * PI / L + w_x * cos(a_wx * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_wx * PI / L + w_y * cos(a_wy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_wy * PI / L - 0.2e1 * w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_wz * PI / L;
  return(Q_w);
}

double MASA::euler_3d::eval_q_e(double x,double y,double z)
{
  double Q_e;
  Q_e = -Gamma * p_x * sin(a_px * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_px * PI / L / (Gamma - 0.1e1) + Gamma * p_y * cos(a_py * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_py * PI / L / (Gamma - 0.1e1) - Gamma * p_z * sin(a_pz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_pz * PI / L / (Gamma - 0.1e1) + (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) * a_rhox * PI * rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) / L / 0.2e1 - (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) * a_rhoy * PI * rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) / L / 0.2e1 + (pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) * a_rhoz * PI * rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) / L / 0.2e1 + a_ux * PI * u_x * cos(a_ux * PI * x / L) * ((0.3e1 * rho_x * sin(a_rhox * PI * x / L) + 0.3e1 * rho_y * cos(a_rhoy * PI * y / L) + 0.3e1 * rho_z * sin(a_rhoz * PI * z / L) + 0.3e1 * rho_0) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) / L / 0.2e1 + Gamma * a_ux * PI * u_x * cos(a_ux * PI * x / L) * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L) + p_0) / L / (Gamma - 0.1e1) - (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * PI * a_uy * u_y * sin(a_uy * PI * y / L) / L - (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * PI * a_uz * u_z * sin(a_uz * PI * z / L) / L - (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * PI * a_vx * v_x * sin(a_vx * PI * x / L) / L + a_vy * PI * v_y * cos(a_vy * PI * y / L) * ((rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + (0.3e1 * rho_x * sin(a_rhox * PI * x / L) + 0.3e1 * rho_y * cos(a_rhoy * PI * y / L) + 0.3e1 * rho_z * sin(a_rhoz * PI * z / L) + 0.3e1 * rho_0) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) / L / 0.2e1 + Gamma * a_vy * PI * v_y * cos(a_vy * PI * y / L) * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L) + p_0) / L / (Gamma - 0.1e1) + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * PI * a_vz * v_z * cos(a_vz * PI * z / L) / L + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * PI * a_wx * w_x * cos(a_wx * PI * x / L) / L + (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * PI * a_wy * w_y * cos(a_wy * PI * y / L) / L - a_wz * PI * w_z * sin(a_wz * PI * z / L) * ((rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L), 0.2e1) + (rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L) + rho_0) * pow(v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L), 0.2e1) + (0.3e1 * rho_x * sin(a_rhox * PI * x / L) + 0.3e1 * rho_y * cos(a_rhoy * PI * y / L) + 0.3e1 * rho_z * sin(a_rhoz * PI * z / L) + 0.3e1 * rho_0) * pow(w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L), 0.2e1)) / L / 0.2e1 - Gamma * a_wz * PI * w_z * sin(a_wz * PI * z / L) * (p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L) + p_0) / L / (Gamma - 0.1e1);
  return(Q_e);
}

double MASA::euler_3d::eval_q_rho(double x,double y,double z)
{
  double Q_rho;
  Q_rho = rho_x * cos(a_rhox * PI * x / L) * (u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L)) * a_rhox * PI / L - rho_y * sin(a_rhoy * PI * y / L) * (v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L)) * a_rhoy * PI / L + rho_z * cos(a_rhoz * PI * z / L) * (w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L)) * a_rhoz * PI / L + u_x * cos(a_ux * PI * x / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_ux * PI / L + v_y * cos(a_vy * PI * y / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_vy * PI / L - w_z * sin(a_wz * PI * z / L) * (rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L)) * a_wz * PI / L;
  return(Q_rho);
}

// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

double MASA::euler_3d::eval_an_u(double x,double y,double z)
{
  double u_an;
  u_an = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_z * cos(a_uz * PI * z / L);
  return u_an;
}

double MASA::euler_3d::eval_an_v(double x,double y,double z)
{
  double v_an;
  v_an = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_z * sin(a_vz * PI * z / L);
  return v_an;
}

double MASA::euler_3d::eval_an_w(double x,double y,double z)
{
  double w_an;
  w_an = w_0 + w_x * sin(a_wx * PI * x / L) + w_y * sin(a_wy * PI * y / L) + w_z * cos(a_wz * PI * z / L);  
  return w_an;
}

double MASA::euler_3d::eval_an_p(double x,double y,double z)
{
  double p_an;
  p_an = p_0 + p_x * cos(a_px * PI * x / L) + p_y * sin(a_py * PI * y / L) + p_z * cos(a_pz * PI * z / L);
  return p_an;
}

double MASA::euler_3d::eval_an_rho(double x,double y,double z)
{
  double rho_an;
  rho_an = rho_0 + rho_x * sin(a_rhox * PI * x / L) + rho_y * cos(a_rhoy * PI * y / L) + rho_z * sin(a_rhoz * PI * z / L);
  return rho_an;
}
