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
// uninit.cpp : program that tests error on uninitialized functions
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

using namespace MASA;
using namespace std;

int main()
{
  double derr;
  int err=0;

  double x=0;
  double y=0;
  double z=0;
  double t=0;
  
  masa_init("masa-test","masa_uninit");

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa.
  freopen("/dev/null","w",stdout);

  // start testing uninitialized functions
  derr = masa_init_param();
  if(derr != 1) 
    err += 1;
  
  // --------------------------------
  // source term(s) -- 1D
  // --------------------------------
  
  derr = masa_eval_t_source(x); 
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_t_source(x,y);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_u_source(x);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_v_source(x);  
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_w_source(x);  
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_e_source(x);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_rho_source(x);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_rho_u_source(x,y);
  if(derr != -1.33) 
    err += 1;
  
  derr = masa_eval_t_an(x); 
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_t_an(x,y);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_u_an(x);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_v_an(x);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_w_an(x);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_p_an(x);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_rho_an(x);
  if(derr != -1.33) 
    err += 1;  

  derr = masa_eval_1d_grad(1,x);
  if(derr != -1.33) 
    err += 1;  
  
  // --------------------------------
  // source term(s) -- 2D
  // --------------------------------
  derr = masa_eval_t_source(x,y,z); 
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_u_source(x,y);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_v_source(x,y);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_w_source(x,y);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_e_source(x,y);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_rho_source(x,y);
  if(derr != -1.33) 
    err += 1;
  
  derr = masa_eval_t_an(x,y,z);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_u_an(x,y);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_v_an(x,y);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_w_an(x,y);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_p_an(x,y);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_rho_an(x,y);
  if(derr != -1.33) 
    err += 1;
  
  derr = masa_eval_2d_grad(1,x,y);
  if(derr != -1.33) 
    err += 1;

  // --------------------------------
  // source term(s) -- 3D
  // --------------------------------

  derr = masa_eval_t_source(x,y,z,t);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_u_source(x,y,z);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_v_source(x,y,z);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_w_source(x,y,z);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_e_source(x,y,z);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_rho_source(x,y,z);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_3d_grad(1,x,y,z);
  if(derr != -1.33) 
    err += 1;
  
  derr = masa_eval_t_an(x,y,z,t);
  if(derr != -1.33) 
  err += 1;

  derr = masa_eval_u_an(x,y,z);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_v_an(x,y,z);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_w_an(x,y,z);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_p_an(x,y,z);
  if(derr != -1.33) 
    err += 1;

  derr = masa_eval_rho_an(x,y,z);
  if(derr != -1.33) 
    err += 1;

  if(err != 0)
    {
      cout << "MASA ERROR: Default output for virtual function failing\n";
      cout << err << " functions failing.\n";
      exit(1);
    }

  return 0; // steady as she goes

}
