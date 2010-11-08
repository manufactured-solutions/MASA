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
// cmasa.cpp: all c-code bindings
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>
#include <string>

using namespace MASA;

extern "C" int cmasa_init(const char* specificname,const char* functionname)
{

  //  printf("sizeofname = %i\n",sizeof(specificname));
  //  printf("functionname = %i\n",sizeof(functionname));
  //  printf(" -> %c%c%c\n",specificname[5],specificname[6],specificname[7]);
  
  std::string sn(specificname);
  std::string fn(functionname);

  //  printf("masa_init, received %s %s\n\n",specificname,functionname);

  masa_init<double>(sn,fn);

  return 0;
}

extern "C" int cmasa_select_mms(const char* function_user_wants)
{
  std::string fuw(function_user_wants);
  masa_select_mms<double>(fuw);
  return 0;
}

extern "C" int cmasa_list_mms()
{
  masa_list_mms<double>();
  return 0;
}

extern "C" int cmasa_init_param()
{
  masa_init_param<double>();
  return 0;
}

extern "C" int cmasa_sanity_check()
{
  masa_sanity_check<double>();
  return 0;
}

extern "C" int cmasa_display_param()
{
  masa_display_param<double>();
  return 0;
}

extern "C" void cmasa_set_param(const char* param,double val)
{
  return(masa_set_param<double>(param,val));
}

extern "C" double cmasa_get_param(const char* param)
{
  return(masa_get_param<double>(param));
}

// --------------------------------
//
//    Source and Analytical Terms
// 
// --------------------------------

// --------------------------------
// source term(s) -- 1D
// --------------------------------

extern "C" double cmasa_eval_1d_t_source  (double x){return(masa_eval_t_source<double>(x));};
extern "C" double cmasa_eval_1d_u_source  (double x){return(masa_eval_u_source<double>(x));};
extern "C" double cmasa_eval_1d_e_source  (double x){return(masa_eval_e_source<double>(x));};
extern "C" double cmasa_eval_1d_rho_source(double x){return(masa_eval_rho_source<double>(x));};

extern "C" double cmasa_eval_1d_t_an      (double x){return(masa_eval_t_an<double>(x));};
extern "C" double cmasa_eval_1d_u_an      (double x){return(masa_eval_u_an<double>(x));};
extern "C" double cmasa_eval_1d_p_an      (double x){return(masa_eval_p_an<double>(x));};
extern "C" double cmasa_eval_1d_rho_an    (double x){return(masa_eval_rho_an<double>(x));};

// --------------------------------
// source term(s) -- 2D
// --------------------------------

extern "C" double cmasa_eval_2d_t_source  (double x,double y){return masa_eval_t_source<double>  (x,y);};
extern "C" double cmasa_eval_2d_u_source  (double x,double y){return(masa_eval_u_source<double>  (x,y));};
extern "C" double cmasa_eval_2d_v_source  (double x,double y){return(masa_eval_v_source<double>  (x,y));};
extern "C" double cmasa_eval_2d_e_source  (double x,double y){return(masa_eval_e_source<double>  (x,y));};
extern "C" double cmasa_eval_2d_rho_source(double x,double y){return(masa_eval_rho_source<double>(x,y));};

extern "C" double cmasa_eval_2d_t_an      (double x,double y){return(masa_eval_t_an<double>  (x,y));};
extern "C" double cmasa_eval_2d_u_an      (double x,double y){return(masa_eval_u_an<double>  (x,y));};
extern "C" double cmasa_eval_2d_v_an      (double x,double y){return(masa_eval_v_an<double>  (x,y));};
extern "C" double cmasa_eval_2d_p_an      (double x,double y){return(masa_eval_p_an<double>  (x,y));};
extern "C" double cmasa_eval_2d_rho_an    (double x,double y){return(masa_eval_rho_an<double>(x,y));};

// --------------------------------
// source term(s) -- 3D
// --------------------------------

extern "C" double cmasa_eval_3d_t_source  (double x,double y,double z){return masa_eval_t_source<double>  (x,y,z); };
extern "C" double cmasa_eval_3d_u_source  (double x,double y,double z){return(masa_eval_u_source<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_v_source  (double x,double y,double z){return(masa_eval_v_source<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_w_source  (double x,double y,double z){return(masa_eval_w_source<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_e_source  (double x,double y,double z){return(masa_eval_e_source<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_rho_source(double x,double y,double z){return(masa_eval_rho_source<double>(x,y,z));};

extern "C" double cmasa_eval_3d_t_an      (double x,double y,double z){return(masa_eval_t_an<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_u_an      (double x,double y,double z){return(masa_eval_u_an<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_v_an      (double x,double y,double z){return(masa_eval_v_an<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_w_an      (double x,double y,double z){return(masa_eval_w_an<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_p_an      (double x,double y,double z){return(masa_eval_p_an<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_rho_an    (double x,double y,double z){return(masa_eval_rho_an<double>(x,y,z));};
