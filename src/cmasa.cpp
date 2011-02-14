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

// This header file contains the public functions designed to be exposed in CMASA
// What follows is the masa.h doxygen documentation headers

/*! \file masa.h
\brief CMASA header file containing all public C/C++ API

CMASA.h is a header file that contains all the public objects and member functions for the C interfaces. 
To use, be sure to #include <masa.h> in any C programs calling MASA routines.

Several simple examples for using these functions are provided in the examples section. 
Functions which have an integer return value return "0" upon success.

*/


extern "C" int cmasa_init(const char* specificname,const char* functionname)
{
  
  std::string sn(specificname);
  std::string fn(functionname);

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

extern "C" int cmasa_purge_default_param()
{
  masa_purge_default_param<double>();
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
// source, analytical and gradient term(s) -- 1D
// --------------------------------

extern "C" double cmasa_eval_1d_source_t     (double x){return(masa_eval_source_t     <double>  (x));};
extern "C" double cmasa_eval_1d_source_u     (double x){return(masa_eval_source_u     <double>  (x));};
extern "C" double cmasa_eval_1d_source_e     (double x){return(masa_eval_source_e     <double>  (x));};
extern "C" double cmasa_eval_1d_source_rho   (double x){return(masa_eval_source_rho   <double>(x));};
extern "C" double cmasa_eval_1d_source_rho_u (double x){return(masa_eval_source_rho_u <double>(x));};
extern "C" double cmasa_eval_1d_source_rho_e (double x){return(masa_eval_source_rho_e <double>(x));};
extern "C" double cmasa_eval_1d_source_rho_N (double x){return(masa_eval_source_rho_N <double>(x));};
extern "C" double cmasa_eval_1d_source_rho_N2(double x){return(masa_eval_source_rho_N2<double>(x));};

extern "C" double cmasa_eval_1d_exact_t      (double x){return(masa_eval_exact_t<double>  (x));};
extern "C" double cmasa_eval_1d_exact_u      (double x){return(masa_eval_exact_u<double>  (x));};
extern "C" double cmasa_eval_1d_exact_p      (double x){return(masa_eval_exact_p<double>  (x));};
extern "C" double cmasa_eval_1d_exact_rho    (double x){return(masa_eval_exact_rho<double>(x));};
extern "C" double cmasa_eval_1d_exact_rho_N  (double x){return(masa_eval_exact_rho_N<double>(x));};
extern "C" double cmasa_eval_1d_exact_rho_N2 (double x){return(masa_eval_exact_rho_N2<double>(x));};

extern "C" double cmasa_eval_grad_1d_u    (double x){return(masa_eval_grad_1d_u<double>  (x));};
extern "C" double cmasa_eval_grad_1d_p    (double x){return(masa_eval_grad_1d_p<double>  (x));};
extern "C" double cmasa_eval_grad_1d_rho  (double x){return(masa_eval_grad_1d_rho<double>(x));};

// --------------------------------
// source, analytical and gradient term(s) -- 2D
// --------------------------------

extern "C" double cmasa_eval_2d_source_t    (double x,double y){return masa_eval_source_t<double>  (x,y);};
extern "C" double cmasa_eval_2d_source_u    (double x,double y){return(masa_eval_source_u<double>  (x,y));};
extern "C" double cmasa_eval_2d_source_v    (double x,double y){return(masa_eval_source_v<double>  (x,y));};
extern "C" double cmasa_eval_2d_source_e    (double x,double y){return(masa_eval_source_e<double>  (x,y));};
extern "C" double cmasa_eval_2d_source_rho  (double x,double y){return(masa_eval_source_rho  <double>(x,y));};
extern "C" double cmasa_eval_2d_source_rho_u(double x,double y){return(masa_eval_source_rho_u<double>(x,y));};
extern "C" double cmasa_eval_2d_source_rho_v(double x,double y){return(masa_eval_source_rho_v<double>(x,y));};
extern "C" double cmasa_eval_2d_source_rho_w(double x,double y){return(masa_eval_source_rho_w<double>(x,y));};
extern "C" double cmasa_eval_2d_source_rho_e(double x,double y){return(masa_eval_source_rho_e<double>(x,y));};

extern "C" double cmasa_eval_2d_exact_t     (double x,double y){return(masa_eval_exact_t<double>  (x,y));};
extern "C" double cmasa_eval_2d_exact_u     (double x,double y){return(masa_eval_exact_u<double>  (x,y));};
extern "C" double cmasa_eval_2d_exact_v     (double x,double y){return(masa_eval_exact_v<double>  (x,y));};
extern "C" double cmasa_eval_2d_exact_p     (double x,double y){return(masa_eval_exact_p<double>  (x,y));};
extern "C" double cmasa_eval_2d_exact_rho   (double x,double y){return(masa_eval_exact_rho<double>(x,y));};

extern "C" double cmasa_eval_grad_2d_u   (double x,double y,int i){return(masa_eval_grad_2d_u<double>  (x,y,i));};
extern "C" double cmasa_eval_grad_2d_v   (double x,double y,int i){return(masa_eval_grad_2d_v<double>  (x,y,i));};
extern "C" double cmasa_eval_grad_2d_w   (double x,double y,int i){return(masa_eval_grad_2d_w<double>  (x,y,i));};
extern "C" double cmasa_eval_grad_2d_p   (double x,double y,int i){return(masa_eval_grad_2d_p<double>  (x,y,i));};
extern "C" double cmasa_eval_grad_2d_rho (double x,double y,int i){return(masa_eval_grad_2d_rho<double>(x,y,i));};

// --------------------------------
// source, analytical and gradient term(s) -- 3D
// --------------------------------

extern "C" double cmasa_eval_3d_source_t    (double x,double y,double z){return masa_eval_source_t<double>  (x,y,z); };
extern "C" double cmasa_eval_3d_source_u    (double x,double y,double z){return(masa_eval_source_u<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_source_v    (double x,double y,double z){return(masa_eval_source_v<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_source_w    (double x,double y,double z){return(masa_eval_source_w<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_source_e    (double x,double y,double z){return(masa_eval_source_e<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_source_rho  (double x,double y,double z){return(masa_eval_source_rho  <double>(x,y,z));};
extern "C" double cmasa_eval_3d_source_rho_u(double x,double y,double z){return(masa_eval_source_rho_u<double>(x,y,z));};
extern "C" double cmasa_eval_3d_source_rho_v(double x,double y,double z){return(masa_eval_source_rho_v<double>(x,y,z));};
extern "C" double cmasa_eval_3d_source_rho_w(double x,double y,double z){return(masa_eval_source_rho_w<double>(x,y,z));};
extern "C" double cmasa_eval_3d_source_rho_e(double x,double y,double z){return(masa_eval_source_rho_e<double>(x,y,z));};

extern "C" double cmasa_eval_3d_exact_t     (double x,double y,double z){return(masa_eval_exact_t<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_exact_u     (double x,double y,double z){return(masa_eval_exact_u<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_exact_v     (double x,double y,double z){return(masa_eval_exact_v<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_exact_w     (double x,double y,double z){return(masa_eval_exact_w<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_exact_p     (double x,double y,double z){return(masa_eval_exact_p<double>  (x,y,z));};
extern "C" double cmasa_eval_3d_exact_rho   (double x,double y,double z){return(masa_eval_exact_rho<double>(x,y,z));};

extern "C" double cmasa_eval_grad_3d_u   (double x,double y,double z,int i){return(masa_eval_grad_3d_u  <double>(x,y,z,i));};
extern "C" double cmasa_eval_grad_3d_v   (double x,double y,double z,int i){return(masa_eval_grad_3d_v  <double>(x,y,z,i));};
extern "C" double cmasa_eval_grad_3d_w   (double x,double y,double z,int i){return(masa_eval_grad_3d_w  <double>(x,y,z,i));};
extern "C" double cmasa_eval_grad_3d_p   (double x,double y,double z,int i){return(masa_eval_grad_3d_p  <double>(x,y,z,i));};
extern "C" double cmasa_eval_grad_3d_rho (double x,double y,double z,int i){return(masa_eval_grad_3d_rho<double>(x,y,z,i));};
