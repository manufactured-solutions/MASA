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
// ns-2d-3d.cpp: program that tests navier-stokes-3d against ns2d, when w=0
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h> // for MASA_STRICT_REGRESSION
#include <masa.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace MASA;

typedef double Real;

const Real pi = acos(-1);

// with strict, between 7-8 decimals
const Real threshold = 1.0e-7; // should be small enough to catch any obvious problems

int main()
{
  //variables
  Real u_0;
  Real u_x;
  Real u_y;
  Real u_z;
  Real v_0;
  Real v_x;
  Real v_y;
  Real v_z;
  Real w_0;
  Real w_x;
  Real w_y;
  Real w_z;
  Real rho_0;
  Real rho_x;
  Real rho_y;
  Real rho_z;
  Real p_0;
  Real p_x;
  Real p_y;
  Real p_z;
  Real a_px;
  Real a_py;
  Real a_pz;
  Real a_rhox;
  Real a_rhoy;
  Real a_rhoz;
  Real a_ux;
  Real a_uy;
  Real a_uz;
  Real a_vx;
  Real a_vy;
  Real a_vz;
  Real a_wx;
  Real a_wy;
  Real a_wz;
  Real mu;
  Real Gamma;
  Real L;    
  Real k;

  Real R;
  Real K;    

  // parameters
  Real x;
  Real y;
  Real z;

  // solutions
  Real ufield,ufield2,ufield3;
  Real vfield,vfield2,vfield3;
  Real wfield,wfield2,wfield3;
  Real efield,efield2,efield3;
  Real rho,rho2,rho3;

  Real u_exact,u_exact2,u_exact3;
  Real v_exact,v_exact2,v_exact3;
  Real w_exact,w_exact2,w_exact3;
  Real p_exact,p_exact2,p_exact3;
  Real rho_exact,rho_exact2,rho_exact3;

  // initalize
  int nx = 20;             // number of points
  int ny = 10;
  int nz = 30;
  int lx=1;                // length
  int ly=2; 
  int lz=3;     
  Real dx=Real(lx)/Real(nx);
  Real dy=Real(ly)/Real(ny);
  Real dz=Real(lz)/Real(nz);

  masa_init<Real>("ns3d","navierstokes_3d_compressible");

  // set params
  masa_init_param<Real>();

  // now set reference values for maple compare
  u_0 = masa_get_param<Real>("u_0");
  u_x = masa_get_param<Real>("u_x");
  u_y = masa_get_param<Real>("u_y");

  v_0 = masa_get_param<Real>("v_0");
  v_x = masa_get_param<Real>("v_x");
  v_y = masa_get_param<Real>("v_y");

  rho_0 = masa_get_param<Real>("rho_0");
  rho_x = masa_get_param<Real>("rho_x");
  rho_y = masa_get_param<Real>("rho_y");

  p_0 = masa_get_param<Real>("p_0");
  p_x = masa_get_param<Real>("p_x");
  p_y = masa_get_param<Real>("p_y");

  a_px = masa_get_param<Real>("a_px");
  a_py = masa_get_param<Real>("a_py");

  a_rhox = masa_get_param<Real>("a_rhox");
  a_rhoy = masa_get_param<Real>("a_rhoy");

  a_ux = masa_get_param<Real>("a_ux");
  a_uy = masa_get_param<Real>("a_uy");

  a_vx = masa_get_param<Real>("a_vx");
  a_vy = masa_get_param<Real>("a_vy");

  Gamma = masa_get_param<Real>("Gamma");
  mu    = masa_get_param<Real>("mu");
  L     = masa_get_param<Real>("L");

  R = masa_get_param<Real>("R");
  K = masa_get_param<Real>("k");

  // check all vars are initialized
  int err = masa_sanity_check<Real>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  // now start up 2d case
  masa_init<Real>("ns2d","navierstokes_2d_compressible");

  // set as reference values
  masa_set_param<Real>("u_0",u_0);
  masa_set_param<Real>("u_x",u_x);
  masa_set_param<Real>("u_y",u_y);

  masa_set_param<Real>("v_0",v_0);
  masa_set_param<Real>("v_x",v_x);
  masa_set_param<Real>("v_y",v_y);

  masa_set_param<Real>("rho_0",rho_0);
  masa_set_param<Real>("rho_x",rho_x);
  masa_set_param<Real>("rho_y",rho_y);

  masa_set_param<Real>("p_0",p_0);
  masa_set_param<Real>("p_x",p_x);
  masa_set_param<Real>("p_y",p_y);

  masa_set_param<Real>("a_px",a_px);
  masa_set_param<Real>("a_py",a_py);

  masa_set_param<Real>("a_rhox",a_rhox);
  masa_set_param<Real>("a_rhoy",a_rhoy);

  masa_set_param<Real>("a_ux",a_ux);
  masa_set_param<Real>("a_uy",a_uy);

  masa_set_param<Real>("a_vx",a_vx);
  masa_set_param<Real>("a_vy",a_vy);

  masa_set_param<Real>("Gamma",Gamma);
  masa_set_param<Real>("mu",mu);
  masa_set_param<Real>("L",L);

  masa_set_param<Real>("R",R);
  masa_set_param<Real>("k",K);

  // check all vars are initialized
  err = masa_sanity_check<Real>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  //masa_display_param<Real>();
  //masa_select_mms<Real>("ns3d");
  //masa_display_param<Real>();
  //exit(1);

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa.
  freopen("/dev/null","w",stdout);

  // evaluate source terms (3D)
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)    
      for(int k=0;k<nz;k++)
	{
	  x=i*dx;
	  y=j*dy;
	  z=k*dz;

	  // switch to 3D
	  masa_select_mms<Real>("ns3d");

	  // evalulate source terms
	  ufield = masa_eval_u_source  <Real>(x,y,z);
	  vfield = masa_eval_v_source  <Real>(x,y,z);
	  efield = masa_eval_e_source  <Real>(x,y,z);
	  rho    = masa_eval_rho_source<Real>(x,y,z);
	  
	  // evaluate analytical terms
	  u_exact = masa_eval_u_exact        <Real>(x,y,z);
	  v_exact = masa_eval_v_exact        <Real>(x,y,z);
	  p_exact = masa_eval_p_exact        <Real>(x,y,z);
	  rho_exact = masa_eval_rho_exact    <Real>(x,y,z);

	  // now switch to 2D and do the same
	  masa_select_mms<Real>("ns2d");

	  // evalulate source terms
	  ufield2 = masa_eval_u_source  <Real>(x,y);
	  vfield2 = masa_eval_v_source  <Real>(x,y);
	  efield2 = masa_eval_e_source  <Real>(x,y);
	  rho2    = masa_eval_rho_source<Real>(x,y);

	  // evaluate analytical terms
	  u_exact2   = masa_eval_u_exact      <Real>(x,y);
	  v_exact2   = masa_eval_v_exact      <Real>(x,y);
	  p_exact2   = masa_eval_p_exact      <Real>(x,y);
	  rho_exact2 = masa_eval_rho_exact    <Real>(x,y);

#ifdef MASA_STRICT_REGRESSION

	  ufield3 = std::abs(ufield-ufield2);
	  vfield3 = std::abs(vfield-vfield2);
	  efield3 = std::abs(efield-efield2);
	  rho3    = std::abs(rho-rho2);

	  u_exact3   = std::abs(u_exact-u_exact2);
	  v_exact3   = std::abs(v_exact-v_exact2);
	  rho_exact3 = std::abs(rho_exact-rho_exact2);
	  p_exact3   = std::abs(p_exact-p_exact2);

#else
	  ufield3 = std::abs(ufield-ufield2)/std::abs(ufield2);
	  vfield3 = std::abs(vfield-vfield2)/std::abs(vfield2);
	  efield3 = std::abs(efield-efield2)/std::abs(efield2);
	  rho3    = std::abs(rho-rho2)/std::abs(rho2);

	  u_exact3   = std::abs(u_exact-u_exact2)/std::abs(u_exact2);
	  v_exact3   = std::abs(v_exact-v_exact2)/std::abs(v_exact2);
	  rho_exact3 = std::abs(rho_exact-rho_exact2)/std::abs(rho_exact2);
	  p_exact3   = std::abs(p_exact-p_exact2)/std::abs(p_exact2);
#endif

	  if(ufield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "U Field Source Term\n";
	      cout << "Exceeded Threshold by: " << ufield << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(u_exact3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "U Field Analytical Term\n";
	      cout << "Exceeded Threshold by: " << u_exact << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(vfield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "V Field Source Term\n";
	      cout << "Exceeded Threshold by: " << vfield << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }
	  
	  if(v_exact3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "V Field Analytical Term\n";
	      cout << "Exceeded Threshold by: " << v_exact << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(efield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "Energy Source Term\n";
	      cout.precision(16);
	      cout << "Exceeded Threshold by: " << efield3 << endl;
	      cout << "Source term is:                   " << efield2 << endl;
	      cout << "MASA term is:                     " << efield << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(p_exact3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "P Field Analytical Term\n";
	      cout << "Exceeded Threshold by: " << p_exact << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(rho3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "RHO Source Term\n";
	      cout << "Exceeded Threshold by: " << rho3 << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(rho_exact3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "RHO Analytical Term\n";
	      cout << "Exceeded Threshold by: " << rho_exact << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	} // done iterating

  // tests passed
  return 0;
}
