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
#include <math.h>

#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace MASA;

const double pi = acos(-1);
const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

int main()
{
  //variables
  double u_0;
  double u_x;
  double u_y;
  double u_z;
  double v_0;
  double v_x;
  double v_y;
  double v_z;
  double w_0;
  double w_x;
  double w_y;
  double w_z;
  double rho_0;
  double rho_x;
  double rho_y;
  double rho_z;
  double p_0;
  double p_x;
  double p_y;
  double p_z;
  double a_px;
  double a_py;
  double a_pz;
  double a_rhox;
  double a_rhoy;
  double a_rhoz;
  double a_ux;
  double a_uy;
  double a_uz;
  double a_vx;
  double a_vy;
  double a_vz;
  double a_wx;
  double a_wy;
  double a_wz;
  double mu;
  double Gamma;
  double L;    
  double k;

  double R;
  double K;    

  // parameters
  double x;
  double y;
  double z;

  // solutions
  double ufield,ufield2,ufield3;
  double vfield,vfield2,vfield3;
  double wfield,wfield2,wfield3;
  double efield,efield2,efield3;
  double rho,rho2,rho3;

  double u_an,u_an2,u_an3;
  double v_an,v_an2,v_an3;
  double w_an,w_an2,w_an3;
  double p_an,p_an2,p_an3;
  double rho_an,rho_an2,rho_an3;

  // initalize
  int nx = 20;             // number of points
  int ny = 10;
  int nz = 30;
  int lx=1;                // length
  int ly=2; 
  int lz=3;     
  double dx=double(lx)/double(nx);
  double dy=double(ly)/double(ny);
  double dz=double(lz)/double(nz);

  masa_init("ns3d","navierstokes_3d_compressible");

  // set params
  masa_init_params();

  // now set reference values for maple compare
  u_0 = masa_get_param("u_0");
  u_x = masa_get_param("u_x");
  u_y = masa_get_param("u_y");

  v_0 = masa_get_param("v_0");
  v_x = masa_get_param("v_x");
  v_y = masa_get_param("v_y");

  rho_0 = masa_get_param("rho_0");
  rho_x = masa_get_param("rho_x");
  rho_y = masa_get_param("rho_y");

  p_0 = masa_get_param("p_0");
  p_x = masa_get_param("p_x");
  p_y = masa_get_param("p_y");

  a_px = masa_get_param("a_px");
  a_py = masa_get_param("a_py");

  a_rhox = masa_get_param("a_rhox");
  a_rhoy = masa_get_param("a_rhoy");

  a_ux = masa_get_param("a_ux");
  a_uy = masa_get_param("a_uy");

  a_vx = masa_get_param("a_vx");
  a_vy = masa_get_param("a_vy");

  Gamma = masa_get_param("Gamma");
  mu    = masa_get_param("mu");
  L     = masa_get_param("L");

  R = masa_get_param("R");
  K = masa_get_param("k");

  // check all vars are initialized
  int err = masa_sanity_check();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  // now start up 2d case
  masa_init("ns2d","navierstokes_2d_compressible");

  // set as reference values
  masa_set_param("u_0",u_0);
  masa_set_param("u_x",u_x);
  masa_set_param("u_y",u_y);

  masa_set_param("v_0",v_0);
  masa_set_param("v_x",v_x);
  masa_set_param("v_y",v_y);

  masa_set_param("rho_0",rho_0);
  masa_set_param("rho_x",rho_x);
  masa_set_param("rho_y",rho_y);

  masa_set_param("p_0",p_0);
  masa_set_param("p_x",p_x);
  masa_set_param("p_y",p_y);

  masa_set_param("a_px",a_px);
  masa_set_param("a_py",a_py);

  masa_set_param("a_rhox",a_rhox);
  masa_set_param("a_rhoy",a_rhoy);

  masa_set_param("a_ux",a_ux);
  masa_set_param("a_uy",a_uy);

  masa_set_param("a_vx",a_vx);
  masa_set_param("a_vy",a_vy);

  masa_set_param("Gamma",Gamma);
  masa_set_param("mu",mu);
  masa_set_param("L",L);

  masa_set_param("R",R);
  masa_set_param("k",k);

  // check all vars are initialized
  int err = masa_sanity_check();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }
  
  // evaluate source terms (3D)
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)    
      for(int k=0;k<nz;k++)
	{
	  x=i*dx;
	  y=j*dy;
	  z=k*dz;

	  // switch to 3D
	  masa_select_mms("ns3d");

	  // evalulate source terms
	  ufield = masa_eval_u_source  (x,y,z);
	  vfield = masa_eval_v_source  (x,y,z);
	  wfield = masa_eval_w_source  (x,y,z);
	  efield = masa_eval_e_source  (x,y,z);
	  rho    = masa_eval_rho_source(x,y,z);
	  
	  // evaluate analytical terms
	  u_an = masa_eval_u_an        (x,y,z);
	  v_an = masa_eval_v_an        (x,y,z);
	  w_an = masa_eval_w_an        (x,y,z);
	  p_an = masa_eval_p_an        (x,y,z);
	  rho_an = masa_eval_rho_an    (x,y,z);

	  // now switch to 2D and do the same
	  masa_select_mms("ns2d");

	  // evalulate source terms
	  ufield2 = masa_eval_u_source  (x,y);
	  vfield2 = masa_eval_v_source  (x,y);
	  wfield2 = masa_eval_w_source  (x,y);
	  efield2 = masa_eval_e_source  (x,y);
	  rho2    = masa_eval_rho_source(x,y);

	  // evaluate analytical terms
	  u_an2   = masa_eval_u_an      (x,y);
	  v_an2   = masa_eval_v_an      (x,y);
	  w_an2   = masa_eval_w_an      (x,y);
	  p_an2   = masa_eval_p_an      (x,y);
	  rho_an2 = masa_eval_rho_an    (x,y);

#ifdef MASA_STRICT_REGRESSION

	  ufield3 = fabs(ufield-ufield2);
	  vfield3 = fabs(vfield-vfield2);
	  wfield3 = fabs(wfield-wfield2);
	  efield3 = fabs(efield-efield2);
	  rho3    = fabs(rho-rho2);

	  u_an3   = fabs(u_an-u_an2);
	  v_an3   = fabs(v_an-v_an2);
	  w_an3   = fabs(w_an-w_an2);
	  rho_an3 = fabs(rho_an-rho_an2);
	  p_an3   = fabs(p_an-p_an2);

#else
	  ufield3 = fabs(ufield-ufield2)/fabs(ufield2);
	  vfield3 = fabs(vfield-vfield2)/fabs(vfield2);
	  wfield3 = fabs(wfield-wfield2)/fabs(wfield2);
	  efield3 = fabs(efield-efield2)/fabs(efield2);
	  rho3    = fabs(rho-rho2)/fabs(rho2);

	  u_an3   = fabs(u_an-u_an2)/fabs(u_an2);
	  v_an3   = fabs(v_an-v_an2)/fabs(v_an2);
	  w_an3   = fabs(w_an-w_an2)/fabs(w_an2);
	  rho_an3 = fabs(rho_an-rho_an2)/fabs(rho_an2);
	  p_an3   = fabs(p_an-p_an2)/fabs(p_an2);
#endif

	  if(ufield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "U Field Source Term\n";
	      cout << "Exceeded Threshold by: " << ufield << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(u_an3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "U Field Analytical Term\n";
	      cout << "Exceeded Threshold by: " << u_an << endl;
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
	  
	  if(v_an3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "V Field Analytical Term\n";
	      cout << "Exceeded Threshold by: " << v_an << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(wfield3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "W Field Source Term\n";
	      cout << "Exceeded Threshold by: " << wfield << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }
	  
	  if(w_an3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "W Field Analytical Term\n";
	      cout << "Exceeded Threshold by: " << w_an3 << endl;
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

	  if(p_an3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "P Field Analytical Term\n";
	      cout << "Exceeded Threshold by: " << p_an << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(rho3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "RHO Source Term\n";
	      cout << "Exceeded Threshold by: " << rho << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	  if(rho_an3 > threshold)
	    {
	      cout << "\nMASA REGRESSION TEST FAILED: Navier-Stokes 3d\n";
	      cout << "RHO Analytical Term\n";
	      cout << "Exceeded Threshold by: " << rho_an << endl;
	      cout << x << " " << y << " " << z << endl;
	      exit(1);
	    }

	} // done iterating

  // tests passed
  return 0;
}
