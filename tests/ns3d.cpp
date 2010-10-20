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
// ns3d.cpp: program that tests navier-stokes-3d against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h> // for MASA_STRICT_REGRESSION
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace MASA;

typedef double Scalar;

const Scalar pi = acos(-1);
const Scalar threshold = 1.0e-15; // should be small enough to catch any obvious problems

Scalar nancheck(Scalar x)
{
  if(isnan(x))
    {
      cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
  return 1;
}

Scalar anQ_p (Scalar x,Scalar y,Scalar z,Scalar p_0,Scalar p_x,Scalar p_y,Scalar p_z,Scalar a_px,Scalar a_py,Scalar a_pz,Scalar L)
{
  Scalar p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L);
  return p_an;
}
 
Scalar anQ_u (Scalar x,Scalar y,Scalar z,Scalar u_0,Scalar u_x,Scalar u_y,Scalar u_z,Scalar a_ux,Scalar a_uy,Scalar a_uz,Scalar L)
{
  Scalar u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L);  
  return u_an;
} 
 
Scalar anQ_v (Scalar x,Scalar y,Scalar z,Scalar v_0,Scalar v_x,Scalar v_y,Scalar v_z,Scalar a_vx,Scalar a_vy,Scalar a_vz,Scalar L)
{
  Scalar v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L);
  return v_an;
}

Scalar anQ_w (Scalar x,Scalar y,Scalar z,Scalar w_0,Scalar w_x,Scalar w_y,Scalar w_z,Scalar a_wx,Scalar a_wy,Scalar a_wz,Scalar L)
{
  Scalar w_an = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L);
  return w_an;
}

Scalar anQ_rho (Scalar x,Scalar y,Scalar z,Scalar rho_0,Scalar rho_x,Scalar rho_y,Scalar rho_z,Scalar a_rhox,Scalar a_rhoy,Scalar a_rhoz,Scalar L)
{ 
  Scalar rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L);
  return rho_an;
}

Scalar SourceQ_e (
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar);


Scalar SourceQ_u (
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar);


Scalar SourceQ_v (
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar);

Scalar SourceQ_w (
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar);

Scalar SourceQ_rho(
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar,
  Scalar);

int main()
{
  //variables
  Scalar u_0;
  Scalar u_x;
  Scalar u_y;
  Scalar u_z;
  Scalar v_0;
  Scalar v_x;
  Scalar v_y;
  Scalar v_z;
  Scalar w_0;
  Scalar w_x;
  Scalar w_y;
  Scalar w_z;
  Scalar rho_0;
  Scalar rho_x;
  Scalar rho_y;
  Scalar rho_z;
  Scalar p_0;
  Scalar p_x;
  Scalar p_y;
  Scalar p_z;
  Scalar a_px;
  Scalar a_py;
  Scalar a_pz;
  Scalar a_rhox;
  Scalar a_rhoy;
  Scalar a_rhoz;
  Scalar a_ux;
  Scalar a_uy;
  Scalar a_uz;
  Scalar a_vx;
  Scalar a_vy;
  Scalar a_vz;
  Scalar a_wx;
  Scalar a_wy;
  Scalar a_wz;
  Scalar mu;
  Scalar Gamma;
  Scalar L;    
  Scalar k;

  Scalar R;
  Scalar K;    

  // parameters
  Scalar x;
  Scalar y;
  Scalar z;

  // solutions
  Scalar ufield,ufield2,ufield3;
  Scalar vfield,vfield2,vfield3;
  Scalar wfield,wfield2,wfield3;
  Scalar efield,efield2,efield3;
  Scalar rho,rho2,rho3;
  Scalar gradx,grady,gradz,gradp,gradrho;

  Scalar u_an,u_an2,u_an3;
  Scalar v_an,v_an2,v_an3;
  Scalar w_an,w_an2,w_an3;
  Scalar p_an,p_an2,p_an3;
  Scalar rho_an,rho_an2,rho_an3;

  // initalize
  int nx = 20;             // number of points
  int ny = 10;
  int nz = 30;
  int lx=1;                // length
  int ly=2; 
  int lz=3;     
  Scalar dx=Scalar(lx)/Scalar(nx);
  Scalar dy=Scalar(ly)/Scalar(ny);
  Scalar dz=Scalar(lz)/Scalar(nz);

  masa_init<Scalar>("navier-stokes-test","navierstokes_3d_compressible");

  // set params -- do not call default values for this problem!
  masa_set_param<Scalar>("R",1.01);
  masa_set_param<Scalar>("k",1.38);
  masa_set_param<Scalar>("u_0",2.27);
  masa_set_param<Scalar>("u_x",6.00);
  masa_set_param<Scalar>("u_y",1.35);
  masa_set_param<Scalar>("u_z",5.34);
  masa_set_param<Scalar>("v_0",3.46);
  masa_set_param<Scalar>("v_x",6.13);
  masa_set_param<Scalar>("v_y",.54);
  masa_set_param<Scalar>("v_z",.30);
  masa_set_param<Scalar>("w_0",.411);
  masa_set_param<Scalar>("w_x",3.14);
  masa_set_param<Scalar>("w_y",5.68);
  masa_set_param<Scalar>("w_z",6.51);
  masa_set_param<Scalar>("rho_0",1.63);
  masa_set_param<Scalar>("rho_x",4.7);
  masa_set_param<Scalar>("rho_y",0.85);
  masa_set_param<Scalar>("rho_z",12.15);
  masa_set_param<Scalar>("p_0",1.135);
  masa_set_param<Scalar>("p_x",.73);
  masa_set_param<Scalar>("p_y",4.9);
  masa_set_param<Scalar>("p_z",1.8);
  masa_set_param<Scalar>("a_px",8.8);
  masa_set_param<Scalar>("a_py",0.1);
  masa_set_param<Scalar>("a_pz",38.5);
  masa_set_param<Scalar>("a_rhox",.82);
  masa_set_param<Scalar>("a_rhoy",.41);
  masa_set_param<Scalar>("a_rhoz",.44);
  masa_set_param<Scalar>("a_ux",.46);
  masa_set_param<Scalar>("a_uy",.425);
  masa_set_param<Scalar>("a_uz",.42);
  masa_set_param<Scalar>("a_vx",.52);
  masa_set_param<Scalar>("a_vy",.23);
  masa_set_param<Scalar>("a_vz",1.2);
  masa_set_param<Scalar>("a_wx",1.05);
  masa_set_param<Scalar>("a_wy",1.8);
  masa_set_param<Scalar>("a_wz",13.6);
  masa_set_param<Scalar>("Gamma",2.5);
  masa_set_param<Scalar>("mu",2.01);
  masa_set_param<Scalar>("L",1.0);

  // now set reference values for maple compare
  u_0 = masa_get_param<Scalar>("u_0");
  u_x = masa_get_param<Scalar>("u_x");
  u_y = masa_get_param<Scalar>("u_y");
  u_z = masa_get_param<Scalar>("u_z");

  v_0 = masa_get_param<Scalar>("v_0");
  v_x = masa_get_param<Scalar>("v_x");
  v_y = masa_get_param<Scalar>("v_y");
  v_z = masa_get_param<Scalar>("v_z");

  w_0 = masa_get_param<Scalar>("w_0");
  w_x = masa_get_param<Scalar>("w_x");
  w_y = masa_get_param<Scalar>("w_y");
  w_z = masa_get_param<Scalar>("w_z");

  rho_0 = masa_get_param<Scalar>("rho_0");
  rho_x = masa_get_param<Scalar>("rho_x");
  rho_y = masa_get_param<Scalar>("rho_y");
  rho_z = masa_get_param<Scalar>("rho_z");

  p_0 = masa_get_param<Scalar>("p_0");
  p_x = masa_get_param<Scalar>("p_x");
  p_y = masa_get_param<Scalar>("p_y");
  p_z = masa_get_param<Scalar>("p_z");

  a_px = masa_get_param<Scalar>("a_px");
  a_py = masa_get_param<Scalar>("a_py");
  a_pz = masa_get_param<Scalar>("a_pz");

  a_rhox = masa_get_param<Scalar>("a_rhox");
  a_rhoy = masa_get_param<Scalar>("a_rhoy");
  a_rhoz = masa_get_param<Scalar>("a_rhoz");

  a_ux = masa_get_param<Scalar>("a_ux");
  a_uy = masa_get_param<Scalar>("a_uy");
  a_uz = masa_get_param<Scalar>("a_uz");

  a_vx = masa_get_param<Scalar>("a_vx");
  a_vy = masa_get_param<Scalar>("a_vy");
  a_vz = masa_get_param<Scalar>("a_vz");

  a_wx = masa_get_param<Scalar>("a_wx");
  a_wy = masa_get_param<Scalar>("a_wy");
  a_wz = masa_get_param<Scalar>("a_wz");

  Gamma = masa_get_param<Scalar>("Gamma");
  mu    = masa_get_param<Scalar>("mu");
  L     = masa_get_param<Scalar>("L");

  R = masa_get_param<Scalar>("R");
  K = masa_get_param<Scalar>("k");

  // check all vars are initialized
  int err = masa_sanity_check<Scalar>();
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

	  // evalulate source terms
	  ufield = masa_eval_u_source  <Scalar>(x,y,z);
	  vfield = masa_eval_v_source  <Scalar>(x,y,z);
	  wfield = masa_eval_w_source  <Scalar>(x,y,z);
	  efield = masa_eval_e_source  <Scalar>(x,y,z);
	  rho    = masa_eval_rho_source<Scalar>(x,y,z);
	  
	  // evaluate analytical terms
	  u_an = masa_eval_u_an        <Scalar>(x,y,z);
	  v_an = masa_eval_v_an        <Scalar>(x,y,z);
	  w_an = masa_eval_w_an        <Scalar>(x,y,z);
	  p_an = masa_eval_p_an        <Scalar>(x,y,z);
	  rho_an = masa_eval_rho_an    <Scalar>(x,y,z);

	  // eval gradient terms
	  gradx = masa_eval_3d_grad_u<Scalar>(x,y,z,1);
	  grady = masa_eval_3d_grad_u<Scalar>(x,y,z,2);
	  gradz = masa_eval_3d_grad_u<Scalar>(x,y,z,3);

	  gradx = masa_eval_3d_grad_v<Scalar>(x,y,z,1);
	  grady = masa_eval_3d_grad_v<Scalar>(x,y,z,2);
	  gradz = masa_eval_3d_grad_v<Scalar>(x,y,z,3);
	  
	  gradx = masa_eval_3d_grad_w<Scalar>(x,y,z,1);
	  grady = masa_eval_3d_grad_w<Scalar>(x,y,z,2);
	  gradz = masa_eval_3d_grad_w<Scalar>(x,y,z,3);

	  gradx = masa_eval_3d_grad_p<Scalar>(x,y,z,1);
	  grady = masa_eval_3d_grad_p<Scalar>(x,y,z,2);
	  gradz = masa_eval_3d_grad_p<Scalar>(x,y,z,3);

	  gradx = masa_eval_3d_grad_rho<Scalar>(x,y,z,1);
	  grady = masa_eval_3d_grad_rho<Scalar>(x,y,z,2);
	  gradz = masa_eval_3d_grad_rho<Scalar>(x,y,z,3);

	  // check against maple output
	  ufield2   = SourceQ_u  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,L,R,K);
	  vfield2   = SourceQ_v  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,L,R,K);
	  wfield2   = SourceQ_w  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,L,R,K);
	  rho2      = SourceQ_rho(x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,L,R,K);
	  efield2   = SourceQ_e  (x,y,z,u_0,u_x,u_y,u_z,v_0,v_x,v_y,v_z,w_0,w_x,w_y,w_z,rho_0,rho_x,rho_y,rho_z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,a_rhox,a_rhoy,a_rhoz,a_ux,a_uy,a_uz,a_vx,a_vy,a_vz,a_wx,a_wy,a_wz,mu,Gamma,L,R,K);

	  u_an2     = anQ_u   (x,y,z,u_0,u_x,u_y,u_z,a_ux,a_uy,a_uz,L);
	  v_an2     = anQ_v   (x,y,z,v_0,v_x,v_y,v_z,a_vx,a_vy,a_vz,L);
	  w_an2     = anQ_w   (x,y,z,w_0,w_x,w_y,w_z,a_wx,a_wy,a_wz,L);
	  rho_an2   = anQ_rho (x,y,z,rho_0,rho_x,rho_y,rho_z,a_rhox,a_rhoy,a_rhoz,L);
	  p_an2     = anQ_p   (x,y,z,p_0,p_x,p_y,p_z,a_px,a_py,a_pz,L);

	  // test the result is roughly zero
	  // choose between abs and rel error
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

	  nancheck(ufield3);
	  nancheck(vfield3);
	  nancheck(wfield3);
	  nancheck(efield3);
	  nancheck(rho3);
	  
	  nancheck(u_an3);
	  nancheck(v_an3);
	  nancheck(w_an3);
	  nancheck(rho_an3);
	  nancheck(p_an3);

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

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.
  freopen("/dev/null","w",stdout);

  // test gradient error terms
  Scalar derr = masa_eval_3d_grad_u<Scalar>(0,0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_3d_grad_u<Scalar>(0,0,0,4);
  if(derr != -1)
    {
      cout << "MASA :: gradient (4) error condition failed!\n";
      exit(1);
    }

  // v
  derr = masa_eval_3d_grad_v<Scalar>(0,0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: v gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_3d_grad_v<Scalar>(0,0,0,4);
  if(derr != -1)
    {
      cout << "MASA :: v gradient (4) error condition failed!\n";
      exit(1);
    }

  // w
  derr = masa_eval_3d_grad_w<Scalar>(0,0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: w gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_3d_grad_w<Scalar>(0,0,0,4);
  if(derr != -1)
    {
      cout << "MASA :: w gradient (4) error condition failed!\n";
      exit(1);
    }

  // p
  derr = masa_eval_3d_grad_p<Scalar>(0,0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: p gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_3d_grad_p<Scalar>(0,0,0,4);
  if(derr != -1)
    {
      cout << "MASA :: p gradient (4) error condition failed!\n";
      exit(1);
    }

  // rho
  derr = masa_eval_3d_grad_rho<Scalar>(0,0,0,0);
  if(derr != -1)
    {
      cout << "MASA :: rho gradient (0) error condition failed!\n";
      exit(1);
    }
  
  derr = masa_eval_3d_grad_rho<Scalar>(0,0,0,4);
  if(derr != -1)
    {
      cout << "MASA :: rho gradient (4) error condition failed!\n";
      exit(1);
    }

  // tests passed
  return 0;
}
