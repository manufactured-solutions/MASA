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
// $Author: nick $
// $Id: ns2d.cpp 12693 2010-08-26 03:35:34Z nick $
//
// axi_cns.cpp: program that tests axisymmetric navier-stokes
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <masa.h>
#include <math.h> 
using namespace std;
using namespace MASA;

const double pi = acos(-1);
const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double anQ_p (double x,double y,double p_0,double p_x,double p_y,double a_px,double a_py,double L)
{
  double p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L);
  return p_an;
}
  
double anQ_u (double x,double y,double u_0,double u_x,double u_y,double a_ux,double a_uy,double L)
{
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  return u_an;
} 
 
double anQ_v (double x,double y,double v_0,double v_x,double v_y,double a_vx,double a_vy,double L)
{
  double v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
  return v_an;
}

double anQ_rho (double x,double y,double rho_0,double rho_x,double rho_y,double a_rhox,double a_rhoy,double L)
{ 
  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L);
  return rho_an;
}

double SourceQ_e (
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double);

double SourceQ_u (
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double);

double SourceQ_v (
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double);

double SourceQ_rho( 
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double,
  double);

int main()
{
  //variables
  double u_0;
  double u_x;
  double u_y;
  double v_0;
  double v_x;
  double v_y;
  double rho_0;
  double rho_x;
  double rho_y;
  double p_0;
  double p_x;
  double p_y;
  double a_px;
  double a_py;
  double a_rhox;
  double a_rhoy;
  double a_ux;
  double a_uy;
  double a_vx;
  double a_vy;
  double Gamma;
  double mu;
  double L;
  double R;
  double k;

  // parameters
  double x;
  double y;
  double z;

  // solutions
  double ufield,ufield2;
  double vfield,vfield2;
  double efield,efield2;
  double rho,rho2;

  // initalize
  int nx = 10;  // number of points
  int ny = 8;  
  int lx=2;     // length
  int ly=1; 
  
  double dx=double(lx)/double(nx);
  double dy=double(ly)/double(ny);

  masa_init("navier-stokes-test","navierstokes_2d_compressible");

  // set params
  masa_init_param();
  
  // get vars
  masa_get_param("u_0",&u_0);
  masa_get_param("u_x",&u_x);
  masa_get_param("u_y",&u_y);
  masa_get_param("v_0",&v_0);
  masa_get_param("v_x",&v_x);
  masa_get_param("v_y",&v_y);

  masa_get_param("rho_0",&rho_0);
  masa_get_param("rho_x",&rho_x);
  masa_get_param("rho_y",&rho_y);

  masa_get_param("p_0",&p_0);
  masa_get_param("p_x",&p_x);
  masa_get_param("p_y",&p_y);

  masa_get_param("a_px",&a_px);
  masa_get_param("a_py",&a_py);

  masa_get_param("a_rhox",&a_rhox);
  masa_get_param("a_rhoy",&a_rhoy);

  masa_get_param("a_ux",&a_ux);
  masa_get_param("a_uy",&a_uy);

  masa_get_param("a_vx",&a_vx);
  masa_get_param("a_vy",&a_vy);

  masa_get_param("Gamma",&Gamma);
  masa_get_param("mu",&mu);
  masa_get_param("L",&L);

  masa_get_param("R",&R);
  masa_get_param("k",&k);

  // check that all terms have been initialized
  masa_sanity_check();

  // evaluate source terms (2D)
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)    
      {
	x=i*dx;
	y=j*dy;

	masa_eval_u_source  (x,y,&ufield);
	masa_eval_v_source  (x,y,&vfield);
	masa_eval_e_source  (x,y,&efield);
	masa_eval_rho_source(x,y,&rho);

	ufield2   = SourceQ_u  (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);
	vfield2   = SourceQ_v  (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);
	rho2      = SourceQ_rho(x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,mu,L,R,k);
	efield2   = SourceQ_e  (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,Gamma,mu,L,R,k);

	// test the result is roughly zero
	ufield = fabs(ufield-ufield2);
	vfield = fabs(vfield-vfield2);
	efield = fabs(efield-efield2);
	rho    = fabs(rho-rho2);
  
	//cout << endl << ufield << endl << vfield << endl << efield << rho << endl;

	if(ufield > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout << "U Field Source Term\n";
	    cout << "Exceeded Threshold by: " << ufield << endl;
	    exit(1);
	  }

	if(vfield > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout << "V Field Source Term\n";
	    cout << "Exceeded Threshold by: " << vfield << endl;
	    exit(1);
	  }

	if(efield > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout << "Energy Source Term\n";
	    cout << "Exceeded Threshold by: " << efield << endl;
	    exit(1);
	  }

	if(rho > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Navier-Stokes\n";
	    cout << "RHO Source Term\n";
	    cout << "Exceeded Threshold by: " << rho << endl;
	    exit(1);
	  }
      } // done iterating

  // tests passed
  return 0;
}
