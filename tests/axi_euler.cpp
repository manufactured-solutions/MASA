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
// $Id: euler2d.cpp 12693 2010-08-26 03:35:34Z nick $
//
// axi_euler.cpp: program that tests axisymmetric euler
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <masa.h>
#include <math.h>

using namespace MASA;
using namespace std;

const double pi = acos(-1);
const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double anQ_p (double x,double y,double p_0,double p_x,double p_y,double a_px,double a_py,double L)
{
  double p_an = p_0 + p_1 * sin(a_pr * PI * r / L) * cos(a_pz * PI * z / L);
  return p_an;
}
  
double anQ_u (double x,double y,double u_0,double u_x,double u_y,double a_ux,double a_uy,double L)
{
  double u_an = u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L);
  return u_an;
} 
 
double anQ_w (double x,double y,double v_0,double v_x,double v_y,double a_vx,double a_vy,double L)
{
  double w_an = w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L);
  return w_an;
}

double anQ_rho (double x,double y,double rho_0,double rho_x,double rho_y,double a_rhox,double a_rhoy,double L)
{ 
  double rho_an = rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L);
  return rho_an;
}

// ----------------------------------------
//   Source Terms
// ----------------------------------------

double SourceQ_e()
{
  double Q_e;
  Q_e = Gamma * cos(a_pz * PI * z / L) * cos(a_pr * PI * r / L) * p_1 * u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L) * a_pr * PI / L / (Gamma - 0.1e1) - Gamma * sin(a_pz * PI * z / L) * sin(a_pr * PI * r / L) * p_1 * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_pz * PI / L / (Gamma - 0.1e1) - sin(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) * sin(a_uz * PI * z / L) * u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * (pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) + pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_rhor * PI * rho_1 / L / 0.2e1 + cos(a_rhor * PI * r / L) * cos(a_rhoz * PI * z / L) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * (pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) + pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_rhoz * PI * rho_1 / L / 0.2e1 - sin(a_uz * PI * z / L) * sin(a_ur * PI * r / L) * (p_0 + p_1 * sin(a_pr * PI * r / L) * cos(a_pz * PI * z / L)) * a_ur * PI * u_1 * Gamma / L / (Gamma - 0.1e1) - sin(a_uz * PI * z / L) * sin(a_ur * PI * r / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) + 0.3e1 * pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_ur * PI * u_1 / L / 0.2e1 + (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * u_1 * u_1 * cos(a_uz * PI * z / L) * sin(a_uz * PI * z / L) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * a_uz * PI / L - sin(a_uz * PI * z / L) * (cos(a_ur * PI * r / L) - 0.1e1) * u_1 * w_1 * sin(a_wr * PI * r / L) * sin(a_wz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_wr * PI / L + cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * (p_0 + p_1 * sin(a_pr * PI * r / L) * cos(a_pz * PI * z / L)) * a_wz * PI * w_1 * Gamma / L / (Gamma - 0.1e1) + cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (0.3e1 * pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) + pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_wz * PI * w_1 / L / 0.2e1 + sin(a_uz * PI * z / L) * (cos(a_ur * PI * r / L) - 0.1e1) * u_1 * (p_0 + p_1 * sin(a_pr * PI * r / L) * cos(a_pz * PI * z / L)) * Gamma / (Gamma - 0.1e1) / r + sin(a_uz * PI * z / L) * (cos(a_ur * PI * r / L) - 0.1e1) * u_1 * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) + pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) / r / 0.2e1;
  return(Q_e);
}

double SourceQ_u()
{
  double Q_u = p_1 * cos(a_pr * PI * r / L) * cos(a_pz * PI * z / L) * a_pr * PI / L - u_1 * u_1 * pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * rho_1 * sin(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) * a_rhor * PI / L + u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * rho_1 * cos(a_rhor * PI * r / L) * cos(a_rhoz * PI * z / L) * sin(a_uz * PI * z / L) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_rhoz * PI / L - 0.2e1 * u_1 * u_1 * pow(sin(a_uz * PI * z / L), 0.2e1) * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_ur * PI * r / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * a_ur * PI / L + u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * cos(a_uz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_uz * PI / L + u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * w_1 * cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * sin(a_uz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * a_wz * PI / L + u_1 * u_1 * pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) / r;
  return(Q_u);  
}

double SourceQ_w()
{
  double Q_w;
  Q_w = -p_1 * sin(a_pr * PI * r / L) * sin(a_pz * PI * z / L) * a_pz * PI / L - u_1 * sin(a_uz * PI * z / L) * rho_1 * sin(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * (cos(a_ur * PI * r / L) - 0.1e1) * a_rhor * PI / L + pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) * rho_1 * cos(a_rhor * PI * r / L) * cos(a_rhoz * PI * z / L) * a_rhoz * PI / L - u_1 * sin(a_uz * PI * z / L) * sin(a_ur * PI * r / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_ur * PI / L - u_1 * sin(a_uz * PI * z / L) * w_1 * sin(a_wr * PI * r / L) * sin(a_wz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (cos(a_ur * PI * r / L) - 0.1e1) * a_wr * PI / L + (0.2e1 * w_0 + 0.2e1 * w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * w_1 * cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * a_wz * PI / L + u_1 * sin(a_uz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * (cos(a_ur * PI * r / L) - 0.1e1) / r;
  return(Q_w);
}

double SourceQ_rho ()
{
  double Q_rho;
  Q_rho = -(cos(a_ur * PI * r / L) - 0.1e1) * a_rhor * PI * rho_1 * u_1 * sin(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) * sin(a_uz * PI * z / L) / L + (w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L) + w_0) * a_rhoz * PI * rho_1 * cos(a_rhor * PI * r / L) * cos(a_rhoz * PI * z / L) / L - (rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) + rho_0) * a_ur * PI * u_1 * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) / L + (rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) + rho_0) * a_wz * PI * w_1 * cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) / L + (rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) + rho_0) * (cos(a_ur * PI * r / L) - 0.1e1) * u_1 * sin(a_uz * PI * z / L) / r;
  return(Q_rho);
}


int main()

{  //variables
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

  // parameters
  double x;
  double y;

  // solutions -- efield is MASA term, efield2 is maple, efield3 is abs error between them
  double ufield,ufield2,ufield3;
  double vfield,vfield2,vfield3;
  double efield,efield2,efield3;
  double rho,rho2;

  double u_an,u_an2;
  double v_an,v_an2,v_an3;
  double p_an,p_an2;
  double rho_an,rho_an2;

  // initalize
  int nx = 115;  // number of points
  int ny = 68;  
  int lx=3;     // length
  int ly=1; 
  
  double dx=double(lx)/double(nx);
  double dy=double(ly)/double(ny);

  masa_init("euler-test","axisymmetric_euler");

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

  // check that all terms have been initialized
  masa_sanity_check();

  // evaluate source terms (2D)
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)    
      {
	x=i*dx;
	y=j*dy;
	
	//evalulate source terms
	masa_eval_u_source  (x,y,&ufield);
	masa_eval_v_source  (x,y,&vfield);
	masa_eval_e_source  (x,y,&efield);
	masa_eval_rho_source(x,y,&rho);

	//evaluate analytical terms
	masa_eval_u_an        (x,y,&u_an);
	masa_eval_v_an        (x,y,&v_an);
	masa_eval_p_an        (x,y,&p_an);
	masa_eval_rho_an      (x,y,&rho_an);
	  
	// check against maple
	ufield2 = SourceQ_u   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);
	vfield2 = SourceQ_v   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);
	rho2    = SourceQ_rho (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,L);  
	efield2 = SourceQ_e   (x,y,u_0,u_x,u_y,v_0,v_x,v_y,rho_0,rho_x,rho_y,p_0,p_x,p_y,a_px,a_py,a_rhox,a_rhoy,a_ux,a_uy,a_vx,a_vy,Gamma,mu,L);
	
	u_an2   = anQ_u   (x,y,u_0,u_x,u_y,a_ux,a_uy,L);
	v_an2   = anQ_v   (x,y,v_0,v_x,v_y,a_vx,a_vy,L);
	rho_an2 = anQ_rho (x,y,rho_0,rho_x,rho_y,a_rhox,a_rhoy,L);
	p_an2   = anQ_p   (x,y,p_0,p_x,p_y,a_px,a_py,L);
	
	// test the result is roughly zero
	ufield3 = fabs(ufield-ufield2);
	vfield3 = fabs(vfield-vfield2);
	efield3 = fabs(efield-efield2);
	rho     = fabs(rho-rho2);
	
	u_an   = fabs(u_an-u_an2);
	v_an3  = fabs(v_an-v_an2);
	rho_an = fabs(rho_an-rho_an2);
	p_an   = fabs(p_an-p_an2);

	if(ufield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout << "U Field Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << ufield3 << endl;
	    cout << "Source term is:                   " << ufield2 << endl;
	    cout << "MASA term is:                     " << ufield << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(u_an > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout << "U Field Analytical Term\n";
	    cout << "Exceeded Threshold by: " << u_an << endl;
	    cout.precision(16);
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(vfield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout << "V Field Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << vfield3 << endl;
	    cout << "Source term is:                   " << vfield2 << endl;
	    cout << "MASA term is:                     " << vfield << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(v_an3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout << "V Field Analytical Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << v_an3 << endl;
	    cout << "Source term is:        " << v_an2 << endl;
	    cout << "MASA term is:          " << v_an << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(efield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout << "Energy Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << efield3 << endl;
	    cout << "Source term is:        " << efield2 << endl;
	    cout << "MASA term is:          " << efield << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(p_an > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout << "P Field Analytical Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << p_an << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(rho > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout.precision(16);
	    cout << "RHO Source Term\n";
	    cout << "Exceeded Threshold by: " << rho << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

	if(rho_an > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout.precision(16);
	    cout << "RHO Analytical Term\n";
	    cout << "Exceeded Threshold by: " << rho_an << endl;
	    cout << x << " " << y << endl;
	    exit(1);
	  }

      } // done iterating

  // tests passed
  return 0;
}
