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
// axi_euler.cpp: program that tests axisymmetric euler
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h> // for MASA_STRICT_REGRESSION
#include <masa.h>
#include <math.h>

#include <iostream>
#include <stdlib.h>

using namespace MASA;
using namespace std;

const double pi = acos(-1);
const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double nancheck(double x)
{
  if(isnan(x))
    {
      cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
  return 1;
}

double anQ_p(double r,double z,double p_0,double p_1,double rho_0,double rho_1,double u_1,double w_0,double w_1,double a_pr,double a_pz,double a_rhor,double a_rhoz,double a_ur,double a_uz,double a_wr,double a_wz,double pi,double L,double Gamma)
{
  double p_an = p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L);
  return p_an;
}
  
double anQ_u (double r,double z,double p_0,double p_1,double rho_0,double rho_1,double u_1,double w_0,double w_1,double a_pr,double a_pz,double a_rhor,double a_rhoz,double a_ur,double a_uz,double a_wr,double a_wz,double pi,double L,double Gamma)
{
  double u_an = u_1 * (cos(a_ur * pi * r / L) - 0.1e1) * sin(a_uz * pi * z / L);
  return u_an;
} 
 
double anQ_w (double r,double z,double p_0,double p_1,double rho_0,double rho_1,double u_1,double w_0,double w_1,double a_pr,double a_pz,double a_rhor,double a_rhoz,double a_ur,double a_uz,double a_wr,double a_wz,double pi,double L,double Gamma)
{
  double w_an = w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L);
  return w_an;
}

double anQ_rho (double r,double z,double p_0,double p_1,double rho_0,double rho_1,double u_1,double w_0,double w_1,double a_pr,double a_pz,double a_rhor,double a_rhoz,double a_ur,double a_uz,double a_wr,double a_wz,double pi,double L,double Gamma)
{ 
  double rho_an = rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L);
  return rho_an;
}

// ----------------------------------------
//   Source Terms
// ----------------------------------------

double SourceQ_e(double r,double z,double p_0,double p_1,double rho_0,double rho_1,double u_1,double w_0,double w_1,double a_pr,double a_pz,double a_rhor,double a_rhoz,double a_ur,double a_uz,double a_wr,double a_wz,double pi,double L,double Gamma)
{
  double Q_e = Gamma * cos(a_pz * pi * z / L) * cos(a_pr * pi * r / L) * p_1 * u_1 * (cos(a_ur * pi * r / L) - 0.1e1) * sin(a_uz * pi * z / L) * a_pr * pi / L / (Gamma - 0.1e1) - Gamma * sin(a_pz * pi * z / L) * sin(a_pr * pi * r / L) * p_1 * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_pz * pi / L / (Gamma - 0.1e1) - sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * sin(a_uz * pi * z / L) * u_1 * (cos(a_ur * pi * r / L) - 0.1e1) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), 0.2e1) + pow(sin(a_uz * pi * z / L), 0.2e1) * pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_rhor * pi * rho_1 / L / 0.2e1 + cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), 0.2e1) + pow(sin(a_uz * pi * z / L), 0.2e1) * pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_rhoz * pi * rho_1 / L / 0.2e1 - sin(a_uz * pi * z / L) * sin(a_ur * pi * r / L) * (p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L)) * a_ur * pi * u_1 * Gamma / L / (Gamma - 0.1e1) - sin(a_uz * pi * z / L) * sin(a_ur * pi * r / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), 0.2e1) + 0.3e1 * pow(sin(a_uz * pi * z / L), 0.2e1) * pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_ur * pi * u_1 / L / 0.2e1 + (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * u_1 * u_1 * cos(a_uz * pi * z / L) * sin(a_uz * pi * z / L) * pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_uz * pi / L - sin(a_uz * pi * z / L) * (cos(a_ur * pi * r / L) - 0.1e1) * u_1 * w_1 * sin(a_wr * pi * r / L) * sin(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_wr * pi / L + cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * (p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L)) * a_wz * pi * w_1 * Gamma / L / (Gamma - 0.1e1) + cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (0.3e1 * pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), 0.2e1) + pow(sin(a_uz * pi * z / L), 0.2e1) * pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_wz * pi * w_1 / L / 0.2e1 + sin(a_uz * pi * z / L) * (cos(a_ur * pi * r / L) - 0.1e1) * u_1 * (p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L)) * Gamma / (Gamma - 0.1e1) / r + sin(a_uz * pi * z / L) * (cos(a_ur * pi * r / L) - 0.1e1) * u_1 * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), 0.2e1) + pow(sin(a_uz * pi * z / L), 0.2e1) * pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) / r / 0.2e1;
  return(Q_e);
}

double SourceQ_u(double r,double z,double p_0,double p_1,double rho_0,double rho_1,double u_1,double w_0,double w_1,double a_pr,double a_pz,double a_rhor,double a_rhoz,double a_ur,double a_uz,double a_wr,double a_wz,double pi,double L,double Gamma)
{
  double Q_u = p_1 * cos(a_pr * pi * r / L) * cos(a_pz * pi * z / L) * a_pr * pi / L - u_1 * u_1 * pow(sin(a_uz * pi * z / L), 0.2e1) * pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * rho_1 * sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * a_rhor * pi / L + u_1 * (cos(a_ur * pi * r / L) - 0.1e1) * rho_1 * cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) * sin(a_uz * pi * z / L) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_rhoz * pi / L - 0.2e1 * u_1 * u_1 * pow(sin(a_uz * pi * z / L), 0.2e1) * (cos(a_ur * pi * r / L) - 0.1e1) * sin(a_ur * pi * r / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_ur * pi / L + u_1 * (cos(a_ur * pi * r / L) - 0.1e1) * cos(a_uz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_uz * pi / L + u_1 * (cos(a_ur * pi * r / L) - 0.1e1) * w_1 * cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * sin(a_uz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_wz * pi / L + u_1 * u_1 * pow(sin(a_uz * pi * z / L), 0.2e1) * pow(cos(a_ur * pi * r / L) - 0.1e1, 0.2e1) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) / r;
  return(Q_u);  
}

double SourceQ_w(double r,double z,double p_0,double p_1,double rho_0,double rho_1,double u_1,double w_0,double w_1,double a_pr,double a_pz,double a_rhor,double a_rhoz,double a_ur,double a_uz,double a_wr,double a_wz,double pi,double L,double Gamma)
{
  double Q_w = -p_1 * sin(a_pr * pi * r / L) * sin(a_pz * pi * z / L) * a_pz * pi / L - u_1 * sin(a_uz * pi * z / L) * rho_1 * sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * (cos(a_ur * pi * r / L) - 0.1e1) * a_rhor * pi / L + pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), 0.2e1) * rho_1 * cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) * a_rhoz * pi / L - u_1 * sin(a_uz * pi * z / L) * sin(a_ur * pi * r / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_ur * pi / L - u_1 * sin(a_uz * pi * z / L) * w_1 * sin(a_wr * pi * r / L) * sin(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (cos(a_ur * pi * r / L) - 0.1e1) * a_wr * pi / L + (0.2e1 * w_0 + 0.2e1 * w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * w_1 * cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_wz * pi / L + u_1 * sin(a_uz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * (cos(a_ur * pi * r / L) - 0.1e1) / r;
  return(Q_w);
}

double SourceQ_rho(double r,double z,double p_0,double p_1,double rho_0,double rho_1,double u_1,double w_0,double w_1,double a_pr,double a_pz,double a_rhor,double a_rhoz,double a_ur,double a_uz,double a_wr,double a_wz,double pi,double L,double Gamma)
{
  double Q_rho = -(cos(a_ur * pi * r / L) - 0.1e1) * a_rhor * pi * rho_1 * u_1 * sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * sin(a_uz * pi * z / L) / L + (w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L) + w_0) * a_rhoz * pi * rho_1 * cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) / L - (rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) + rho_0) * a_ur * pi * u_1 * sin(a_ur * pi * r / L) * sin(a_uz * pi * z / L) / L + (rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) + rho_0) * a_wz * pi * w_1 * cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) / L + (rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) + rho_0) * (cos(a_ur * pi * r / L) - 0.1e1) * u_1 * sin(a_uz * pi * z / L) / r;
  return(Q_rho);
}

int main()
{  
  //variables
  double R;  
  double p_0;
  double p_1;
  double rho_0;
  double rho_1;
  double u_1;
  double w_0;
  double w_1;
  double a_pr;
  double a_pz;
  double a_rhor;
  double a_rhoz;
  double a_ur;
  double a_uz;
  double a_wr;
  double a_wz;
  double L;
  double mu;
  double Gamma;    
  
  // parameters
  double r;
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
  int nx = 115;  // number of points
  int ny = 68;  
  int lx=3;     // length
  int ly=1; 
  
  double dx=double(lx)/double(nx);
  double dy=double(ly)/double(ny);

  masa_init("axisymmetric_euler","axisymmetric_euler");

  // set params
  masa_init_param();
  
  // get vars
  R      = masa_get_param("R");
  p_0    = masa_get_param("p_0");
  p_1    = masa_get_param("p_1");
  rho_0  = masa_get_param("rho_0");
  rho_1  = masa_get_param("rho_1");
  u_1    = masa_get_param("u_1");
  w_0    = masa_get_param("w_0");
  w_1    = masa_get_param("w_1");
  a_pr   = masa_get_param("a_pr");
  a_pz   = masa_get_param("a_pz");
  a_rhor = masa_get_param("a_rhor");
  a_rhoz = masa_get_param("a_rhoz");
  a_ur   = masa_get_param("a_ur");
  a_uz   = masa_get_param("a_uz");
  a_wr   = masa_get_param("a_wr");
  a_wz   = masa_get_param("a_wz");
  L      = masa_get_param("L");
  Gamma  = masa_get_param("Gamma");
  mu     = masa_get_param("mu");

  // check that all terms have been initialized
  int err = masa_sanity_check();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }  

  // evaluate source terms (2D)
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)    
      {
	r=i*dx;
	z=j*dy;

	//evalulate source terms
	ufield = masa_eval_u_source  (r,z);
	wfield = masa_eval_w_source  (r,z);
	efield = masa_eval_e_source  (r,z);
	rho    = masa_eval_rho_source(r,z);

	//evaluate analytical terms
	u_an = masa_eval_u_an        (r,z);
	w_an = masa_eval_w_an        (r,z);
	p_an = masa_eval_p_an        (r,z);
	rho_an = masa_eval_rho_an    (r,z);
	  	  
	// check against maple
	ufield2 = SourceQ_u   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	wfield2 = SourceQ_w   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	rho2    = SourceQ_rho (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	efield2 = SourceQ_e   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	
	u_an2   = anQ_u   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	w_an2   = anQ_w   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	rho_an2 = anQ_rho (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	p_an2   = anQ_p   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);

	// test the result is roughly zero
	// choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	ufield3 = fabs(ufield-ufield2);
	wfield3 = fabs(wfield-wfield2);
	efield3 = fabs(efield-efield2);
	rho3    = fabs(rho-rho2);
	
	u_an3   = fabs(u_an-u_an2);
	w_an3   = fabs(w_an-w_an2);
	rho_an3 = fabs(rho_an-rho_an2);
	p_an3   = fabs(p_an-p_an2);


#else

	ufield3 = fabs(ufield-ufield2)/fabs(ufield2);
	wfield3 = fabs(wfield-wfield2)/fabs(wfield2);
	efield3 = fabs(efield-efield2)/fabs(efield2);
	rho3    = fabs(rho-rho2)/fabs(rho2);
	
	u_an3   = fabs(u_an-u_an2)/fabs(u_an2);
	w_an3   = fabs(v_an-v_an2)/fabs(w_an2);
	rho_an3 = fabs(rho_an-rho_an2)/fabs(rho_an2);
	p_an3   = fabs(p_an-p_an2)/fabs(p_an2);

#endif

	/*
	nancheck(ufield3);
	nancheck(wfield3);
	nancheck(efield3);
	nancheck(rho3);

	nancheck(u_an3);
	nancheck(w_an3);
	nancheck(rho_an3);
	nancheck(p_an3);
	*/

	if(ufield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout << "U Field Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << ufield3 << endl;
	    cout << "Source term is:                   " << ufield2 << endl;
	    cout << "MASA term is:                     " << ufield << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	if(u_an3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout << "U Field Analytical Term\n";
	    cout << "Exceeded Threshold by: " << u_an << endl;
	    cout.precision(16);
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	if(wfield3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout << "W Field Source Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by:            " << wfield3 << endl;
	    cout << "Source term is:                   " << wfield2 << endl;
	    cout << "MASA term is:                     " << wfield << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }
	
	// this is the 'bad' term: rho_0 is corrupting somehow
	if(w_an3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout << "W Field Analytical Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << w_an3 << endl;
	    cout << "Source term is:        " << w_an2 << endl;
	    cout << "MASA term is:          " << w_an << endl;
	    cout << r << " " << z << endl;
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
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	if(p_an3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout << "P Field Analytical Term\n";
	    cout.precision(16);
	    cout << "Exceeded Threshold by: " << p_an << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	if(rho3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout.precision(16);
	    cout << "RHO Source Term\n";
	    cout << "Exceeded Threshold by: " << rho << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }

	if(rho_an3 > threshold)
	  {
	    cout << "\nMASA REGRESSION TEST FAILED: Axisymmetric Euler\n";
	    cout.precision(16);
	    cout << "RHO Analytical Term\n";
	    cout << "Exceeded Threshold by: " << rho_an << endl;
	    cout << r << " " << z << endl;
	    exit(1);
	  }

      } // done iterating

  // tests passed
  return 0;
}
