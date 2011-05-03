// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011 The PECOS Development Team
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

#include <tests.h>
#include <cmath>

using namespace MASA;
using namespace std;

template<typename Scalar>
Scalar anQ_p(Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar exact_p = p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L);
  return exact_p;
}
  
template<typename Scalar>
Scalar anQ_u (Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar exact_u = u_1 * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * sin(a_uz * pi * z / L);
  return exact_u;
} 
 
template<typename Scalar>
Scalar anQ_w (Scalar r,Scalar z,Scalar w_0,Scalar w_1,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L)
{
  Scalar exact_w = w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L);
  return exact_w;
}

template<typename Scalar>
Scalar anQ_rho (Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{ 
  Scalar exact_rho = rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L);
  return exact_rho;
}

// ----------------------------------------
//   Source Terms
// ----------------------------------------

template<typename Scalar>
Scalar SourceQ_e(Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar Q_e = Gamma * cos(a_pz * pi * z / L) * cos(a_pr * pi * r / L) * p_1 * u_1 * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * sin(a_uz * pi * z / L) * a_pr * pi / L / (Gamma - Scalar(0.1e1)) - Gamma * sin(a_pz * pi * z / L) * sin(a_pr * pi * r / L) * p_1 * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_pz * pi / L / (Gamma - Scalar(0.1e1)) - sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * sin(a_uz * pi * z / L) * u_1 * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(0.2e1)) + pow(sin(a_uz * pi * z / L), Scalar(0.2e1)) * pow(cos(a_ur * pi * r / L) - 0.1e1, Scalar(0.2e1)) * u_1 * u_1) * a_rhor * pi * rho_1 / L / Scalar(0.2e1) + cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(0.2e1)) + pow(sin(a_uz * pi * z / L), Scalar(0.2e1)) * pow(cos(a_ur * pi * r / L) - 0.1e1, Scalar(0.2e1)) * u_1 * u_1) * a_rhoz * pi * rho_1 / L / Scalar(0.2e1) - sin(a_uz * pi * z / L) * sin(a_ur * pi * r / L) * (p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L)) * a_ur * pi * u_1 * Gamma / L / (Gamma - Scalar(0.1e1)) - sin(a_uz * pi * z / L) * sin(a_ur * pi * r / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(0.2e1)) + Scalar(0.3e1) * pow(sin(a_uz * pi * z / L), Scalar(0.2e1)) * pow(cos(a_ur * pi * r / L) - 0.1e1, Scalar(0.2e1)) * u_1 * u_1) * a_ur * pi * u_1 / L / Scalar(0.2e1) + (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * u_1 * u_1 * cos(a_uz * pi * z / L) * sin(a_uz * pi * z / L) * pow(cos(a_ur * pi * r / L) - 0.1e1, Scalar(0.2e1)) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_uz * pi / L - sin(a_uz * pi * z / L) * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * u_1 * w_1 * sin(a_wr * pi * r / L) * sin(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_wr * pi / L + cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * (p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L)) * a_wz * pi * w_1 * Gamma / L / (Gamma - Scalar(0.1e1)) + cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (Scalar(0.3e1) * pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(0.2e1)) + pow(sin(a_uz * pi * z / L), Scalar(0.2e1)) * pow(cos(a_ur * pi * r / L) - 0.1e1, Scalar(0.2e1)) * u_1 * u_1) * a_wz * pi * w_1 / L / Scalar(0.2e1) + sin(a_uz * pi * z / L) * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * u_1 * (p_0 + p_1 * sin(a_pr * pi * r / L) * cos(a_pz * pi * z / L)) * Gamma / (Gamma - Scalar(0.1e1)) / r + sin(a_uz * pi * z / L) * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * u_1 * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(0.2e1)) + pow(sin(a_uz * pi * z / L), Scalar(0.2e1)) * pow(cos(a_ur * pi * r / L) - 0.1e1, Scalar(0.2e1)) * u_1 * u_1) / r / Scalar(0.2e1);
  return(Q_e);
}

template<typename Scalar>
Scalar SourceQ_u(Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar Q_u = p_1 * cos(a_pr * pi * r / L) * cos(a_pz * pi * z / L) * a_pr * pi / L - u_1 * u_1 * pow(sin(a_uz * pi * z / L), Scalar(0.2e1)) * pow(cos(a_ur * pi * r / L) - 0.1e1, Scalar(0.2e1)) * rho_1 * sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * a_rhor * pi / L + u_1 * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * rho_1 * cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) * sin(a_uz * pi * z / L) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_rhoz * pi / L - Scalar(0.2e1) * u_1 * u_1 * pow(sin(a_uz * pi * z / L), Scalar(0.2e1)) * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * sin(a_ur * pi * r / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_ur * pi / L + u_1 * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * cos(a_uz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_uz * pi / L + u_1 * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * w_1 * cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * sin(a_uz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_wz * pi / L + u_1 * u_1 * pow(sin(a_uz * pi * z / L), Scalar(0.2e1)) * pow(cos(a_ur * pi * r / L) - 0.1e1, Scalar(0.2e1)) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) / r;
  return(Q_u);  
}

template<typename Scalar>
Scalar SourceQ_w(Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar Q_w = -p_1 * sin(a_pr * pi * r / L) * sin(a_pz * pi * z / L) * a_pz * pi / L - u_1 * sin(a_uz * pi * z / L) * rho_1 * sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * a_rhor * pi / L + pow(w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L), Scalar(0.2e1)) * rho_1 * cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) * a_rhoz * pi / L - u_1 * sin(a_uz * pi * z / L) * sin(a_ur * pi * r / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * a_ur * pi / L - u_1 * sin(a_uz * pi * z / L) * w_1 * sin(a_wr * pi * r / L) * sin(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * a_wr * pi / L + (Scalar(0.2e1) * w_0 + Scalar(0.2e1) * w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * w_1 * cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * a_wz * pi / L + u_1 * sin(a_uz * pi * z / L) * (rho_0 + rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L)) * (w_0 + w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L)) * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) / r;
  return(Q_w);
}

template<typename Scalar>
Scalar SourceQ_rho(Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma)
{
  Scalar Q_rho = -(cos(a_ur * pi * r / L) - Scalar(0.1e1)) * a_rhor * pi * rho_1 * u_1 * sin(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) * sin(a_uz * pi * z / L) / L + (w_1 * cos(a_wr * pi * r / L) * sin(a_wz * pi * z / L) + w_0) * a_rhoz * pi * rho_1 * cos(a_rhor * pi * r / L) * cos(a_rhoz * pi * z / L) / L - (rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) + rho_0) * a_ur * pi * u_1 * sin(a_ur * pi * r / L) * sin(a_uz * pi * z / L) / L + (rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) + rho_0) * a_wz * pi * w_1 * cos(a_wr * pi * r / L) * cos(a_wz * pi * z / L) / L + (rho_1 * cos(a_rhor * pi * r / L) * sin(a_rhoz * pi * z / L) + rho_0) * (cos(a_ur * pi * r / L) - Scalar(0.1e1)) * u_1 * sin(a_uz * pi * z / L) / r;
  return(Q_rho);
}

template<typename Scalar>
int run_regression()
{  
  //variables
  Scalar p_0;
  Scalar p_1;
  Scalar rho_0;
  Scalar rho_1;
  Scalar u_1;
  Scalar w_0;
  Scalar w_1;
  Scalar a_pr;
  Scalar a_pz;
  Scalar a_rhor;
  Scalar a_rhoz;
  Scalar a_ur;
  Scalar a_uz;
  Scalar a_wr;
  Scalar a_wz;
  Scalar L;
  Scalar Gamma;    
  
  // parameters
  Scalar r;
  Scalar z;

  // solutions
  Scalar ufield,ufield2,ufield3;
  Scalar wfield,wfield2,wfield3;
  Scalar efield,efield2,efield3;
  Scalar rho,rho2,rho3;

  Scalar exact_u,exact_u2,exact_u3;
  Scalar exact_w,exact_w2,exact_w3;
  Scalar exact_p,exact_p2,exact_p3;
  Scalar exact_rho,exact_rho2,exact_rho3;

  const Scalar pi = acos(Scalar(-1));

  // initalize
  int nx = 115;  // number of points
  int ny = 68;  
  int lx=3;     // length
  int ly=1; 
  
  Scalar dx=Scalar(lx)/Scalar(nx);
  Scalar dy=Scalar(ly)/Scalar(ny);

  masa_init<Scalar>("axisymmetric_euler","axisymmetric_euler");

  // set params
  masa_init_param<Scalar>();
  
  // get vars
  p_0    = masa_get_param<Scalar>("p_0");
  p_1    = masa_get_param<Scalar>("p_1");
  rho_0  = masa_get_param<Scalar>("rho_0");
  rho_1  = masa_get_param<Scalar>("rho_1");
  u_1    = masa_get_param<Scalar>("u_1");
  w_0    = masa_get_param<Scalar>("w_0");
  w_1    = masa_get_param<Scalar>("w_1");
  a_pr   = masa_get_param<Scalar>("a_pr");
  a_pz   = masa_get_param<Scalar>("a_pz");
  a_rhor = masa_get_param<Scalar>("a_rhor");
  a_rhoz = masa_get_param<Scalar>("a_rhoz");
  a_ur   = masa_get_param<Scalar>("a_ur");
  a_uz   = masa_get_param<Scalar>("a_uz");
  a_wr   = masa_get_param<Scalar>("a_wr");
  a_wz   = masa_get_param<Scalar>("a_wz");
  L      = masa_get_param<Scalar>("L");
  Gamma  = masa_get_param<Scalar>("Gamma");

  // check that all terms have been initialized
  int err = masa_sanity_check<Scalar>();
  if(err != 0)
    {
      cout << "MASA :: Sanity Check Failed!\n";
      exit(1);
    }  

  // evaluate source terms (2D)
  for(int i=1;i<nx;i++)    // this is the radial term -- thus, do not start at 0!
    for(int j=1;j<ny;j++)  // z component  
      {
	r=i*dx;
	z=j*dy;

	//evalulate source terms
	ufield = masa_eval_source_rho_u<Scalar>  (r,z);
	wfield = masa_eval_source_rho_w<Scalar>  (r,z);
	efield = masa_eval_source_rho_e<Scalar>  (r,z);
	rho    = masa_eval_source_rho<Scalar>(r,z);

	//evaluate analytical terms
	exact_u = masa_eval_exact_u<Scalar>        (r,z);
	exact_w = masa_eval_exact_w<Scalar>        (r,z);
	exact_p = masa_eval_exact_p<Scalar>        (r,z);
	exact_rho = masa_eval_exact_rho<Scalar>    (r,z);
	  	  
	// check against maple
	ufield2 = SourceQ_u   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	wfield2 = SourceQ_w   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	rho2    = SourceQ_rho (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	efield2 = SourceQ_e   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	
	exact_u2   = anQ_u   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	exact_rho2 = anQ_rho (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	exact_p2   = anQ_p   (r, z, p_0, p_1, rho_0, rho_1, u_1, w_0, w_1, a_pr, a_pz, a_rhor, a_rhoz, a_ur, a_uz, a_wr, a_wz, pi, L, Gamma);
	exact_w2   = anQ_w   (r, z, w_0, w_1, a_wr, a_wz, pi, L);

	// test the result is roughly zero
	// choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	ufield3 = std::abs(ufield-ufield2);
	wfield3 = std::abs(wfield-wfield2);
	efield3 = std::abs(efield-efield2);
	rho3    = std::abs(rho-rho2);
	
	exact_u3   = std::abs(exact_u-exact_u2);
	exact_w3   = std::abs(exact_w-exact_w2);
	exact_rho3 = std::abs(exact_rho-exact_rho2);
	exact_p3   = std::abs(exact_p-exact_p2);

#else

	ufield3 = std::abs(ufield-ufield2)/std::abs(ufield2);
	wfield3 = std::abs(wfield-wfield2)/std::abs(wfield2);
	efield3 = std::abs(efield-efield2)/std::abs(efield2);
	rho3    = std::abs(rho-rho2)/std::abs(rho2);
	
	exact_u3   = std::abs(exact_u-exact_u2)/std::abs(exact_u2);
	exact_w3   = std::abs(exact_w-exact_w2)/std::abs(exact_w2);
	exact_rho3 = std::abs(exact_rho-exact_rho2)/std::abs(exact_rho2);
	exact_p3   = std::abs(exact_p-exact_p2)/std::abs(exact_p2);

#endif

	threshcheck(ufield3);
	threshcheck(wfield3);
	threshcheck(efield3);
	threshcheck(rho3);

	threshcheck(exact_u3);
	threshcheck(exact_w3);
	threshcheck(exact_rho3);
	threshcheck(exact_p3);

      } // done iterating

  // tests passed
  return 0;
}

int main()
{
  int err=0;

  err += run_regression<double>();
  err += run_regression<long double>();

  return err;
}
