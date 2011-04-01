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
// $Author:
// $Id:
//
// euler_chem.cpp: program that tests euler with chemistry
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <tests.h>

using namespace MASA;
using namespace std;

/* ------------------------------------------------
 *
 *
 *   Source Terms
 *
 * -----------------------------------------------
 */ 

template<typename Scalar>
Scalar SourceQ_e (
  Scalar x,
  Scalar t,
  Scalar u_0,
  Scalar u_x,
  Scalar u_t,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_t,
  Scalar p_0,
  Scalar p_x,
  Scalar p_t,
  Scalar a_px,
  Scalar a_pt,
  Scalar a_rhox,
  Scalar a_rhot,
  Scalar a_ux,
  Scalar a_ut,
  Scalar L,
  Scalar Gamma)
{
  Scalar pi = acos(-1);
  Scalar Q_e_t;
  Scalar RHO;
  Scalar U;
  Scalar P;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  P = p_0 + p_x * cos(a_px * pi * x / L) + p_t * cos(a_pt * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);

  Q_e_t = cos(a_rhot * pi * t / L) * a_rhot * pi * rho_t * U * U / L / 0.2e1 - sin(a_ut * pi * t / L) * a_ut * pi * u_t * RHO * U / L - sin(a_pt * pi * t / L) * a_pt * pi * p_t / (Gamma - 0.1e1) / L + cos(a_rhox * pi * x / L) * pow(U, 0.3e1) * a_rhox * pi * rho_x / L / 0.2e1 + cos(a_ux * pi * x / L) * P * a_ux * pi * u_x * Gamma / (Gamma - 0.1e1) / L + 0.3e1 / 0.2e1 * cos(a_ux * pi * x / L) * RHO * U * U * a_ux * pi * u_x / L - sin(a_px * pi * x / L) * U * a_px * pi * p_x * Gamma / (Gamma - 0.1e1) / L;
  return(Q_e_t);

}

template<typename Scalar>
Scalar SourceQ_u (
  Scalar x,
  Scalar t,
  Scalar u_0,
  Scalar u_x,
  Scalar u_t,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_t,
  Scalar p_x,
  Scalar p_t,
  Scalar a_px,
  Scalar a_rhox,
  Scalar a_rhot,
  Scalar a_ux,
  Scalar a_ut,
  Scalar L)
{
  Scalar pi = acos(-1);
  Scalar Q_u_t;
  Scalar RHO;
  Scalar U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);

  Q_u_t = cos(a_rhox * pi * x / L) * a_rhox * pi * rho_x * U * U / L + 0.2e1 * cos(a_ux * pi * x / L) * RHO * a_ux * pi * u_x * U / L + cos(a_rhot * pi * t / L) * a_rhot * pi * rho_t * U / L - sin(a_ut * pi * t / L) * a_ut * pi * u_t * RHO / L - sin(a_px * pi * x / L) * a_px * pi * p_x / L;
  return(Q_u_t);

}

template<typename Scalar>
Scalar SourceQ_rho (
  Scalar x,
  Scalar t,
  Scalar u_0,
  Scalar u_x,
  Scalar u_t,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_t,
  Scalar p_x,
  Scalar p_t,
  Scalar a_px,
  Scalar a_rhox,
  Scalar a_rhot,
  Scalar a_ux,
  Scalar a_ut,
  Scalar L)
{
  Scalar pi = acos(-1);
  Scalar Q_rho_t;
  Scalar RHO;
  Scalar U;

  RHO = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);

  Q_rho_t = cos(a_rhot * pi * t / L) * a_rhot * pi * rho_t / L + cos(a_ux * pi * x / L) * RHO * a_ux * pi * u_x / L + cos(a_rhox * pi * x / L) * U * a_rhox * pi * rho_x / L;
  return(Q_rho_t);

}

/* ------------------------------------------------
 *
 *
 *   Analytical terms
 *
 * -----------------------------------------------
 */ 

template<typename Scalar>
Scalar anQ_p (Scalar x,Scalar t,Scalar p_0,Scalar p_x,Scalar p_t,Scalar a_px,Scalar a_pt,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar exact_p;
  exact_p = p_0 + p_x * cos(a_px * pi * x / L) + p_t * cos(a_pt * pi * t / L);
  return exact_p;
}
  
template<typename Scalar>
Scalar anQ_u (Scalar x,Scalar t,Scalar u_0,Scalar u_x,Scalar a_ux,Scalar u_t,Scalar a_ut,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar exact_u;
  exact_u = u_0 + u_x * sin(a_ux * pi * x / L) + u_t * cos(a_ut * pi * t / L);
  return exact_u;
} 
 
template<typename Scalar>
Scalar anQ_rho (Scalar x,Scalar t,Scalar rho_0,Scalar rho_x,Scalar a_rhox,Scalar rho_t,Scalar a_rhot,Scalar L)
{ 
  Scalar pi = acos(-1);
  Scalar exact_rho;
  exact_rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_t * sin(a_rhot * pi * t / L);
  return exact_rho;
}

template<typename Scalar>
int run_regression()
{

  Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();

  Scalar u_0;
  Scalar u_x;
  Scalar u_t;    
  Scalar a_ux;
  Scalar a_ut;

  Scalar rho_0;
  Scalar rho_x;
  Scalar rho_t;
  Scalar a_rhox;
  Scalar a_rhot;

  Scalar p_0;
  Scalar p_x;
  Scalar p_t;
  Scalar a_px;
  Scalar a_pt;

  Scalar Gamma;
  Scalar L;    
  
  // solutions
  Scalar ufield,ufield2,ufield3;
  Scalar efield,efield2,efield3;
  Scalar rho,rho2,rho3;
  //Scalar gradx,grady,gradz,gradp,gradrho;

  Scalar exact_u,exact_u2,exact_u3;
  Scalar exact_p,exact_p2,exact_p3;
  Scalar exact_rho,exact_rho2,exact_rho3;

  // parameters
  Scalar x;
  Scalar t;

  //problem size
  int nx    = 200;  // number of points
  int lx    = 10;     // length
  Scalar dx = Scalar(lx)/Scalar(nx);

  int nt    = 100;  // number of points
  int lt    = 10;     // length
  Scalar dt = Scalar(lt)/Scalar(nt);

  // initalize
  masa_init<Scalar>("euler-transient-test","euler_transient_1d");

  // initialize the default parameters
  masa_init_param<Scalar>();

  // get defaults for comparison to source terms
  u_0 = masa_get_param<Scalar>("u_0");
  u_x = masa_get_param<Scalar>("u_x");
  u_t = masa_get_param<Scalar>("u_t");

  rho_0 = masa_get_param<Scalar>("rho_0");
  rho_x = masa_get_param<Scalar>("rho_x");
  rho_t = masa_get_param<Scalar>("rho_t");

  p_0 = masa_get_param<Scalar>("p_0");
  p_x = masa_get_param<Scalar>("p_x");
  p_t = masa_get_param<Scalar>("p_t");

  a_px = masa_get_param<Scalar>("a_px");
  a_pt = masa_get_param<Scalar>("a_pt");

  a_rhox = masa_get_param<Scalar>("a_rhox");
  a_rhot = masa_get_param<Scalar>("a_rhot");

  a_ux = masa_get_param<Scalar>("a_ux");
  a_ut = masa_get_param<Scalar>("a_ut");

  Gamma = masa_get_param<Scalar>("Gamma");
  L     = masa_get_param<Scalar>("L");

  // check that all terms have been initialized
  masa_sanity_check<Scalar>();

  // evaluate MMS (1D+time)
  for(int i=0;i<nx;i++)
    for(int j=0;j<nt;j++)
      {
	t=j*dt;
	x=i*dx;
	
	// evalulate source terms
	ufield = masa_eval_source_rho_u  <Scalar>(x,t);
	efield = masa_eval_source_rho_e  <Scalar>(x,t);
	rho    = masa_eval_source_rho    <Scalar>(x,t);

	exact_u = masa_eval_exact_u      <Scalar>(x,t);
	exact_p = masa_eval_exact_p      <Scalar>(x,t);
	exact_rho = masa_eval_exact_rho  <Scalar>(x,t);

	efield2 = SourceQ_e (x,t,u_0,
			     u_x,
			     u_t,
			     rho_0,
			     rho_x,
			     rho_t,
			     p_0,
			     p_x,
			     p_t,
			     a_px,
			     a_pt,
			     a_rhox,
			     a_rhot,
			     a_ux,
			     a_ut,
			     L,
			     Gamma);

	ufield2 = SourceQ_u (x,
			     t,
			     u_0,
			     u_x,
			     u_t,
			     rho_0,
			     rho_x,
			     rho_t,
			     p_x,
			     p_t,
			     a_px,
			     a_rhox,
			     a_rhot,
			     a_ux,
			     a_ut,
			     L);

	rho2 = SourceQ_rho (x,
			    t,
			    u_0,
			    u_x,
			    u_t,
			    rho_0,
			    rho_x,
			    rho_t,
			    p_x,
			    p_t,
			    a_px,
			    a_rhox,
			    a_rhot,
			    a_ux,
			    a_ut,
			    L);


	exact_p2   = anQ_p (x,t,p_0,p_x,p_t,a_px,a_pt,L);
	exact_u2   = anQ_u (x,t,u_0,u_x,a_ux,u_t,a_ut, L);
	exact_rho2 = anQ_rho(x,t,rho_0,rho_x,a_rhox,rho_t,a_rhot,L);

	// test the result is roughly zero
	// choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

	ufield3 = fabs(ufield-ufield2);
	efield3 = fabs(efield-efield2);
	rho3    = fabs(rho-rho2);

	exact_u3   = fabs(exact_u-exact_u2);
	exact_rho3 = fabs(exact_rho-exact_rho2);
	exact_p3   = fabs(exact_p-exact_p2);

#else

	ufield3 = fabs(ufield-ufield2)/fabs(ufield2);
	efield3 = fabs(efield-efield2)/fabs(efield2);
	rho3    = fabs(rho-rho2)/fabs(rho2);

	exact_u3   = fabs(exact_u-exact_u2)/fabs(exact_u2);
	exact_rho3 = fabs(exact_rho-exact_rho2)/fabs(exact_rho2);
	exact_p3   = fabs(exact_p-exact_p2)/fabs(exact_p2);

#endif

	threshcheck(ufield3,threshold);
	threshcheck(efield3,threshold);
	threshcheck(rho3,threshold);
	threshcheck(exact_u3,threshold);
	threshcheck(exact_rho3,threshold);
	threshcheck(exact_p3,threshold);

      }// done with time and space loop
  
  return 0;

}

int main()
{
  int err=0;

  err += run_regression<double>();
  //err += run_regression<long double>();

  return err;
}

