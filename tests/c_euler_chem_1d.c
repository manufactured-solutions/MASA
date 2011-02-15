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
// $Author: karl $
// $Id: c_misc.c 13969 2010-09-23 13:55:56Z karl $
//
// c_purge.c :program that tests masa purge function
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double MASA_DEFAULT = -12345.67;
const double threshold = 1.0e-15; // should be small enough to catch any obvious problems

double threshcheck(double x,double thresh)
{
  
  if(x > thresh)
    {
      printf("\nCMASA REGRESSION TEST FAILED: Euler-1d + chemistry\n");
      printf("Exceeded Threshold by: %g \n",x);
      exit(1);
    }
  return 1;  
}

double SourceQ_rho_u(double x, 
		     double R_N, 
		     double rho_N_0, 
		     double rho_N_x,
		     double a_rho_N_x, 
		     double rho_N2_0, 
		     double rho_N2_x,
		     double a_rho_N2_x,
		     double L,
		     double u_0,
		     double u_x,
		     double a_ux,
		     double T_0,
		     double T_x,
		     double a_Tx)
{
  double Q_u;
  double RHO;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;

  double pi = acos(-1);

  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  T = T_0 + T_x * cos(a_Tx * pi * x / L);

  Q_u = a_rho_N_x * pi * rho_N_x * U * U * cos(a_rho_N_x * pi * x / L) / L - a_rho_N2_x * pi * rho_N2_x * U * U * sin(a_rho_N2_x * pi * x / L) / L - a_Tx * pi * T_x * R_N * RHO_N * sin(a_Tx * pi * x / L) / L - a_Tx * pi * T_x * R_N * RHO_N2 * sin(a_Tx * pi * x / L) / L / 0.2e1 + 0.2e1 * a_ux * pi * u_x * RHO * U * cos(a_ux * pi * x / L) / L - (-0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * pi * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * pi * x / L)) * pi * R_N * T / L / 0.2e1;
  return(Q_u);
}

double SourceQ_rho_e(double x,
		     double R_N,
		     double R_N2,
		     double h0_N,
		     double h0_N2,
		     double theta_v_N2,
		     double rho_N_0, 
		     double rho_N_x,
		     double a_rho_N_x, 
		     double rho_N2_0, 
		     double rho_N2_x,
		     double a_rho_N2_x,
		     double L,
		     double u_0,
		     double u_x,
		     double a_ux,
		     double T_0,
		     double T_x,
		     double a_Tx)
{
  double Q_e;
  double RHO;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  double alpha;
  double E_vib_N2;

  double pi = acos(-1);

  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  T = T_0 + T_x * cos(a_Tx * pi * x / L);
  alpha = exp(theta_v_N2 / T);
  E_vib_N2 = R_N2 * theta_v_N2 / (alpha - 0.1e1);

  Q_e = -a_Tx * pi * T_x * alpha * theta_v_N2 * E_vib_N2 * RHO_N2 * U * sin(a_Tx * pi * x / L) / L / (alpha - 0.1e1) * pow(T, -0.2e1) - 0.5e1 / 0.2e1 * a_Tx * pi * T_x * R_N * RHO_N * U * sin(a_Tx * pi * x / L) / L - 0.7e1 / 0.4e1 * a_Tx * pi * T_x * R_N * RHO_N2 * U * sin(a_Tx * pi * x / L) / L + 0.5e1 / 0.2e1 * a_ux * pi * u_x * R_N * RHO_N * T * cos(a_ux * pi * x / L) / L + 0.7e1 / 0.4e1 * a_ux * pi * u_x * R_N * RHO_N2 * T * cos(a_ux * pi * x / L) / L + 0.3e1 / 0.2e1 * a_ux * pi * u_x * RHO * U * U * cos(a_ux * pi * x / L) / L - a_rho_N2_x * pi * rho_N2_x * E_vib_N2 * U * sin(a_rho_N2_x * pi * x / L) / L + a_ux * pi * u_x * E_vib_N2 * RHO_N2 * cos(a_ux * pi * x / L) / L + (h0_N * RHO_N + h0_N2 * RHO_N2) * a_ux * pi * u_x * cos(a_ux * pi * x / L) / L - (-0.10e2 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * pi * x / L) + 0.7e1 * a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * pi * x / L)) * pi * R_N * U * T / L / 0.4e1 - (-a_rho_N_x * rho_N_x * cos(a_rho_N_x * pi * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * pi * x / L)) * pi * pow(U, 0.3e1) / L / 0.2e1 - (-a_rho_N_x * rho_N_x * h0_N * cos(a_rho_N_x * pi * x / L) + a_rho_N2_x * rho_N2_x * h0_N2 * sin(a_rho_N2_x * pi * x / L)) * pi * U / L;

  return(Q_e);
}

double SourceQ_rho_N(double x,
		     double M_N,
		     double h0_N,
		     double h0_N2,
		     double Cf1_N,
		     double Cf1_N2,
		     double etaf1_N,
		     double etaf1_N2,
		     double Ea_N,
		     double Ea_N2,
		     //double Function_to_Calculate_K,
		     double (*in_func)(double),
		     double R_N,
		     double R_N2,
		     double theta_v_N2,
		     double rho_N_0, 
		     double rho_N_x,
		     double a_rho_N_x, 
		     double rho_N2_0, 
		     double rho_N2_x,
		     double a_rho_N2_x,
		     double L,
		     double u_0,
		     double u_x,
		     double a_ux,
		     double T_0,
		     double T_x,
		     double a_Tx,
		     double R)

{
  double Q_rho_N;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  double kf1_N;
  double kf1_N2;
  double K_eq;

  double pi = acos(-1);

  // Equilibrium Konstant
  T = T_0 + T_x * cos(a_Tx * pi * x / L);
  K_eq = in_func(T);

  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  kf1_N = Cf1_N * pow(T, etaf1_N) * exp(-Ea_N / R / T);
  kf1_N2 = Cf1_N2 * pow(T, etaf1_N2) * exp(-Ea_N2 / R / T);

  Q_rho_N = a_rho_N_x * pi * rho_N_x * U * cos(a_rho_N_x * pi * x / L) / L + a_ux * pi * u_x * RHO_N * cos(a_ux * pi * x / L) / L + (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N * RHO_N * pow(M_N, -0.2e1) / K_eq - (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N2 / M_N / 0.2e1;

  return(Q_rho_N);
}

double SourceQ_rho_N2(double x,
		      double M_N,
		      double h0_N,
		      double h0_N2,
		      double Cf1_N,
		      double Cf1_N2,
		      double etaf1_N,
		      double etaf1_N2,
		      double Ea_N,
		      double Ea_N2,
		      //double Function_to_Calculate_K,
		      double (*in_func)(double),
		      double K,
		      double R_N,
		      double R_N2,
		      double theta_v_N2,
		      double rho_N_0, 
		      double rho_N_x,
		      double a_rho_N_x, 
		      double rho_N2_0, 
		      double rho_N2_x,
		      double a_rho_N2_x,
		      double L,
		      double u_0,
		      double u_x,
		      double a_ux,
		      double T_0,
		      double T_x,
		      double a_Tx,
		      double R)
{
  
  double Q_rho_N2;
  double RHO_N;
  double RHO_N2;
  double U;
  double T;
  double kf1_N;
  double kf1_N2;
  double K_eq;

  double pi = acos(-1);

  // Equilibrium Konstant
  T = T_0 + T_x * cos(a_Tx * pi * x / L);
  K_eq = in_func(T);

  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  kf1_N = Cf1_N * pow(T, etaf1_N) * exp(-Ea_N / R / T);
  kf1_N2 = Cf1_N2 * pow(T, etaf1_N2) * exp(-Ea_N2 / R / T);

  Q_rho_N2 = -a_rho_N2_x * pi * rho_N2_x * U * sin(a_rho_N2_x * pi * x / L) / L + a_ux * pi * u_x * RHO_N2 * cos(a_ux * pi * x / L) / L - (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N * RHO_N * pow(M_N, -0.2e1) / K_eq + (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N2 / M_N / 0.2e1;
  return(Q_rho_N2);
}

// ----------------------------------------
//   Analytical Terms
// ----------------------------------------

double anQ_t(double x,double T_0,double T_x,double a_Tx,double L)
{
  double pi = acos(-1);
  double T_an = T_0 + T_x * cos(a_Tx * pi * x / L);
  return T_an;
}

double anQ_u(double x,double u_0,double u_x,double a_ux,double L)
{
  double pi = acos(-1);
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L);
  return u_an;
}

double anQ_rho(double x,double rho_N_0,double rho_N_x,double a_rho_N_x,
	       double L,double rho_N2_0,double rho_N2_x,double a_rho_N2_x)
{  
  double pi = acos(-1);
  double rho_an = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  return rho_an;
}

double anQ_rho_N(double x,double rho_N_0,double rho_N_x,double a_rho_N_x,double L)
{
  double pi = acos(-1);
  double rho_an_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L);
  return rho_an_N;
}

double anQ_rho_N2(double x,double rho_N2_0,double rho_N2_x,double a_rho_N2_x,double L)
{
  double pi = acos(-1);
  double rho_an_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  return rho_an_N2;
}

// ----------------------------------------
//   Regresssion
// ----------------------------------------

double temp_function(double T)
{
  // hackish functional here
  // This is an eyeballed fit (focusing on the 5000K-6000K range) 
  // for the equilibrium constant for N2->N+N dissociation
  double K = exp(4+(T-6000)/500);
  return K;
}


int main()
{

  //variables 
  double u_0;
  double u_x;
  double a_ux;
  double L;

  double R;
  double Cf1_N;
  double Cf1_N2;
  double etaf1_N;
  double etaf1_N2;
  double Ea_N;
  double Ea_N2;
  double R_N;
  double R_N2;
  double theta_v_N2;
  double M_N;
  double h0_N;
  double h0_N2;
  double K;

  double rho_N_0;
  double rho_N_x;
  double a_rho_N_x;

  double rho_N2_0;
  double rho_N2_x;
  double a_rho_N2_x;

  double T_0;
  double T_x;
  double a_Tx;

  // parameters
  double x;
  int i;

  //problem size
  int nx = 200;  // number of points
  int lx=10;     // length
  double dx=(double)(lx)/(double)(nx);

  // solutions
  double ufield,ufield2,ufield3;
  double efield,efield2,efield3;
  double N,N2,N3;
  double Ntwo,Ntwo2,Ntwo3;
  double gradx,grady,gradz,gradp,gradrho;

  double exact_t,exact_t2,exact_t3;
  double exact_u,exact_u2,exact_u3;
  double exact_rho,exact_rho2,exact_rho3;
  double exact_N,exact_N2,exact_N3;
  double exact_Ntwo,exact_Ntwo2,exact_Ntwo3;

  // initalize
  cmasa_init("euler-chemistry-test","euler_chem_1d");

  // initialize the default parameters
  cmasa_init_param();

  // get defaults for comparison to source terms
  u_0  = cmasa_get_param("u_0");
  u_x  = cmasa_get_param("u_x");
  a_ux = cmasa_get_param("a_ux");
  L    = cmasa_get_param("L");
  R    = cmasa_get_param("R");

  Cf1_N   = cmasa_get_param("Cf1_N");
  Cf1_N2  = cmasa_get_param("Cf1_N2");
  etaf1_N  = cmasa_get_param("etaf1_N");
  etaf1_N2 = cmasa_get_param("etaf1_N2");

  Ea_N  = cmasa_get_param("Ea_N");
  Ea_N2 = cmasa_get_param("Ea_N2");

  R_N   = cmasa_get_param("R_N");
  R_N2  = cmasa_get_param("R_N2");

  theta_v_N2 = cmasa_get_param("theta_v_N2");
  M_N   = cmasa_get_param("M_N");
  h0_N  = cmasa_get_param("h0_N");
  h0_N2 = cmasa_get_param("h0_N2");
  K     = cmasa_get_param("K");

  rho_N_0   = cmasa_get_param("rho_N_0");
  rho_N_x   = cmasa_get_param("rho_N_x");
  a_rho_N_x = cmasa_get_param("a_rho_N_x");

  rho_N2_0   = cmasa_get_param("rho_N2_0");  
  rho_N2_x   = cmasa_get_param("rho_N2_x");
  a_rho_N2_x = cmasa_get_param("a_rho_N2_x");

  T_0  = cmasa_get_param("T_0");
  T_x  = cmasa_get_param("T_x");
  a_Tx = cmasa_get_param("a_Tx");

  // check that all terms have been initialized
  cmasa_sanity_check();

  // evaluate MMS (1D)
  for(i=0;i<nx;i++)
    {
      x=i*dx;

      // evalulate source terms
      ufield = cmasa_eval_1d_source_rho_u  (x);
      efield = cmasa_eval_1d_source_rho_e  (x);
      N      = cmasa_eval_1d_source_rho_N  (x,&temp_function);
      Ntwo   = cmasa_eval_1d_source_rho_N2 (x,&temp_function);

      // evaluate analytical solution terms
      exact_t    = cmasa_eval_1d_exact_t     (x);
      exact_u    = cmasa_eval_1d_exact_u     (x);
      exact_rho  = cmasa_eval_1d_exact_rho   (x);
      exact_N    = cmasa_eval_1d_exact_rho_N (x);
      exact_Ntwo = cmasa_eval_1d_exact_rho_N2(x);

      // get comparison solution
      ufield2   = SourceQ_rho_u  (x,R_N,rho_N_0,rho_N_x,a_rho_N_x,
				  rho_N2_0,rho_N2_x,a_rho_N2_x,L,
				  u_0,u_x,a_ux,T_0,T_x,a_Tx);

      efield2   = SourceQ_rho_e  (x,R_N,R_N2,h0_N,h0_N2,theta_v_N2,
				  rho_N_0,rho_N_x,a_rho_N_x,rho_N2_0,
				  rho_N2_x,a_rho_N2_x,L,u_0,u_x,a_ux,
				  T_0,T_x,a_Tx);

      N2        = SourceQ_rho_N  (x,M_N,h0_N,h0_N2,Cf1_N,Cf1_N2,
				  etaf1_N,etaf1_N2,Ea_N,Ea_N2,
				  //Function_to_Calculate_K,
				  &temp_function,
				  R_N,R_N2,theta_v_N2,
				  rho_N_0,rho_N_x,a_rho_N_x,rho_N2_0,
				  rho_N2_x,a_rho_N2_x,L,u_0,u_x,a_ux,
				  T_0,T_x,a_Tx,R);


      Ntwo2     = SourceQ_rho_N2 (x,
				  M_N,
				  h0_N,
				  h0_N2,
				  Cf1_N,
				  Cf1_N2,
				  etaf1_N,
				  etaf1_N2,
				  Ea_N,
				  Ea_N2,
				  //Function_to_Calculate_K,
				  &temp_function,
				  K,
				  R_N,
				  R_N2,
				  theta_v_N2,
				  rho_N_0,
				  rho_N_x,
				  a_rho_N_x,
				  rho_N2_0,
				  rho_N2_x,
				  a_rho_N2_x,
				  L,
				  u_0,
				  u_x,
				  a_ux,
				  T_0,
				  T_x,
				  a_Tx,
				  R);

      exact_t2    = anQ_t      (x,T_0,T_x,a_Tx,L);
      exact_u2    = anQ_u      (x,u_0,u_x,a_ux,L);
      exact_rho2  = anQ_rho    (x,rho_N_0,rho_N_x,a_rho_N_x,L,
				rho_N2_0,rho_N2_x,a_rho_N2_x);

      exact_N2    = anQ_rho_N  (x,rho_N_0,rho_N_x,a_rho_N_x,L);
      exact_Ntwo2 = anQ_rho_N2 (x,rho_N2_0,rho_N2_x,a_rho_N2_x,L);

      // test the result is roughly zero
      ufield3 = fabs(ufield-ufield2);
      efield3 = fabs(efield-efield2);
      N3      = fabs(N-N2);
      Ntwo3   = fabs(Ntwo-Ntwo2);

      exact_t3    = fabs(exact_t2-exact_t);
      exact_u3    = fabs(exact_u2-exact_u);
      exact_rho3  = fabs(exact_rho2-exact_rho);
      exact_N3    = fabs(exact_N2-exact_N);
      exact_Ntwo3 = fabs(exact_Ntwo2-exact_Ntwo);

      // check threshold has not been exceeded
      threshcheck(ufield3,threshold);
      threshcheck(exact_u3   ,threshold);
      threshcheck(efield3,threshold);
      threshcheck(exact_rho3 ,threshold);
      threshcheck(exact_N3   ,threshold);
      threshcheck(exact_Ntwo3,threshold);
      threshcheck(exact_t3   ,threshold);
      threshcheck(  Ntwo3,threshold);
      threshcheck(     N3,threshold);

    } // done w/ spatial interations
  
  // steady as she goes
  return 0;

}
