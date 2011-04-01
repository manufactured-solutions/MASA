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

template<typename Scalar>
Scalar temp_function(Scalar T)
{
  // hackish functional here
  // This is an eyeballed fit (focusing on the 5000K-6000K range) 
  // for the equilibrium constant for N2->N+N dissociation
  Scalar K = exp(4+(T-6000)/500);
  return K;
}

template <typename Scalar>
Scalar SourceQ_rho_u(Scalar x, 
		     Scalar R_N, 
		     Scalar rho_N_0, 
		     Scalar rho_N_x,
		     Scalar a_rho_N_x, 
		     Scalar rho_N2_0, 
		     Scalar rho_N2_x,
		     Scalar a_rho_N2_x,
		     Scalar L,
		     Scalar u_0,
		     Scalar u_x,
		     Scalar a_ux,
		     Scalar T_0,
		     Scalar T_x,
		     Scalar a_Tx)
{
  Scalar Q_u;
  Scalar RHO;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar U;
  Scalar T;

  Scalar pi = acos(-1);

  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  RHO = RHO_N + RHO_N2;
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  T = T_0 + T_x * cos(a_Tx * pi * x / L);

  Q_u = a_rho_N_x * pi * rho_N_x * U * U * cos(a_rho_N_x * pi * x / L) / L - a_rho_N2_x * pi * rho_N2_x * U * U * sin(a_rho_N2_x * pi * x / L) / L - a_Tx * pi * T_x * R_N * RHO_N * sin(a_Tx * pi * x / L) / L - a_Tx * pi * T_x * R_N * RHO_N2 * sin(a_Tx * pi * x / L) / L / 0.2e1 + 0.2e1 * a_ux * pi * u_x * RHO * U * cos(a_ux * pi * x / L) / L - (-0.2e1 * a_rho_N_x * rho_N_x * cos(a_rho_N_x * pi * x / L) + a_rho_N2_x * rho_N2_x * sin(a_rho_N2_x * pi * x / L)) * pi * R_N * T / L / 0.2e1;
  return(Q_u);
}

template <typename Scalar>
Scalar SourceQ_rho_e(Scalar x,
		     Scalar R_N,
		     Scalar R_N2,
		     Scalar h0_N,
		     Scalar h0_N2,
		     Scalar theta_v_N2,
		     Scalar rho_N_0, 
		     Scalar rho_N_x,
		     Scalar a_rho_N_x, 
		     Scalar rho_N2_0, 
		     Scalar rho_N2_x,
		     Scalar a_rho_N2_x,
		     Scalar L,
		     Scalar u_0,
		     Scalar u_x,
		     Scalar a_ux,
		     Scalar T_0,
		     Scalar T_x,
		     Scalar a_Tx)
{
  Scalar Q_e;
  Scalar RHO;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar U;
  Scalar T;
  Scalar alpha;
  Scalar E_vib_N2;

  Scalar pi = acos(-1);

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

template <typename Scalar>
Scalar SourceQ_rho_N(Scalar x,
		     Scalar M_N,
		     Scalar Cf1_N,
		     Scalar Cf1_N2,
		     Scalar etaf1_N,
		     Scalar etaf1_N2,
		     Scalar Ea_N,
		     Scalar Ea_N2,
		     Scalar R_N,
		     Scalar R_N2,
		     Scalar theta_v_N2,
		     Scalar rho_N_0, 
		     Scalar rho_N_x,
		     Scalar a_rho_N_x, 
		     Scalar rho_N2_0, 
		     Scalar rho_N2_x,
		     Scalar a_rho_N2_x,
		     Scalar L,
		     Scalar u_0,
		     Scalar u_x,
		     Scalar a_ux,
		     Scalar T_0,
		     Scalar T_x,
		     Scalar a_Tx,
		     Scalar R)

{
  Scalar Q_rho_N;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar U;
  Scalar T;
  Scalar kf1_N;
  Scalar kf1_N2;
  Scalar K_eq;

  Scalar pi = acos(-1);

  T = T_0 + T_x * cos(a_Tx * pi * x / L);
  K_eq = temp_function(T);

  RHO_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L);
  RHO_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  U = u_0 + u_x * sin(a_ux * pi * x / L);
  kf1_N = Cf1_N * pow(T, etaf1_N) * exp(-Ea_N / R / T);
  kf1_N2 = Cf1_N2 * pow(T, etaf1_N2) * exp(-Ea_N2 / R / T);

  Q_rho_N = a_rho_N_x * pi * rho_N_x * U * cos(a_rho_N_x * pi * x / L) / L + a_ux * pi * u_x * RHO_N * cos(a_ux * pi * x / L) / L + (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N * RHO_N * pow(M_N, -0.2e1) / K_eq - (0.2e1 * kf1_N * RHO_N + kf1_N2 * RHO_N2) * RHO_N2 / M_N / 0.2e1;

  return(Q_rho_N);
}

template <typename Scalar>
Scalar SourceQ_rho_N2(Scalar x,
		      Scalar M_N,
		      Scalar Cf1_N,
		      Scalar Cf1_N2,
		      Scalar etaf1_N,
		      Scalar etaf1_N2,
		      Scalar Ea_N,
		      Scalar Ea_N2,
		      //Scalar (*in_func)(Scalar),
		      Scalar K,
		      Scalar R_N,
		      Scalar R_N2,
		      Scalar theta_v_N2,
		      Scalar rho_N_0, 
		      Scalar rho_N_x,
		      Scalar a_rho_N_x, 
		      Scalar rho_N2_0, 
		      Scalar rho_N2_x,
		      Scalar a_rho_N2_x,
		      Scalar L,
		      Scalar u_0,
		      Scalar u_x,
		      Scalar a_ux,
		      Scalar T_0,
		      Scalar T_x,
		      Scalar a_Tx,
		      Scalar R)
{
  
  Scalar Q_rho_N2;
  Scalar RHO_N;
  Scalar RHO_N2;
  Scalar U;
  Scalar T;
  Scalar kf1_N;
  Scalar kf1_N2;
  Scalar K_eq;

  Scalar pi = acos(-1);
  T = T_0 + T_x * cos(a_Tx * pi * x / L);
  K_eq = temp_function(T);

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

template <typename Scalar>
Scalar anQ_t(Scalar x,Scalar T_0,Scalar T_x,Scalar a_Tx,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar T_an = T_0 + T_x * cos(a_Tx * pi * x / L);
  return T_an;
}

template <typename Scalar>
Scalar anQ_u(Scalar x,Scalar u_0,Scalar u_x,Scalar a_ux,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar u_an = u_0 + u_x * sin(a_ux * pi * x / L);
  return u_an;
}

template <typename Scalar>
Scalar anQ_rho(Scalar x,Scalar rho_N_0,Scalar rho_N_x,Scalar a_rho_N_x,
	       Scalar L,Scalar rho_N2_0,Scalar rho_N2_x,Scalar a_rho_N2_x)
{  
  Scalar pi = acos(-1);
  Scalar rho_an = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L) + rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  return rho_an;
}

template <typename Scalar>
Scalar anQ_rho_N(Scalar x,Scalar rho_N_0,Scalar rho_N_x,Scalar a_rho_N_x,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar rho_an_N = rho_N_0 + rho_N_x * sin(a_rho_N_x * pi * x / L);
  return rho_an_N;
}

template <typename Scalar>
Scalar anQ_rho_N2(Scalar x,Scalar rho_N2_0,Scalar rho_N2_x,Scalar a_rho_N2_x,Scalar L)
{
  Scalar pi = acos(-1);
  Scalar rho_an_N2 = rho_N2_0 + rho_N2_x * cos(a_rho_N2_x * pi * x / L);
  return rho_an_N2;
}

// ----------------------------------------
//   Regresssion
// ----------------------------------------

template<typename Scalar>
int run_regression()
{

  Scalar threshold = 5 * numeric_limits<Scalar>::epsilon();

  //variables 
  Scalar u_0;
  Scalar u_x;
  Scalar a_ux;
  Scalar L;

  Scalar R;
  Scalar Cf1_N;
  Scalar Cf1_N2;
  Scalar etaf1_N;
  Scalar etaf1_N2;
  Scalar Ea_N;
  Scalar Ea_N2;
  Scalar R_N;
  Scalar R_N2;
  Scalar theta_v_N2;
  Scalar M_N;
  Scalar h0_N;
  Scalar h0_N2;
  Scalar K;

  Scalar rho_N_0;
  Scalar rho_N_x;
  Scalar a_rho_N_x;

  Scalar rho_N2_0;
  Scalar rho_N2_x;
  Scalar a_rho_N2_x;

  Scalar T_0;
  Scalar T_x;
  Scalar a_Tx;

  // parameters
  Scalar x;

  //problem size
  int nx = 200;  // number of points
  int lx=10;     // length
  Scalar dx=Scalar(lx)/Scalar(nx);

  // solutions
  Scalar ufield,ufield2,ufield3;
  Scalar efield,efield2,efield3;
  Scalar N,N2,N3;
  Scalar Ntwo,Ntwo2,Ntwo3;
  //Scalar gradx,grady,gradz,gradp,gradrho;

  Scalar exact_t,exact_t2,exact_t3;
  Scalar exact_u,exact_u2,exact_u3;
  Scalar exact_rho,exact_rho2,exact_rho3;
  Scalar exact_N,exact_N2,exact_N3;
  Scalar exact_Ntwo,exact_Ntwo2,exact_Ntwo3;

  // initalize
  masa_init<Scalar>("euler-chemistry-test","euler_chem_1d");

  // initialize the default parameters
  masa_init_param<Scalar>();

  // get defaults for comparison to source terms
  u_0  = masa_get_param<Scalar>("u_0");
  u_x  = masa_get_param<Scalar>("u_x");
  a_ux = masa_get_param<Scalar>("a_ux");
  L    = masa_get_param<Scalar>("L");
  R    = masa_get_param<Scalar>("R");

  Cf1_N   = masa_get_param<Scalar>("Cf1_N");
  Cf1_N2  = masa_get_param<Scalar>("Cf1_N2");
  etaf1_N  = masa_get_param<Scalar>("etaf1_N");
  etaf1_N2 = masa_get_param<Scalar>("etaf1_N2");

  Ea_N  = masa_get_param<Scalar>("Ea_N");
  Ea_N2 = masa_get_param<Scalar>("Ea_N2");

  R_N   = masa_get_param<Scalar>("R_N");
  R_N2  = masa_get_param<Scalar>("R_N2");

  theta_v_N2 = masa_get_param<Scalar>("theta_v_N2");
  M_N   = masa_get_param<Scalar>("M_N");
  h0_N  = masa_get_param<Scalar>("h0_N");
  h0_N2 = masa_get_param<Scalar>("h0_N2");
  K     = masa_get_param<Scalar>("K");

  rho_N_0   = masa_get_param<Scalar>("rho_N_0");
  rho_N_x   = masa_get_param<Scalar>("rho_N_x");
  a_rho_N_x = masa_get_param<Scalar>("a_rho_N_x");

  rho_N2_0   = masa_get_param<Scalar>("rho_N2_0");  
  rho_N2_x   = masa_get_param<Scalar>("rho_N2_x");
  a_rho_N2_x = masa_get_param<Scalar>("a_rho_N2_x");

  T_0  = masa_get_param<Scalar>("T_0");
  T_x  = masa_get_param<Scalar>("T_x");
  a_Tx = masa_get_param<Scalar>("a_Tx");

  // check that all terms have been initialized
  masa_sanity_check<Scalar>();

  // evaluate MMS (1D)
  for(int i=0;i<nx;i++)
    {
      x=i*dx;

      // evalulate source terms
      ufield = masa_eval_source_rho_u  <Scalar>(x);
      efield = masa_eval_source_rho_e  <Scalar>(x);
      N      = masa_eval_source_rho_N  <Scalar>(x,&temp_function);
      Ntwo   = masa_eval_source_rho_N2 <Scalar>(x,&temp_function);

      // evaluate analytical solution terms
      exact_t    = masa_eval_exact_t     <Scalar>(x);
      exact_u    = masa_eval_exact_u     <Scalar>(x);
      exact_rho  = masa_eval_exact_rho   <Scalar>(x);
      exact_N    = masa_eval_exact_rho_N <Scalar>(x);
      exact_Ntwo = masa_eval_exact_rho_N2<Scalar>(x);

      // get comparison solution
      ufield2   = SourceQ_rho_u  (x,R_N,rho_N_0,rho_N_x,a_rho_N_x,
				  rho_N2_0,rho_N2_x,a_rho_N2_x,L,
				  u_0,u_x,a_ux,T_0,T_x,a_Tx);

      efield2   = SourceQ_rho_e  (x,R_N,R_N2,h0_N,h0_N2,theta_v_N2,
				  rho_N_0,rho_N_x,a_rho_N_x,rho_N2_0,
				  rho_N2_x,a_rho_N2_x,L,u_0,u_x,a_ux,
				  T_0,T_x,a_Tx);

      N2        = SourceQ_rho_N  (x,M_N,Cf1_N,Cf1_N2,
				  etaf1_N,etaf1_N2,Ea_N,Ea_N2,
				  R_N,R_N2,theta_v_N2,
				  rho_N_0,rho_N_x,a_rho_N_x,rho_N2_0,
				  rho_N2_x,a_rho_N2_x,L,u_0,u_x,a_ux,
				  T_0,T_x,a_Tx,R);


      Ntwo2     = SourceQ_rho_N2 (x,M_N,Cf1_N,Cf1_N2,etaf1_N,etaf1_N2,
				  Ea_N,Ea_N2,
				  K,R_N,R_N2,theta_v_N2,
				  rho_N_0,rho_N_x,a_rho_N_x,rho_N2_0,rho_N2_x,
				  a_rho_N2_x,L,u_0,u_x,a_ux,T_0,T_x,a_Tx,R);

      exact_t2    = anQ_t      (x,T_0,T_x,a_Tx,L);
      exact_u2    = anQ_u      (x,u_0,u_x,a_ux,L);
      exact_rho2  = anQ_rho    (x,rho_N_0,rho_N_x,a_rho_N_x,L,
				rho_N2_0,rho_N2_x,a_rho_N2_x);

      exact_N2    = anQ_rho_N  (x,rho_N_0,rho_N_x,a_rho_N_x,L);
      exact_Ntwo2 = anQ_rho_N2 (x,rho_N2_0,rho_N2_x,a_rho_N2_x,L);

      // test the result is roughly zero
      // choose between abs and rel error
      
      // #ifdef MASA_STRICT_REGRESSION

      ufield3 = fabs(ufield-ufield2);
      efield3 = fabs(efield-efield2);
      N3      = fabs(N-N2);
      Ntwo3   = fabs(Ntwo-Ntwo2);

      exact_t3    = fabs(exact_t2-exact_t);
      exact_u3    = fabs(exact_u2-exact_u);
      exact_rho3  = fabs(exact_rho2-exact_rho);
      exact_N3    = fabs(exact_N2-exact_N);
      exact_Ntwo3 = fabs(exact_Ntwo2-exact_Ntwo);

      /* #else

      ufield3 = fabs(ufield-ufield2)/fabs(ufield2);
      efield3 = fabs(efield-efield2)/fabs(efield2);
      N3      = fabs(N-N2)/fabs(N2);
      Ntwo3   = fabs(Ntwo-Ntwo2)/fabs(Ntwo2);

      exact_t3    = fabs(exact_t2-exact_t)/fabs(exact_t3);
      exact_u3    = fabs(exact_u2-exact_u)/fabs(exact_u3);
      exact_rho3  = fabs(exact_rho2-exact_rho)/fabs(exact_rho3);
      exact_N3    = fabs(exact_N2-exact_N)/fabs(exact_N3);
      exact_Ntwo3 = fabs(exact_Ntwo2-exact_Ntwo)/fabs(exact_Ntwo3);

      #endif */

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

  return 0;

} // done with tests



int main()
{
  int err=0;

  err += run_regression<double>();
  //err += run_regression<long double>();

  return err;
}
