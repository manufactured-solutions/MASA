// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
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
//
// tests.h: helper routines for testing
//
// $Id: masa.h.in 19231 2011-03-30 01:16:36Z nick $
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h> 
#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

double nancheck(double x)
{

  if(isnan(x))
    {
      printf("MASA REGRESSION FAILURE:: nan found!\n");
      exit(1);
    }
  
  return 0;

}

double threshcheck(double x, double thresh)
{
  if(isnan(x))
    {
      printf("MASA REGRESSION FAILURE:: nan found!\n");
      exit(1);
    }
  
  if(x > thresh)
    {
      printf("\nMASA REGRESSION TEST FAILED!\n");
      printf("Difference of: %g\n",x);
      printf("Exceeded Threshold: %g\n",thresh);
      exit(1);
    }
  return 0;  
}



double test_default(double x)
{
  double MASA_VAR_DEFAULT = -12345.67;
  //double uninit = -1.33;
  double thresh = 5 * 1e-15;

  if(isnan(x))
    {
      printf("MASA REGRESSION FAILURE:: nan found!\n");
      exit(1);
    }
  
  if((fabs(x-MASA_VAR_DEFAULT)/MASA_VAR_DEFAULT) > thresh)
    {
      printf("\nMASA REGRESSION TEST FAILED!\n");
      printf("Not set to default! %g",x);
      exit(1);
    }
  return 0;  
}



#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
#include <limits>
#include <iostream>

template<typename Scalar>
Scalar nancheck(Scalar x)
{

  if(isnan(x))
    {
      std::cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
}

template<typename Scalar>
Scalar threshcheck(Scalar x, Scalar thresh)
{
  //Scalar MASA_VAR_DEFAULT = -12345.67;
  //Scalar uninit = -1.33;
  //Scalar thresh = 5 * std::numeric_limits<Scalar>::epsilon();

  if(isnan(x))
    {
      std::cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
  
  if(x > thresh)
    {
      std::cout << "\nMASA REGRESSION TEST FAILED!\n";
      std::cout << "Difference of: " << x << std::endl;
      std::cout << "Exceeded Threshold: " << thresh << std::endl;
      exit(1);
    }
  return 0;  
}

template<typename Scalar>
Scalar test_default(Scalar x)
{
  Scalar MASA_VAR_DEFAULT = -12345.67;
  //Scalar uninit = -1.33;
  Scalar thresh = 5 * std::numeric_limits<Scalar>::epsilon();

  if(isnan(x))
    {
      std::cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
  
  if((fabs(x-MASA_VAR_DEFAULT)/MASA_VAR_DEFAULT) > thresh)
    {
      std::cout << "\nMASA REGRESSION TEST FAILED!\n";
      std::cout << "Exceeded Threshold by: " << x << std::endl;
      exit(1);
    }
  return 0;  
}

template<typename Scalar>
Scalar threshcheck(Scalar x)
{
  //Scalar MASA_VAR_DEFAULT = -12345.67;
  //Scalar uninit = -1.33;
  Scalar thresh = 5 * std::numeric_limits<Scalar>::epsilon();

  if(isnan(x))
    {
      std::cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
  
  if(x > thresh)
    {
      std::cout << "\nMASA REGRESSION TEST FAILED!\n";
      std::cout << "Exceeded Threshold by: " << x << std::endl;
      exit(1);
    }
  return 0;  
}

template<typename Scalar>
Scalar tester(Scalar a)
{
  
  a = 4.4;
  return a;

}

template<typename Scalar>
Scalar test_grad(Scalar derr)
{
  //Scalar MASA_VAR_DEFAULT = -12345.67;
  //Scalar uninit = -1.33;
  Scalar thresh = 5 * std::numeric_limits<Scalar>::epsilon();

  if(fabs(derr+1) > thresh)
    {
      std::cout << "MASA :: gradient error condition failed!\n";
      exit(1);
    }
  
}

template<typename Scalar>
Scalar temp_function(Scalar T)
{
  // hackish functional here
  // This is an eyeballed fit (focusing on the 5000K-6000K range) 
  // for the equilibrium constant for N2->N+N dissociation
  Scalar K = exp(4+(T-6000)/500);
  return K;
};

template<typename Scalar> int run_regression();

// heat sv
template<typename Scalar>
Scalar SourceQ_t_1d(
  Scalar x,
  Scalar A_x,
  Scalar k_0,
  Scalar k_1,
  Scalar k_2);

template<typename Scalar>
Scalar SourceQ_t_2d (
  Scalar x,
  Scalar y,
  Scalar A_x,
  Scalar B_y,
  Scalar k_0,
  Scalar k_1,
  Scalar k_2);

template<typename Scalar>
Scalar SourceQ_t_3d (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar A_x,
  Scalar B_y,
  Scalar C_z,
  Scalar k_0,
  Scalar k_1,
  Scalar k_2);

// heat sc
template<typename Scalar> Scalar SourceQ_t_1d(Scalar, Scalar, Scalar);
template<typename Scalar> Scalar Source_t_1d_exact(Scalar A_x,Scalar x);
template<typename Scalar> Scalar SourceQ_t_2d (Scalar x,Scalar y,Scalar A_x,Scalar B_y,Scalar k_0);
template<typename Scalar> Scalar SourceQ_t_3d (Scalar x,Scalar y,Scalar z,Scalar A,Scalar B,Scalar C,Scalar k_0);

// heat uc
template<typename Scalar>Scalar SourceQ_t_1d(
  Scalar x,
  Scalar t,
  Scalar A_x,
  Scalar A_t,
  Scalar D_t,
  Scalar k_0,
  Scalar cp_0,
  Scalar rho);

template<typename Scalar>
Scalar SourceQ_t_2d (
  Scalar x,
  Scalar y,
  Scalar t,
  Scalar A_x,
  Scalar A_t,
  Scalar B_y,
  Scalar B_t,
  Scalar D_t,
  Scalar rho,
  Scalar k_0,
  Scalar cp_0);

template<typename Scalar>
Scalar SourceQ_t_3d (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar t,
  Scalar A_x,
  Scalar A_t,
  Scalar B_y,
  Scalar B_t,
  Scalar C_z,
  Scalar C_t,
  Scalar D_t,
  Scalar k_0,
  Scalar cp_0,
  Scalar rho);


// heat uv
template<typename Scalar>
Scalar SourceQ_t_1d (
  Scalar x,
  Scalar t,
  Scalar A_x,
  Scalar A_t,
  Scalar D_t,
  Scalar rho,
  Scalar k_0,
  Scalar k_1,
  Scalar k_2,
  Scalar cp_0,
  Scalar cp_1,
  Scalar cp_2);

template<typename Scalar>
Scalar SourceQ_t_2d (
  Scalar x,
  Scalar y,
  Scalar t,
  Scalar A_x,
  Scalar A_t,
  Scalar B_y,
  Scalar B_t,
  Scalar D_t,
  Scalar rho,
  Scalar k_0,
  Scalar k_1,
  Scalar k_2,
  Scalar cp_0,
  Scalar cp_1,
  Scalar cp_2);

template<typename Scalar>
Scalar SourceQ_t_3d (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar t,
  Scalar A_x,
  Scalar A_t,
  Scalar B_y,
  Scalar B_t,
  Scalar C_z,
  Scalar C_t,
  Scalar D_t,
  Scalar rho,
  Scalar k_0,
  Scalar k_1,
  Scalar k_2,
  Scalar cp_0,
  Scalar cp_1,
  Scalar cp_2);

// euler 3d
template<typename Scalar>
Scalar anQ_p (Scalar x,Scalar y,Scalar z,Scalar p_0,Scalar p_x,Scalar p_y,Scalar p_z,Scalar a_px,Scalar a_py,Scalar a_pz,Scalar L);

template<typename Scalar>
Scalar anQ_u (Scalar x,Scalar y,Scalar z,Scalar u_0,Scalar u_x,Scalar u_y,Scalar u_z,Scalar a_ux,Scalar a_uy,Scalar a_uz,Scalar L);

template<typename Scalar>
Scalar anQ_v (Scalar x,Scalar y,Scalar z,Scalar v_0,Scalar v_x,Scalar v_y,Scalar v_z,Scalar a_vx,Scalar a_vy,Scalar a_vz,Scalar L);

template<typename Scalar>
Scalar anQ_w (Scalar x,Scalar y,Scalar z,Scalar w_0,Scalar w_x,Scalar w_y,Scalar w_z,Scalar a_wx,Scalar a_wy,Scalar a_wz,Scalar L);

template<typename Scalar>
Scalar anQ_rho (Scalar x,Scalar y,Scalar z,Scalar rho_0,Scalar rho_x,Scalar rho_y,Scalar rho_z,Scalar a_rhox,Scalar a_rhoy,Scalar a_rhoz,Scalar L);

template<typename Scalar>
Scalar SourceQ_e (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar mu,
  Scalar Gamma,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_u ( // should be 38
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_v ( // 38
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_w ( // 38
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_rho( // 39
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar mu,
  Scalar L);


// euler 2d
template<typename Scalar>
Scalar anQ_p (Scalar x,Scalar y,Scalar p_0,Scalar p_x,Scalar p_y,Scalar a_px,Scalar a_py,Scalar L);

template<typename Scalar>
Scalar anQ_u (Scalar x,Scalar y,Scalar u_0,Scalar u_x,Scalar u_y,Scalar a_ux,Scalar a_uy,Scalar L);

template<typename Scalar> 
Scalar anQ_v (Scalar x,Scalar y,Scalar v_0,Scalar v_x,Scalar v_y,Scalar a_vx,Scalar a_vy,Scalar L);

template<typename Scalar>
Scalar anQ_rho (Scalar x,Scalar y,Scalar rho_0,Scalar rho_x,Scalar rho_y,Scalar a_rhox,Scalar a_rhoy,Scalar L);

template<typename Scalar>
Scalar SourceQ_e (
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar Gamma,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_u ( // 23 variables
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_x,
  Scalar a_px,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_v (
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_y,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_rho(
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar L);

//euler 1d
template<typename Scalar>
Scalar anQ_u (Scalar x,Scalar u_0,Scalar u_x,Scalar a_ux,Scalar L);

template<typename Scalar>
Scalar anQ_p (Scalar x,Scalar p_0,Scalar p_x,Scalar a_px,Scalar L);

template<typename Scalar>
Scalar anQ_rho (Scalar x,Scalar rho_0,Scalar rho_x,Scalar a_rhox,Scalar L);

template<typename Scalar>
Scalar SourceQ_e ( // 12
  Scalar x,
  Scalar u_0,
  Scalar u_x,
  Scalar rho_0,
  Scalar rho_x,
  Scalar p_0,
  Scalar p_x,
  Scalar a_px,
  Scalar a_rhox,
  Scalar a_ux,
  Scalar Gamma,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_u ( // should be 10
  Scalar x,
  Scalar u_0,
  Scalar u_x,
  Scalar rho_0,
  Scalar rho_x,
  Scalar p_x,
  Scalar a_px,
  Scalar a_rhox,
  Scalar a_ux,
  Scalar L);


template<typename Scalar>
Scalar SourceQ_rho ( 
  Scalar x,
  Scalar u_0,
  Scalar u_x,
  Scalar rho_0,
  Scalar rho_x,
  Scalar a_rhox,
  Scalar a_ux,
  Scalar L);

//euler transient

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
  Scalar Gamma);

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
  Scalar a_px,
  Scalar a_rhox,
  Scalar a_rhot,
  Scalar a_ux,
  Scalar a_ut,
  Scalar L);

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
  Scalar a_rhox,
  Scalar a_rhot,
  Scalar a_ux,
  Scalar a_ut,
  Scalar L);

// euler chemistry

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
		     Scalar a_Tx);

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
		     Scalar a_Tx);

template <typename Scalar>
Scalar SourceQ_rho_N(Scalar x,
		     Scalar M_N,
		     Scalar Cf1_N,
		     Scalar Cf1_N2,
		     Scalar etaf1_N,
		     Scalar etaf1_N2,
		     Scalar Ea_N,
		     Scalar Ea_N2,
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
		     Scalar R);

template <typename Scalar>
Scalar SourceQ_rho_N2(Scalar x,
		      Scalar M_N,
		      Scalar Cf1_N,
		      Scalar Cf1_N2,
		      Scalar etaf1_N,
		      Scalar etaf1_N2,
		      Scalar Ea_N,
		      Scalar Ea_N2,
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
		      Scalar R);

template <typename Scalar>
Scalar anQ_t(Scalar x,Scalar T_0,Scalar T_x,Scalar a_Tx,Scalar L);

template <typename Scalar>
Scalar anQ_rho_N(Scalar x,Scalar rho_N_0,Scalar rho_N_x,Scalar a_rho_N_x,Scalar L);

template <typename Scalar>
Scalar anQ_rho_N2(Scalar x,Scalar rho_N2_0,Scalar rho_N2_x,Scalar a_rho_N2_x,Scalar L);


// euler3d
template<typename Scalar>
Scalar SourceQ_e (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar mu,
  Scalar Gamma,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_u ( // should be 38
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_x,
  Scalar a_px,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_v ( // 38
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_y,
  Scalar p_z,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_w ( // 38
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_z,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar L);

template<typename Scalar>
Scalar SourceQ_rho( // 39
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_y,
  Scalar p_z,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar mu,
  Scalar L);

// cns 2d

template<typename Scalar>
Scalar SourceQ_e (
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar Gamma,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k);

template<typename Scalar>
Scalar SourceQ_u (
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k);

template<typename Scalar>
Scalar SourceQ_v (
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k);

template<typename Scalar>
Scalar SourceQ_rho( 
  Scalar x,
  Scalar y,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar a_px,
  Scalar a_py,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_vx,
  Scalar a_vy,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k);

// cns 3d

template<typename Scalar>
Scalar SourceQ_e (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar mu,
  Scalar Gamma,
  Scalar L,
  Scalar R,
  Scalar k);

template<typename Scalar>
Scalar SourceQ_u (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k);

template<typename Scalar>
Scalar SourceQ_v (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k);

template<typename Scalar>
Scalar SourceQ_w (
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k);

template<typename Scalar>
Scalar SourceQ_rho(
  Scalar x,
  Scalar y,
  Scalar z,
  Scalar u_0,
  Scalar u_x,
  Scalar u_y,
  Scalar u_z,
  Scalar v_0,
  Scalar v_x,
  Scalar v_y,
  Scalar v_z,
  Scalar w_0,
  Scalar w_x,
  Scalar w_y,
  Scalar w_z,
  Scalar rho_0,
  Scalar rho_x,
  Scalar rho_y,
  Scalar rho_z,
  Scalar p_0,
  Scalar p_x,
  Scalar p_y,
  Scalar p_z,
  Scalar a_px,
  Scalar a_py,
  Scalar a_pz,
  Scalar a_rhox,
  Scalar a_rhoy,
  Scalar a_rhoz,
  Scalar a_ux,
  Scalar a_uy,
  Scalar a_uz,
  Scalar a_vx,
  Scalar a_vy,
  Scalar a_vz,
  Scalar a_wx,
  Scalar a_wy,
  Scalar a_wz,
  Scalar mu,
  Scalar L,
  Scalar R,
  Scalar k);


// axi euler
template<typename Scalar>
Scalar anQ_p(Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma);

template<typename Scalar>
Scalar anQ_u (Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma);

template<typename Scalar>
Scalar anQ_w (Scalar r,Scalar z,Scalar w_0,Scalar w_1,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L);

template<typename Scalar>
Scalar anQ_rho (Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma);

template<typename Scalar>
Scalar SourceQ_e(Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma);

template<typename Scalar>
Scalar SourceQ_w(Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma);

template<typename Scalar>
Scalar SourceQ_rho(Scalar r,Scalar z,Scalar p_0,Scalar p_1,Scalar rho_0,Scalar rho_1,Scalar u_1,Scalar w_0,Scalar w_1,Scalar a_pr,Scalar a_pz,Scalar a_rhor,Scalar a_rhoz,Scalar a_ur,Scalar a_uz,Scalar a_wr,Scalar a_wz,Scalar pi,Scalar L,Scalar Gamma);





#endif // __cplusplus
