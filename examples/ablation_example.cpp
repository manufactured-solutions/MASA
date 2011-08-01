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
// $Author: nick $
// $Id: 
//
// ablation_example.cpp:
// this is an example of the API used for calling the 1d ablation 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>

using namespace MASA;
using namespace std;

typedef double Scalar;

template<typename Scalar>
Scalar funct(Scalar T)
{
  // hackish functional here
  // This is an eyeballed fit (focusing on the 5000K-6000K range) 
  // for the equilibrium constant for N2->N+N dissociation
  Scalar K = exp(4+(T-6000)/500);
  return K;
}

/* The following are the variables needed by the ablation solution:

    Scalar R; 
    Scalar k; 

    Scalar u_0;
    Scalar u_x;
    Scalar a_ux;
    Scalar Gamma;
    Scalar mu;
    Scalar L;
    
    Scalar T_0;
    Scalar T_x;
    Scalar a_Tx;

    Scalar W_C;
    Scalar W_C3;    

    Scalar rho_N_0;
    Scalar rho_N_x;
    Scalar rho_N2_0;
    Scalar rho_N2_x;
    Scalar a_rho_N2_x;

    Scalar rho_C3_0;
    Scalar rho_C3_x;
    Scalar a_rho_C3_x;

    Scalar rho_C_0;
    Scalar rho_C_x;
    Scalar a_rho_C_x;

    Scalar rho_an_C3;

    Scalar k_B;
    Scalar beta_C3;
    Scalar A_C3Enc;
    Scalar D_C;
    Scalar D_C3;
    Scalar m_C3;
    Scalar E_aC3nc;

    Scalar qr;
    Scalar alpha;
    Scalar a_rho_N_x;
    Scalar sigma;
    Scalar epsilon;
*/

int main()
{
  // declarations
  Scalar tempx,tempy,tempz;

  Scalar ufield;
  Scalar efield;
  Scalar cfield;
  Scalar c3field;
  Scalar rho_cfield;
  Scalar rho_c3field;
  Scalar boundary;

  Scalar e_rho_field;

  Scalar exact_u;
  Scalar exact_t;
  Scalar exact_rho;
  Scalar exact_c;
  Scalar exact_c3;

  //problem size
  Scalar lx;
  Scalar dx;
  int nx;

  // initialize
  nx = 10;  // number of points
  lx=1;     // length

  // error handling
  int err=0;

  dx=(Scalar)lx/(Scalar)nx;

  // initialize the problem
  err += masa_init<Scalar>("ablation-example","navierstokes_ablation_1d_steady");

  // test that all variables have been initialized
  err += masa_sanity_check<Scalar>();

  // evaluate source terms over the domain (0<x<1)
  for(int i=0;i<nx;i++)
    {  
      tempx=i*dx;
      tempy=0;
      tempz=0;

      /* evaluate source terms */
      ufield = masa_eval_source_rho_u<Scalar>       (tempx);
      e_rho_field = masa_eval_source_rho_e<Scalar>  (tempx);     
      efield  = masa_eval_source_e <Scalar>         (tempx,&funct);
      cfield  = masa_eval_source_C <Scalar>         (tempx);
      c3field = masa_eval_source_C3<Scalar>         (tempx);
      rho_cfield  = masa_eval_source_rho_C <Scalar> (tempx);
      rho_c3field = masa_eval_source_rho_C3<Scalar> (tempx);
      boundary   = masa_eval_source_boundary<Scalar> (tempx);
	
      /* evaluate manufactured analytical solution */
      exact_u   = masa_eval_exact_u  <Scalar>   (tempx);
      exact_t   = masa_eval_exact_t  <Scalar>   (tempx);
      exact_rho = masa_eval_exact_rho<Scalar>   (tempx);
      exact_c   = masa_eval_exact_rho_C<Scalar> (tempx);
      exact_c3  = masa_eval_exact_rho_C3<Scalar>(tempx);
    
      masa_test_default(ufield);
      masa_test_default(efield);
      masa_test_default(cfield);
      masa_test_default(c3field);
      masa_test_default(rho_c3field);
      masa_test_default(rho_cfield);
      masa_test_default(rho_c3field);
      masa_test_default(boundary);
      masa_test_default(e_rho_field);

      masa_test_default(exact_u);
      masa_test_default(exact_t);
      masa_test_default(exact_rho);
      masa_test_default(exact_c);
      masa_test_default(exact_c3);

    }

  return err;

}// end program
