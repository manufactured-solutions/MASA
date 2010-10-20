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
// euler1d.cpp :program that tests masa against known source term
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <math.h>

#include <iostream>
#include <stdlib.h>

using namespace MASA;
using namespace std;

typedef double FP;

const FP pi = acos(-1);
const FP threshold = 1.0e-15; // should be small enough to catch any obvious problems

FP nancheck(FP x)
{
  if(isnan(x))
    {
      cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
  return 1;
}

FP anQ_p (FP x,FP p_0,FP p_x,FP a_px,FP L)
{
  FP p_an = p_0 + p_x * cos(a_px * pi * x / L);
  return p_an;
}
  
FP anQ_u (FP x,FP u_0,FP u_x,FP a_ux,FP L)
{
  FP u_an = u_0 + u_x * sin(a_ux * pi * x / L);
  return u_an;
} 
 
FP anQ_rho (FP x,FP rho_0,FP rho_x,FP a_rhox,FP L)
{ 
  FP rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L);
  return rho_an;
}

FP SourceQ_e ( // 12
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP);

FP SourceQ_u ( // should be 10
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP);

FP SourceQ_rho ( // 10
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP,
  FP);

int main()
{
  //variables 
  FP u_0;
  FP u_x;
  FP rho_0;
  FP rho_x;
  FP p_0;
  FP p_x;
  FP a_px;
  FP a_rhox;
  FP a_ux;
  FP Gamma;
  FP mu;
  FP L;

  // parameters
  FP x;

  //problem size
  int nx = 200;  // number of points
  int lx=10;     // length
  FP dx=FP(lx)/FP(nx);

  // solutions
  FP ufield,ufield2,ufield3;
  FP efield,efield2,efield3;
  FP rho,rho2,rho3;
  FP gradx,grady,gradz,gradp,gradrho;

  FP u_an,u_an2,u_an3;
  FP p_an,p_an2,p_an3;
  FP rho_an,rho_an2,rho_an3;

  // initalize
  masa_init<FP>("euler-test","euler_1d");

  // initialize the default parameters
  masa_init_param<FP>();

  // get defaults for comparison to source terms
  // get vars
  u_0 = masa_get_param<FP>("u_0");
  u_x = masa_get_param<FP>("u_x");
  rho_0 = masa_get_param<FP>("rho_0");
  rho_x = masa_get_param<FP>("rho_x");

  p_0 = masa_get_param<FP>("p_0");
  p_x = masa_get_param<FP>("p_x");

  a_px = masa_get_param<FP>("a_px");
  a_rhox = masa_get_param<FP>("a_rhox");
  a_ux = masa_get_param<FP>("a_ux");

  Gamma = masa_get_param<FP>("Gamma");
  mu    = masa_get_param<FP>("mu");
  L     = masa_get_param<FP>("L");

  // check that all terms have been initialized
  masa_sanity_check<FP>();

  // evaluate source terms (1D)
  for(int i=0;i<nx;i++)
    {
      x=i*dx;
      	
      //evalulate source terms
      ufield = masa_eval_u_source  <FP>(x);
      efield = masa_eval_e_source  <FP>(x);
      rho    = masa_eval_rho_source<FP>(x);
      
      //evaluate analytical terms
      u_an = masa_eval_u_an        <FP>(x);
      p_an = masa_eval_p_an        <FP>(x);
      rho_an = masa_eval_rho_an    <FP>(x);

      // eval gradient terms
      gradx   = masa_eval_1d_grad_u  (x);
      gradp   = masa_eval_1d_grad_p  (x);
      gradrho = masa_eval_1d_grad_rho(x);

      // get fundamental source term solution
      ufield2   = SourceQ_u  (x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L);
      rho2      = SourceQ_rho(x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,L);
      efield2   = SourceQ_e  (x,u_0,u_x,rho_0,rho_x,p_0,p_x,a_px,a_rhox,a_ux,Gamma,mu,L);
  
      u_an2   = anQ_u   (x,u_0,u_x,a_ux,L);
      rho_an2 = anQ_rho (x,rho_0,rho_x,a_rhox,L);
      p_an2   = anQ_p   (x,p_0,p_x,a_px,L);

      // test the result is roughly zero
      // choose between abs and rel error
#ifdef MASA_STRICT_REGRESSION

      ufield3 = fabs(ufield-ufield2);
      efield3 = fabs(efield-efield2);
      rho3    = fabs(rho-rho2);

      u_an3   = fabs(u_an-u_an2);
      rho_an3 = fabs(rho_an-rho_an2);
      p_an3   = fabs(p_an-p_an2);

#else

      ufield3 = fabs(ufield-ufield2)/fabs(ufield2);
      efield3 = fabs(efield-efield2)/fabs(efield2);
      rho3    = fabs(rho-rho2)/fabs(rho2);

      u_an3   = fabs(u_an-u_an2)/fabs(u_an2);
      rho_an3 = fabs(rho_an-rho_an2)/fabs(rho_an2);
      p_an3   = fabs(p_an-p_an2)/fabs(p_an2);

#endif

      nancheck(ufield3);
      nancheck(efield3);
      nancheck(rho3);
      
      nancheck(u_an3);
      nancheck(rho_an3);
      nancheck(p_an3);
      
      if(ufield3 > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "U Field Source Term\n";
	  cout << "Exceeded Threshold by: " << ufield << endl;
	  cout << x << " " << endl;
	  exit(1);
	}

      if(u_an3 > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "U Field Analytical Term\n";
	  cout << "Exceeded Threshold by: " << u_an << endl;
	  cout << x << " " << endl;
	  exit(1);
	}
      
      if(efield3 > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "Threshold Exceeded: " << efield3 << endl;
	  cout << "MASA:               " <<  efield << endl;
	  cout << "Maple:              " <<  efield2 << endl;
	  cout << x << endl;
	  exit(1);
	}

      if(p_an3 > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "P Field Analytical Term\n";
	  cout << "Exceeded Threshold by: " << p_an << endl;
	  cout << x << endl;
	  exit(1);
	}
      
      if(rho3 > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "RHO Source Term\n";
	  cout << "Exceeded Threshold by: " << rho << endl;
	  cout << x << endl;
	  exit(1);
	}
      
      if(rho_an3 > threshold)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
	  cout << "RHO Analytical Term\n";
	  cout << "Exceeded Threshold by: " << rho_an << endl;
	  cout << x << endl;
	  exit(1);
	}

      // adding a new error check: ensure physical results are coming out!
      /*
	if(0 > rho)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
   	  cout << "Initial Variables are returning non-physical results!\n";
	  cout << "RHO\n";
	  exit(1);
	}

      if(0 > rho_an)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
   	  cout << "Initial Variables are returning non-physical results!\n";
	  cout << "RHO analytical\n";
	  exit(1);
	}

      if(0 > p_an)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
   	  cout << "Initial Variables are returning non-physical results!\n";
	  cout << "Pressure is negative!\n";
	  exit(1);
	}

      if(0 > efield)
	{
	  cout << "\nMASA REGRESSION TEST FAILED: Euler-1d\n";
   	  cout << "Initial Variables are returning non-physical results!\n";
	  cout << "Energy is negative!\n";
	  exit(1);
	}
      */
    } // done interating 

  // tests passed
  return 0;
}
