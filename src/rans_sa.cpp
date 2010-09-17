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
// $Id: cns.cpp 13034 2010-09-08 16:17:42Z nick $
//
// rans_sa.cpp: These are the MASA class member functions and constructors
//              For a Reynold Averaged Navier Stokes (RANS) model
//              The Spelart Alamaras
//              For Channel Flow only
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         Spelart Alamaras
 *
 *
 * -----------------------------------------------
 */ 

MASA::rans_sa::rans_sa()
{
  mmsname = "rans_sa";
  dimension=1;
  
  // value of 'peak' in eddy viscosity
  etam=.6;

  // u parameter
  a1=2;

  // eddy viscosity parameter
  b1=1;

  register_var("cb1",  &cb1);
  register_var("cb2",  &cb2);
  register_var("cv1",  &cv1);
  register_var("cw2",  &cw2);
  register_var("cw3",  &cw3);
  register_var("sigma",&sigma);
  register_var("kappa",&kappa);
  register_var("re_tau",&re_tau);

}// done with constructor

void MASA::rans_sa::init_var()
{

  // currently randomly generated -- these are placeholders for now!
  double twothirds=(double)2/3;

  set_var("cb1",      0.1355);
  set_var("cb2",       0.622);
  set_var("cv1",         7.1);
  set_var("cw2",         0.3);
  set_var("cw3",           2);
  set_var("sigma", twothirds);
  set_var("kappa",      0.41);
  set_var("re_tau",     1000);

} // done with variable initializer

/* ------------------------------------------------
 *
 *         source and analytical terms
 *
 * -----------------------------------------------
 */ 

double MASA::rans_sa::eval_q_u(double eta)
{
  // this is just the momentum equation simplified for channel
  double source_u;
  // achtung:: do we want to add one here?
  source_u = d2u(eta)/re_tau + dnu(eta)*du(eta) + nu(eta)*d2u(eta)+1;
  return source_u;
}

double MASA::rans_sa::eval_q_v(double eta)
{
  double source_v;
  source_v = production(eta) - destruction(eta) + transport(eta);
  return source_v;
}

double MASA::rans_sa::eval_an_u(double eta)
{
  double u_an;
  u_an = u(eta);
  return u_an;
}

double MASA::rans_sa::eval_an_v(double eta)
{
  double v_an;
  v_an = nu(eta);
  return v_an;

}

/* ------------------------------------------------
 *
 *         manufactured solutions and derivatives
 *
 * -----------------------------------------------
 */ 

//velocity
double MASA::rans_sa::u(double eta)
{

  return a1*eta*(1-.5*eta);
  
}

double MASA::rans_sa::du(double eta)
{
  return a1-a1*eta;
}

double MASA::rans_sa::d2u(double eta)
{
  return -a1;
}

// eddy viscosity
double MASA::rans_sa::nu(double eta)
{
  return b1*eta - .5*(etam-1)*b1*eta*eta/etam + b1*pow(eta,3)/(3*etam);
}

double MASA::rans_sa::dnu(double eta)
{
  return b1 - (etam-1)*b1*eta/etam + b1*eta*eta/etam;
}

double MASA::rans_sa::d2nu(double eta)
{
  return -(etam-1)*b1/etam + 2*b1*eta/etam;
}



/* ------------------------------------------------
 *
 *         SA model:
 *         production,destruction,transport
 *
 * -----------------------------------------------
 */ 

// pieces of SA model:
double MASA::rans_sa::production(double eta)
{

  return cb1*s(eta)*nu(eta);
  
}

double MASA::rans_sa::destruction(double eta)
{

  return cw1()*fw(eta)*pow((nu(eta)/eta),2);

}

// some rearranging here, for transport. might look funky to SA specialists.
double MASA::rans_sa::transport(double eta)
{

  return (1/sigma)* ((1/re_tau + nu(eta)) * d2nu(eta) + (1+cb2)*pow(dnu(eta),2));

}


/* ------------------------------------------------
 *
 *         composite functions 
 *
 * composite functions needed to assemble a SA model
 * -----------------------------------------------
 */ 

double MASA::rans_sa::cw1()
{

  return cb1/kappa*kappa + (1+ cb2)/sigma;

}

double MASA::rans_sa::chi(double eta)
{
  
  return nu(eta)*re_tau;

}

double MASA::rans_sa::fv1(double eta)
{

  return pow(chi(eta),3)/(pow(chi(eta),3) + pow(cv1,3));

}

double MASA::rans_sa::fv2(double eta)
{

  return 1- (chi(eta)/(1+chi(eta)*fv1(eta)));

}

// model term for magnitude of mean vorticity 
double MASA::rans_sa::s(double eta)
{
  
  return du(eta) + nu(eta)*fv2(eta)/(kappa*kappa*eta*eta);
  
}

double MASA::rans_sa::r(double eta)
{

  return nu(eta)/(s(eta)*kappa*kappa*eta*eta);

}

double MASA::rans_sa::g(double eta)
{

  return r(eta) + cw3*(pow(r(eta),6) - r(eta));

}

// wall function
double MASA::rans_sa::fw(double eta)
{

  return g(eta)*pow((1+pow(cw3,6))/(pow(g(eta),6)+pow(cw3,6)),1/6.);

}

double MASA::rans_sa::vt(double eta)
{

  return nu(eta)*fv1(eta);

}
