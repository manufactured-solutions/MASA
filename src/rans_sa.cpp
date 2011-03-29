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

template <typename Scalar>
MASA::rans_sa<Scalar>::rans_sa()
{
  this->mmsname = "rans_sa";
  this->dimension=1;
  
  // value of 'peak' in eddy viscosity
  etam=Scalar(6)/Scalar(10);

  // u parameter
  a1=2;

  // eddy viscosity parameter
  b1=1;

  this->register_var("cb1",  &cb1);
  this->register_var("cb2",  &cb2);
  this->register_var("cv1",  &cv1);
  this->register_var("cw2",  &cw2);
  this->register_var("cw3",  &cw3);
  this->register_var("sigma",&sigma);
  this->register_var("kappa",&kappa);
  this->register_var("re_tau",&re_tau);

  this->register_var("cv2",&cv2);
  this->register_var("cv3",&cv3);

  this->init_var();

}// done with constructor

template <typename Scalar>
int MASA::rans_sa<Scalar>::init_var()
{
  int err = 0;

  // Default values - set with double precision, but they're independent
  // variables so that's okay
  Scalar twothirds=Scalar(2)/Scalar(3);

  err += this->set_var("cb1",      0.1355);
  err += this->set_var("cb2",       0.622);
  err += this->set_var("cv1",         7.1);
  err += this->set_var("cw2",         0.3);
  err += this->set_var("cw3",           2);
  err += this->set_var("sigma", twothirds);
  err += this->set_var("kappa",      0.41);
  //this->set_var("re_tau",     1000);

  // use lower re_tau for default.  this
  // lengthens in the inner region and
  // allows us to achieve asymptotic
  // results with coarser meshes.
  err += this->set_var("re_tau",     100); 

  err += this->set_var("cv2",         0.7);
  err += this->set_var("cv3",         0.9);

  return err;

} // done with variable initializer

/* ------------------------------------------------
 *
 *         source and analytical terms
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::eval_q_u(Scalar eta)
{
  // this is just the momentum equation simplified for channel
  Scalar source_u;
  // achtung:: do we want to add one here?
  source_u = d2u()/re_tau + dvt(eta)*du(eta) + vt(eta)*d2u()+1;
  return source_u;
}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::eval_q_v(Scalar eta)
{
  Scalar source_v;
  source_v = production(eta) - destruction(eta) + transport(eta);
  return source_v;
}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::eval_exact_u(Scalar eta)
{
  Scalar exact_u;
  exact_u = u(eta);
  return exact_u;
}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::eval_exact_v(Scalar eta)
{
  Scalar exact_v;
  exact_v = nu(eta);
  return exact_v;

}

/* ------------------------------------------------
 *
 *         manufactured solutions and derivatives
 *
 * -----------------------------------------------
 */ 

//velocity
template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::u(Scalar eta)
{

  return a1*eta*(1-.5*eta);
  
}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::du(Scalar eta)
{
  return a1-a1*eta;
}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::d2u()
{
  return -a1;
}

// eddy viscosity
template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::nu(Scalar eta)
{
  return b1*eta - .5*(etam+1)*b1*eta*eta/etam + b1*pow(eta,3)/(3*etam);
}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::dnu(Scalar eta)
{
  return b1 - (etam+1)*b1*eta/etam + b1*eta*eta/etam;
}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::d2nu(Scalar eta)
{
  return -(etam+1)*b1/etam + 2*b1*eta/etam;
}



/* ------------------------------------------------
 *
 *         SA model:
 *         production,destruction,transport
 *
 * -----------------------------------------------
 */ 

// pieces of SA model:
template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::production(Scalar eta)
{

  return cb1*s(eta)*nu(eta);
  
}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::destruction(Scalar eta)
{

  return cw1()*fw(eta)*pow((nu(eta)/eta),2);

}

// some rearranging here, for transport. might look funky to SA specialists.
template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::transport(Scalar eta)
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

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::cw1()
{

  return cb1/(kappa*kappa) + (1+ cb2)/sigma;

}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::chi(Scalar eta)
{
  
  return nu(eta)*re_tau;

}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::fv1(Scalar eta)
{

  return pow(chi(eta),3)/(pow(chi(eta),3) + pow(cv1,3));

}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::fv2(Scalar eta)
{

  return 1- (chi(eta)/(1+chi(eta)*fv1(eta)));

}

// model term for magnitude of mean vorticity 
template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::s(Scalar eta)
{

  // original SA
  //return du(eta) + nu(eta)*fv2(eta)/(kappa*kappa*eta*eta);


  // modified SA.  this is modified formulation of S_{sa} due to
  // Johnson and Allmaras.  See ChIPS model document for more details.

  Scalar Sbar = nu(eta)*fv2(eta)/(kappa*kappa*eta*eta);

  if( Sbar >= - cv2*du(eta) )
    return du(eta) + Sbar;
  else
    return du(eta) + du(eta)*(cv2*cv2*du(eta) + cv3*Sbar)/((cv3-2*cv2)*du(eta) - Sbar);
    
}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::r(Scalar eta)
{
  // Note: Limiting r to maximum of 10.0 as recommended in original SA
  // paper (and implemented in ChIPS).

  Scalar rtmp = nu(eta)/(s(eta)*kappa*kappa*eta*eta);

  if (rtmp>10.0)
    return 10.0;
  else
    return rtmp;
}

template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::g(Scalar eta)
{

  return r(eta) + cw2*(pow(r(eta),6) - r(eta));

}

// wall function
template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::fw(Scalar eta)
{

  return g(eta)*pow((1.+pow(cw3,6))/(pow(g(eta),6)+pow(cw3,6)),1/Scalar(6));

}

// first derivative of vt
template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::dvt(Scalar eta)
{

  // oliver: I can't reproduce this expresssion... tried again below

//   // see modeling document for this beast
//   Scalar a = 3* pow(cv1,3) * pow(re_tau,3) * pow(nu(eta),4);
//   Scalar b = pow(re_tau,3)*pow(nu(eta),4) + nu(eta)*pow(cv1,3);
  
//   return a/(b*b);

  Scalar n = nu(eta);
  Scalar a = cv1/re_tau;

  Scalar num = n*n*n * ( n*n*n + 4.0*a*a*a ) * dnu(eta);
  Scalar rde = n*n*n + a*a*a;

  return num/(rde*rde);

}

//this is turbulent eddy viscosity
template <typename Scalar>
Scalar MASA::rans_sa<Scalar>::vt(Scalar eta)
{

  return nu(eta)*fv1(eta);

}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::rans_sa);
