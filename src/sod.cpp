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
// sod.cpp:  These are the MASA class member functions and constructors
//           For the 1D Sod Solutions
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

//---------------------------------------------------------------------
// Koomie Note: found this online somewhere.  Am leaving it mostly in 
// its original form but converting to C for convenience. Also
// converting to a form which provides output as a function of 
// time and spatial location.  
//
// Output: out_soln[0] = rho
//         out_soln[1] = rho*u
//         out_soln[2] = rho*e
//
// 1/26/09
//---------------------------------------------------------------------

//  Program calculates the exact solution to Sod-shock tube class
//  problems -- namely shock tubes which produce shocks, contact
//  discontinuities, and rarefraction waves.
//
//  Solution is computed at locations x at time t. (Though   
//  due to self-similarity, the exact solution is identical for
//  identical values of x/t).
//
//  NOTE : Since the post-shock flow is nonadiabatic, whereas
//  the flow inside the rarefraction fan is adiabatic, the problem
//  is not left-right symmetric. In particular, the high-density
//  initial state MUST BE input on the left side.
// 

#include <masa_internal.h> 
#include <iostream>
#include <limits>

using namespace MASA;

template <typename Scalar>
MASA::sod_1d<Scalar>::sod_1d()
{
  this->mmsname = "sod_1d";
  this->dimension=1;

  this->register_var("Gamma",&Gamma);
  this->register_var("mu",&mu);

}//done with constructor

template <typename Scalar>
int MASA::sod_1d<Scalar>::init_var()
{
  int err = 0;

  err += this->set_var("Gamma",1.4e0);
  Scalar temp = ((Gamma - 1.e0) / (Gamma + 1.e0));
  err += this->set_var("mu",temp);

  return err;

} // done with variable initializer

template <typename Scalar>
Scalar MASA::sod_1d<Scalar>::eval_q_t(Scalar x)
{
  // Define the Sod problem initial conditions for the left and right states.

  pl = 1.e0;
  pr = 0.125;

  rhol = 1.e0;
  rhor = 0.125;

  // Define sound speeds for the left and right sides of tube.

  cl = sqrt (Gamma * pl / rhol);
  cr = sqrt (Gamma * pr / rhor);

  // will hit masa_exit --root not bracketed!
  return rtbis(1.,20.,1.,100);

}

template <typename Scalar>
Scalar MASA::sod_1d<Scalar>::eval_q_t(Scalar x,Scalar t)
{
  // Define the Sod problem initial conditions for the left and right states.

  pl = 1.e0;
  pr = 0.125;

  rhol = 1.e0;
  rhor = 0.125;

  // Define sound speeds for the left and right sides of tube.

  cl = sqrt (Gamma * pl / rhol);
  cr = sqrt (Gamma * pr / rhor);

  // will hit 'too many bisections' 
  return rtbis(-1,2,1,1);

}

template <typename Scalar>
Scalar MASA::sod_1d<Scalar>::eval_q_rho(Scalar x,Scalar t)
{
  // xmax determines the size of the computational domain (-xmax, +xmax).
  // numcells determines the number of cells in the output table.          
  // Scalar xmax 	= 5.e0;

  Scalar pm;
  Scalar rhoml, vs, vt, rhomr, vm, density;

  // Define the Sod problem initial conditions for the left and right states.

  pl = 1.e0;
  pr = 0.125;

  rhol = 1.e0;
  rhor = 0.125;

  // Define sound speeds for the left and right sides of tube.

  cl = sqrt (Gamma * pl / rhol);
  cr = sqrt (Gamma * pr / rhor);

  // Solve for the postshock pressure pm.

  pm = rtbis (pr, pl, std::numeric_limits<Scalar>::epsilon(),100);

  // Define the density to the left of the contact discontinuity rhoml.
 
  rhoml = pow(rhol * (pm / pl),(1.e0 / Gamma));

  // Define the postshock fluid velocity vm.

  vm = 2.e0 * cl / (Gamma - 1.e0) * (1.e0 - pow((pm / pl),( (Gamma - 1.e0) / (2.e0 * Gamma) )));

  // Define the postshock density rhomr.

  rhomr = rhor *  ( (pm + mu * pr) / (pr + mu * pm) );

  // Define the shock velocity vs.

  vs = vm / (1.e0 - rhor / rhomr);

  // Define the velocity of the rarefraction tail, vt.

  vt = cl - vm / (1.e0 - mu);

  // Output tables of density, velocity, and pressure at time t.

  //  for(i=0;i<numnodes;i++)
  //    {
  //x = - xmax +  2.e0 * xmax * i / numnodes;
 
      if (x <= (-cl*t) )
	  density = rhol;
      else if (x <= (-vt*t) )
	  density = pow(rhol * (-mu * (x / (cl * t) ) + (1 - mu) ),(2.e0 / (Gamma - 1.e0)));
      else if (x <= (vm*t) )
	  density = rhoml;
      else if (x <= (vs*t) )
	density = rhomr;
      else
	density = rhor;
      
      Scalar out_soln = density;
      return out_soln;

}

// return product of rho and u
// 
template <typename Scalar>
Scalar MASA::sod_1d<Scalar>::eval_q_rho_u(Scalar x,Scalar t)
{
  // xmax determines the size of the computational domain (-xmax, +xmax).
  // numcells determines the number of cells in the output table.          
  // Scalar xmax 	= 5.e0;

  Scalar pm;
  Scalar rhoml, vs, vt, rhomr, vm, density, velocity;

  // Define the Sod problem initial conditions for the left and right states.

  pl = 1.e0;
  pr = 0.125;

  rhol = 1.e0;
  rhor = 0.125;

  // Define sound speeds for the left and right sides of tube.

  cl = sqrt (Gamma * pl / rhol);
  cr = sqrt (Gamma * pr / rhor);

  // Solve for the postshock pressure pm.

  pm = rtbis (pr, pl, std::numeric_limits<Scalar>::epsilon(),100);

  // Define the density to the left of the contact discontinuity rhoml.
 
  rhoml = pow(rhol * (pm / pl),(1.e0 / Gamma));

  // Define the postshock fluid velocity vm.

  vm = 2.e0 * cl / (Gamma - 1.e0) * (1.e0 - pow((pm / pl),( (Gamma - 1.e0) / (2.e0 * Gamma) )));

  // Define the postshock density rhomr.

  rhomr = rhor *  ( (pm + mu * pr) / (pr + mu * pm) );

  // Define the shock velocity vs.

  vs = vm / (1.e0 - rhor / rhomr);

  // Define the velocity of the rarefraction tail, vt.

  vt = cl - vm / (1.e0 - mu);

  // Output tables of density, velocity, and pressure at time t.

  //  for(i=0;i<numnodes;i++)
  //    {
  //x = - xmax +  2.e0 * xmax * i / numnodes;
 
      if (x <= (-cl*t) )
	  density = rhol;
      else if (x <= (-vt*t) )
	  density = pow(rhol * (-mu * (x / (cl * t) ) + (1 - mu) ),(2.e0 / (Gamma - 1.e0))) ;
      else if (x <= (vm*t) )
	  density = rhoml;
      else if (x <= (vs*t) )
	density = rhomr;
      else
	density = rhor;
      
      if (x <= (-cl*t) )
	velocity = 0.0;
      else if (x <= (-vt*t) )
	velocity = (1 - mu) * (x / t + cl);
      else if (x <= (vs*t) )
	velocity = vm;
      else 
	velocity = 0.0;

      Scalar out_soln = density*velocity;
      return out_soln;
}

// return pressure
// 
template <typename Scalar>
Scalar MASA::sod_1d<Scalar>::eval_q_p(Scalar x,Scalar t)
{

  Scalar pm, pressure;
  Scalar vs, vt, rhomr, vm;

  // Define the Sod problem initial conditions for the left and right states.

  pl = 1.e0;
  pr = 0.125;

  rhol = 1.e0;
  rhor = 0.125;

  // Define sound speeds for the left and right sides of tube.

  cl = sqrt (Gamma * pl / rhol);
  cr = sqrt (Gamma * pr / rhor);

  // Solve for the postshock pressure pm.

  pm = rtbis (pr, pl, std::numeric_limits<Scalar>::epsilon(),100);

  // Define the postshock fluid velocity vm.

  vm = 2.e0 * cl / (Gamma - 1.e0) * (1.e0 - pow((pm / pl),( (Gamma - 1.e0) / (2.e0 * Gamma) )));

  // Define the postshock density rhomr.

  rhomr = rhor *  ( (pm + mu * pr) / (pr + mu * pm) );

  // Define the shock velocity vs.

  vs = vm / (1.e0 - rhor / rhomr);

  // Define the velocity of the rarefraction tail, vt.

  vt = cl - vm / (1.e0 - mu);

  // Output tables of pressure at time t.
  //  for(i=0;i<numnodes;i++)
  //    {
  //x = - xmax +  2.e0 * xmax * i / numnodes;
  
  if (x <= (-cl*t) )
    pressure = pl;
  else if (x <= (-vt*t) )
    pressure = pow(pl * (-mu * (x / (cl * t) ) + (1 - mu) ),(2.e0 * Gamma / (Gamma - 1.e0)));
  else if (x <= vs * t)
    pressure = pm;
  else            
    pressure = pr;
  
  Scalar out_soln = pressure;
  return out_soln;
}

template <typename Scalar>
Scalar MASA::sod_1d<Scalar>::func(Scalar pm)
{

///////////////////////////////////////////////////////////////////////
//
//  func is obtained from an identity matching the post-shocked      
//  pressure to the post-rarefraction pressure (true since there is
//  no pressure jump across the contact discontinuity). We use it to
//  numerically solve for pm given the left and right initial states.
//
///////////////////////////////////////////////////////////////////////
 
  Scalar myval;
 
  myval = -2*cl*(1 - pow((pm/pl),((-1 + Gamma)/(2*Gamma))))/
    (cr*(-1 + Gamma)) +
    (-1 + pm/pr)*sqrt((1 - mu)/(Gamma*(mu + pm/pr)));
  
  return(myval);
}

//////////////////////////////////////////////////////////////////////////
//
// rtbis is borrowed from Numerical Recipes. It is a bisection algorithm,
// which we use to solve for pm using a call to func.
//
// Note that the arguments to rtbis have been altered and the value of
// JMAX is now passed to the function. 
// 
// Otherwise, it is identical to the NR version.
//
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>

template <typename Scalar>
Scalar MASA::sod_1d<Scalar>::rtbis(Scalar x1,Scalar x2,Scalar xacc,int JMAX)
{
  int j;
  Scalar dx,f,fmid,xmid;
  Scalar myval;

  fmid=func(x2);
  f=func(x1);
  
  if(f*fmid >= 0.)
    {
      printf("MASA ERROR:: root must be bracketed in rtbis (sod)\n");
      masa_exit(1);
    }
     
  if(f < 0.)
    {
      myval=x1;
      dx=x2-x1;
    }
  else
    {
      myval=x2;
      dx=x1-x2;
    }
  
  for(j=0;j<JMAX;j++)
    {
      dx=dx*5.e-1;
      xmid=myval+dx;
      fmid=func(xmid);
      if(fmid <= 0.) myval=xmid;
      if(fabs(dx) < xacc || fmid == 0.) 
	return(myval);
    }
  
  printf("MASA Error:: Too many bisection in rtbis (sod)\n");
  return(-1);

}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::sod_1d);
