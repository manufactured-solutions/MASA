 /*--------------------------------------------------------------------------
  *--------------------------------------------------------------------------
  *
  * Copyright (C) 2010 The PECOS Development Team
  *
  * Please see http://pecos.ices.utexas.edu for more information.
  *
  * This file is part of MASA.
  *
  * MASA is free software: you can redistribute it and/or modify it under
  * the terms of the GNU Lesser General Public License as published by the Free
  * Software Foundation, either version 3 of the License, or (at your option)
  * any later version.
  *
  * MASA is distributed in the hope that it will be useful, but WITHOUT ANY
  * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  * details.
  *
  * You should have received a copy of the GNU Lesser General Public License along 
  * with MASA.  If not, see <http://www.gnu.org/licenses/>.
  *
  *--------------------------------------------------------------------------
  
  MASA -- Manufactured Analytical Solutions Abstraction Library

  A software interface that provides access to all manufactured solutions to 
  be used by various models throughout the center.
  
  *--------------------------------------------------------------------------
  */  

//
//   These are the MASA class member functions and constructors
//   For the 1D Sod Solutions

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

using namespace MASA;

double pl, pr, rhol, rhor, cl, cr;

MASA::sod_1d::sod_1d()
{
  mmsname = "sod_1d";
  dimension=1;

  register_var("Gamma",&Gamma);
  register_var("mu",&mu);

}//done with constructor

void MASA::sod_1d::init_var()
{
  set_var("Gamma",1.4e0);
  double temp = ((Gamma - 1.e0) / (Gamma + 1.e0));
  set_var("mu",temp);

} // done with variable initializer

double MASA::sod_1d::eval_q_rho(double x,double t)
{
  double xm;

  // xmax determines the size of the computational domain (-xmax, +xmax).
  // numcells determines the number of cells in the output table.          
  // double xmax 	= 5.e0;

  double pm, pressure;
  double rhoml, vs, vt, rhomr, vm, density, velocity;
  int i;
  FILE *fp;

  // Define the Sod problem initial conditions for the left and right states.

  pl = 1.e0;
  pr = 0.125;

  rhol = 1.e0;
  rhor = 0.125;

  // Define sound speeds for the left and right sides of tube.

  cl = sqrt (Gamma * pl / rhol);
  cr = sqrt (Gamma * pr / rhor);

  // Solve for the postshock pressure pm.

  pm = rtbis (pr, pl, 1.e-16);

  // Define the density to the left of the contact discontinuity rhoml.
 
  rhoml = pow(rhol * (pm / pl),(1.e0 / Gamma));

  // Define the postshock fluid velocity vm.

  vm = 2.e0 * cl / (Gamma - 1.e0) * (1.e0 - pow((pm / pl),
               ( (Gamma - 1.e0) / (2.e0 * Gamma) )));

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
	  density = pow(rhol * (-mu * (x / (cl * t) ) + (1 - mu) ),
			  (2.e0 / (Gamma - 1.e0))) ;
      else if (x <= (vm*t) )
	  density = rhoml;
      else if (x <= (vs*t) )
	density = rhomr;
      else
	density = rhor;
      
      if (x <= (-cl*t) )
	pressure = pl;
      else if (x <= (-vt*t) )
	pressure = pow(pl * (-mu * (x / (cl * t) ) + (1 - mu) ),
			 (2.e0 * Gamma / (Gamma - 1.e0)));
      else if (x <= vs * t)
	pressure = pm;
      else            
	pressure = pr;
      
      if (x <= (-cl*t) )
	velocity = 0.0;
      else if (x <= (-vt*t) )
	velocity = (1 - mu) * (x / t + cl);
      else if (x <= (vs*t) )
	velocity = vm;
      else 
	velocity = 0.0;

      double out_soln = density;
      return out_soln;

}

// return product of rho and u
// 
double MASA::sod_1d::eval_q_rho_u(double x,double t)
{
  double xm;

  // xmax determines the size of the computational domain (-xmax, +xmax).
  // numcells determines the number of cells in the output table.          
  // double xmax 	= 5.e0;

  double pm, pressure;
  double rhoml, vs, vt, rhomr, vm, density, velocity;
  int i;
  FILE *fp;

  // Define the Sod problem initial conditions for the left and right states.

  pl = 1.e0;
  pr = 0.125;

  rhol = 1.e0;
  rhor = 0.125;

  // Define sound speeds for the left and right sides of tube.

  cl = sqrt (Gamma * pl / rhol);
  cr = sqrt (Gamma * pr / rhor);

  // Solve for the postshock pressure pm.

  pm = rtbis (pr, pl, 1.e-16);

  // Define the density to the left of the contact discontinuity rhoml.
 
  rhoml = pow(rhol * (pm / pl),(1.e0 / Gamma));

  // Define the postshock fluid velocity vm.

  vm = 2.e0 * cl / (Gamma - 1.e0) * (1.e0 - pow((pm / pl),
               ( (Gamma - 1.e0) / (2.e0 * Gamma) )));

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
	  density = pow(rhol * (-mu * (x / (cl * t) ) + (1 - mu) ),
			  (2.e0 / (Gamma - 1.e0))) ;
      else if (x <= (vm*t) )
	  density = rhoml;
      else if (x <= (vs*t) )
	density = rhomr;
      else
	density = rhor;
      
      if (x <= (-cl*t) )
	pressure = pl;
      else if (x <= (-vt*t) )
	pressure = pow(pl * (-mu * (x / (cl * t) ) + (1 - mu) ),
			 (2.e0 * Gamma / (Gamma - 1.e0)));
      else if (x <= vs * t)
	pressure = pm;
      else            
	pressure = pr;
      
      if (x <= (-cl*t) )
	velocity = 0.0;
      else if (x <= (-vt*t) )
	velocity = (1 - mu) * (x / t + cl);
      else if (x <= (vs*t) )
	velocity = vm;
      else 
	velocity = 0.0;

      double out_soln = density*velocity;
      return out_soln;
}

double MASA::sod_1d::func(double pm)
{

///////////////////////////////////////////////////////////////////////
//
//  func is obtained from an identity matching the post-shocked      
//  pressure to the post-rarefraction pressure (true since there is
//  no pressure jump across the contact discontinuity). We use it to
//  numerically solve for pm given the left and right initial states.
//
///////////////////////////////////////////////////////////////////////
 
  double myval;
 
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
// JMAX increased. Otherwise, it is identical to the NR version.
//
//////////////////////////////////////////////////////////////////////////

double MASA::sod_1d::rtbis(double x1,double x2,double xacc)
{
  int JMAX=100;
  int j;
  double dx,f,fmid,xmid;
  double myval;

  fmid=func(x2);

  f=func(x1);
  if(f*fmid >= 0.)
    {
      printf("root must be bracketed in rtbis\n");
      exit(1);
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
  
  printf("** Error: Too many bisection in rtbis\n");
  return(-1);
  

}
