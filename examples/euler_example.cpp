// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// $Author$
// $Id$
//
// euler_example.cpp:
// this is an example of the API used for calling the 2D euler equation
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>
#include <iostream>
#include <fstream>

using namespace MASA;

int main()
{
  // declarations
  double x,y;
  double tempx,tempy;

  double ufield;
  double vfield;
  double efield;
  double rho;

  double u_an;
  double v_an;
  double p_an;
  double rho_an;

  //problem size
  double lx,ly;
  double dx,dy;
  int nx,ny;

  // initialize
  nx = 10;  // number of points
  ny = 10;  
  lx=1;     // length
  ly=1; 

  dx=double(lx/nx);
  dy=double(ly/ny);

  // initialize the problem
  masa_init("euler-example","euler_2d");

  // initialize the default parameters
  masa_init_param();

  // intialize the various parameters required for Euler 2D
  // call the sanity check routine 
  // (tests that all variables have been initialized)
  masa_sanity_check();

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(int i=0;i<nx;i++)
    for(int j=0;j<nx;j++)
      {  
	tempx=i*dx;
	tempy=j*dy;

	// evaluate source terms
	masa_eval_u_source  (tempx,tempy,&ufield);
	masa_eval_v_source  (tempx,tempy,&vfield);
	masa_eval_e_source  (tempx,tempy,&efield);
	masa_eval_rho_source(tempx,tempy,&rho);
	
	//evaluate analytical solution
	masa_eval_u_an        (tempx,tempy,&u_an);
	masa_eval_v_an        (tempx,tempy,&v_an);
	masa_eval_p_an        (tempx,tempy,&p_an);
	masa_eval_rho_an      (tempx,tempy,&rho_an);

      }

}// end program
