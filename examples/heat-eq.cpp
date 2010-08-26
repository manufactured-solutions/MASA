// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// $Author$
// $Id$
//
// heat-eq.cpp: example of the API used for the 2D heat equation
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
  double tfield;
  double tempx,tempy;

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
  masa_init("heat equation example","heateq_2d_steady_const");

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

	masa_eval_t_source  (tempx,tempy,&tfield);
	
      }

}// end program
