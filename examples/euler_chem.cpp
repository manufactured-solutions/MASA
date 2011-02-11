// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
// $Author: nick $
// $Id: euler_example.cpp 17232 2011-02-07 23:35:22Z nick $
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

typedef double Scalar;

int main()
{

  int    nx = 10;  // number of points
  int    lx = 1;     // length
  Scalar dx=lx/nx;

  // initialize the problem
  masa_init<Scalar>("euler-chemistry-example","euler_chem_1d");

  // initialize the default parameters
  masa_init_param<Scalar>();

  // intialize the various parameters required for Euler 2D
  // call the sanity check routine 
  // (tests that all variables have been initialized)
  masa_sanity_check<Scalar>();

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(int i=0;i<nx;i++)
    {
      

    } // done with loop

}// done 
