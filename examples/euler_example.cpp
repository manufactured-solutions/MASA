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

// this is an example of the MASA API used for calling the 2D euler equation

#include <masa.h>
#include <iostream>
#include <fstream>

using namespace MASA;

int main()
{
  // declarations
  double param =.8;
  double x,y;
  double ufield;
  double vfield;
  double efield;
  double rho;
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
  masa_init("euler-example","euler_2d");

  // intialize the various parameters required for Euler 2D
  masa_set_param("u_0",param);
  masa_set_param("u_x",param);
  masa_set_param("u_y",param);

  masa_set_param("v_0",param);
  masa_set_param("v_x",param);
  masa_set_param("v_y",param);

  masa_set_param("rho_0",param);
  masa_set_param("rho_x",param);
  masa_set_param("rho_y",param);

  masa_set_param("p_0",param);
  masa_set_param("p_x",param);
  masa_set_param("p_y",param);

  masa_set_param("a_px",param);
  masa_set_param("a_py",param);

  masa_set_param("a_rhox",param);
  masa_set_param("a_rhoy",param);

  masa_set_param("a_ux",param);
  masa_set_param("a_uy",param);

  masa_set_param("a_vx",param);
  masa_set_param("a_vy",param);

  masa_set_param("Gamma",param);
  masa_set_param("mu",param);
  masa_set_param("L",param);


  // call the sanity check routine 
  // (tests that all variables have been initialized)
  masa_sanity_check();

  // evaluate source terms over the domain (0<x<1, 0<y<1) 
  for(int i=0;i<nx;i++)
    for(int j=0;j<nx;j++)
      {  
	tempx=i*dx;
	tempy=j*dy;
	masa_eval_u_source  (tempx,tempy,&ufield);
	masa_eval_v_source  (tempx,tempy,&vfield);
	masa_eval_e_source  (tempx,tempy,&efield);
	masa_eval_rho_source(tempx,tempy,&rho);
      }


}// end program
