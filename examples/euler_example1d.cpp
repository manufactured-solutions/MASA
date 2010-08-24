// $License$
// $Author$
// $Id$ 

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

// this is an example of the MASA API used for calling the 1D euler equation

#include <masa.h>
#include <iostream>
#include <fstream>

using namespace MASA;

int main()
{
  // declarations
  double tempx,tempy;
  double ufield,vfield,efield,rho;	
  double u_an,v_an,p_an,rho_an;

  //problem size
  double lx,dx;
  int nx;

  // initialize
  nx = 10;  // number of points
  lx=1;     // length
  dx=double(lx/nx);

  // initialize the problem
  masa_init("euler-example","euler_1d");

  // initialize the default parameters
  masa_init_param();
  masa_sanity_check();

  // evaluate source terms over the domain (0<x<1)
  for(int i=0;i<nx;i++)
    {  
      tempx=i*dx;

      // evaluate source terms
      masa_eval_u_source  (tempx,tempy,&ufield);
      masa_eval_e_source  (tempx,tempy,&efield);
      masa_eval_rho_source(tempx,tempy,&rho);
      
      //evaluate analytical solution
      masa_eval_u_an        (tempx,tempy,&u_an);
      masa_eval_p_an        (tempx,tempy,&p_an);
      masa_eval_rho_an      (tempx,tempy,&rho_an);
      
    }
  
}// end program
