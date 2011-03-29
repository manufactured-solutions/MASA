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
// $Id: uninit.cpp 18117 2011-02-28 05:10:35Z nick $
//
// c_array.c : program that tests array capability in masa
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{
  int i;
  int err=0;
  int n = 6;
  double n2 = 6;
  double vec1[n];
  double vec2[n];

  for(i=0;i<n;i++)
    {
      vec1[i]=2;
      
      vec2[i]=1;
    }
  // initialize the problem
  err += masa_init("radiation","radiation_integrated_intensity");

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.
  freopen("/dev/null","w",stdout);

  //err += masa_display_array();
  printf("n is: %i\n",n);
  printf("n2 is: %g\n",n2);
  //masa_get_array("vec_mean",n,&vec1);
  err += masa_display_array();

  masa_set_array("vec_mean",&n,&vec1);
  
  err += masa_display_array();

  //masa_get_array("vec_mean",&n,&vec2);

  //printf("n is: %i\n",n);
  //for(i=0;i<n;i++)
  //    {
  //        printf("vec is: %g\n",vec2[i]);
  //}

  return 0;

}
