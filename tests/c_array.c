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
  double vec1[n];
  double vec2[n*10];

  for(i=0;i<n;i++)
    {
      vec1[i]=2;      
    }

  // initialize the problem
  err += masa_init("radiation","radiation_integrated_intensity");
  err += masa_init_param();

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa; these tests currently just verify functions
  // run successfully.
  freopen("/dev/null","w",stdout);

  err += masa_display_array();
  masa_set_array("vec_amp",&n,vec1);  
  //err += masa_display_array();

  masa_get_array("vec_mean",&n,vec2);
  
  if(n != 25) 
    {
      printf("masa regression error in c-vector interface!\n");
      printf("size was %i\n",n);
      return 1;
    }
  
  return 0;

}
