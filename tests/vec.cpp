// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
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
// vec.cpp : program that tests vector capability in masa
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>
#include <masa.h>
#include <cstdio>
#include <iostream>

using namespace MASA;
using namespace std;

typedef double Scalar;

int main()
{
  int err=0;
  std::vector<Scalar> vec1;
  std::vector<Scalar> vec2;

  vec2.resize(7);

  // initialize the problem
  err += masa_init<Scalar>("radiation","radiation_integrated_intensity");
  err += masa_init_param<Scalar>();

  // reroute stdout for regressions: TODO remove when logger mechanism
  // is used inside masa.
  freopen("/dev/null","w",stdout);

  err += masa_display_vec<Scalar>();

  masa_set_vec<Scalar>("vec_mean",vec1);
  masa_set_vec<Scalar>("vec_stdev",vec2);
  
  // deliberately ensure that the sanity check catches the vector
  // has not been initialized!
  err += masa_sanity_check<Scalar>();

  if(err != 1)
    {
      cout << "regression test failed: masa_sanity check not catching uninitialized vectors!\n";
      return 1;
    }

  err += masa_display_vec<Scalar>();

  // need to add test checking that the get_vec with a bad string works too!
  masa_get_vec<Scalar>("vec_amp",vec2);
  
  if(vec2.size() != 25)
    {
      cout << "masa regression failure! masa_get_vec failed!\n";
      cout << "size is: " << vec2.size() << endl;
      return 1;
    }

  return 0;

}
