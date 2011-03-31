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
//
// tests.h: helper routines for testing
//
// $Id: masa.h.in 19231 2011-03-30 01:16:36Z nick $
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <cstdio>
#include <config.h>
#include <masa.h>
#include <math.h> 
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <limits>

template<typename Scalar>
Scalar threshcheck(Scalar x, Scalar thresh)
{
  if(isnan(x))
    {
      std::cout << "MASA REGRESSION FAILURE:: nan found!\n";
      exit(1);
    }
  
  if(x > thresh)
    {
      std::cout << "\nMASA REGRESSION TEST FAILED!\n";
      std::cout << "Exceeded Threshold by: " << x << endl;
      exit(1);
    }
  return 0;  
}

template<typename Scalar>
Scalar tester(Scalar a)
{
  
  a = 4.4;
  return a;

}


// heat sc
//template<typename Scalar> Scalar SourceQ_t_1d(Scalar, Scalar, Scalar);
//template<typename Scalar> Scalar Source_t_1d_exact(Scalar A_x,Scalar x)
//template<typename Scalar> Scalar SourceQ_t_2d (Scalar x,Scalar y,Scalar A_x,Scalar B_y,Scalar k_0)
