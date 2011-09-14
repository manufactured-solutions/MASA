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
// $Id$
//--------------------------------------------------------------------------

/*!
 * \file laplacian.c
 *
 * \brief A simple standalone tutorial program which solves a
 * Laplacian using finite-differening and uses the MASA library to aid
 * in verification.  For clarity, this program is self contained with
 * all required source files present in this directory (laplacian.h,
 * laplacian.c, and laplacian_utils.c) Note that for clarity, we use
 * simple matrix data structures and form a full dense system. In
 * addition, no particular performance considerations are applied.
 *
 */

#include "laplacian.h"

int main(int argc, char *argv[])
{
  int n;
  double length;
  pstruct model;		/* primary model data structure */
  int solver;

  /* Parse command-line */

  if(argc < 3)
    {
      printf("\nUsage: laplacian [num_pts] [length] [gauss/cg]\n\n");
      printf("\"num_pts\" is the desired number of mesh points \n");
      printf("\"length\" is the physical length-scale dimension in one direction\n");
      printf("\"gauss/cg\" is which method is used: \n");
      printf("                                    0 --> Gauss-Seidel \n"); 
      printf("                                    1 --> Conjugate Gradient. \n\n");
      exit(1);
    }
  else if(argc == 3)    
    {
      n      = atoi(argv[1]);
      length = (double) atof(argv[2]);
      solver = 0; // default to gauss
    }
  else 
    {
      n      = atoi(argv[1]);
      length = (double) atof(argv[2]);
      solver = atoi(argv[3]);
      
      if(solver > 1 || solver < 0)
	{
	  printf("solver only accepts 0,1 for input");
	  exit(1);
	}      
    }

  /* Problem Initialization */

  problem_initialize (n,length,&model);  
  assemble_matrix    (central_2nd_order,&model);		
  //assemble_matrix    (central_4th_order,&model);		
  init_masa          (&model);
  apply_bcs          (&model);

  print_matrix(&model);

  /* Solve */

  if(solver == 1)
    {
      solve_cg   (&model);  
    }
  else
    {
      solve_gauss(&model);
    }

  /* Compute Error */

  printf("\n** Error Analysis\n");
  printf("   --> npts     = %i\n",model.npts);
  printf("   --> h        = %12.5e\n",model.h);
  printf("   --> l2 error = %12.5e\n",compute_l2_error(&model));

  return 0;
}

