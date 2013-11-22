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
// $Id$
//--------------------------------------------------------------------------

#include "laplacian.h"

/*! 
 * \file laplacian_utils.c
 * \brief Utility functions used with laplacian.c
 */					

/*!
 * \fn init_masa(pstruct *model)
 * \brief Initializes MASA library and applies right-hand side forcing function 
 * at all solution points.
 * \param model Pointer to the primary model data-structure
 */

void init_masa(pstruct *model)
{
  int i,j;
  int index;
  double xval, yval;

  printf("\n** Initializing MASA\n");

  masa_init("A C Laplacian Example","laplace_2d");
  
  printf("   --> applying RHS forcing function\n");

  /* Define RHS for solution domain */

  for(i=0;i<model->npts;i++)
    {
      for(j=0;j<model->npts;j++)
	{
	  xval  = (i)*model->h;
	  yval  = (j)*model->h;
	  index = j+(i*model->npts);

	  model->rhs[index] = masa_eval_2d_source_f(xval,yval);
	}
    }
  printf("   --> complete\n");
  return;
}

/*!
 * \fn enforce_dirichlet_bc(const int row_id, const int col_id, const double value, pstruct *model)
 *
 * \brief Enforces a Dirichlet BC in the model system by nullifying
 * all but the diagonal entry for the row corresponding to a given
 * solution value. The model RHS forcing function is set to the
 * desired Dirchlet value and the diagonal is set to unity to enforce
 * the desired constraint.
 *
 * \param row_id Row index to enforce the constraint
 * \param col_id Column index to enforce the constraint
 * \param value  Desired solution value
 * \param model  Pointer to the primary model data-structure
 */

void enforce_dirichlet_bc(const int row_id, const int col_id, const double value, pstruct *model)
{
  int i,j;
  int index;

  index = col_id+(row_id*model->npts);

  for(j=0;j<model->n;j++)
    model->A[index][j] = 0.0;

  model->A[index][index] = 1.0;
  model->rhs[index]      = value;

  return;
}

/*!
 * \fn print_matrix(pstruct *model)
 *
 * \brief Prints model system matrix and right-hand side vector to stdout.
 *
 * \param model  Pointer to the primary model data-structure
 */

void print_matrix(pstruct *model)
{
  int i,j;

  printf("     j =");
  for(j=0;j<model->n;j++)
    printf("%6i ",j);
  printf("\n\n");
  
  for(i=0;i<model->n;i++)
    {
      printf("i = %2i: ",i);
      for(j=0;j<model->n;j++)
	printf("%6.2f ",model->A[i][j]);
      printf("\n");
    }

  printf("\nRHS:\n");
  for(j=0;j<model->n;j++)
    printf("%6.2f\n",model->rhs[j]);

}

/*!
 * \fn assemble_matrix(const int fd_method, pstruct *model)
 *
 * \brief Assembles system matrix entries using a finite-difference approximation. 
 * The fd_method parameter controls the choice of underlying FD stencil.
 *
 * \param fd_method Desired finite-difference stencil
 * \param model     Pointer to the primary model data-structure
 */

void assemble_matrix(const int fd_method, pstruct *model)
{
  int i,j,index,index2;
  const double h_squared = model->h*model->h;

  printf("\n** Assembling linear system\n");

  switch(fd_method)
    {
    case(central_2nd_order):
      {
	model->pad = 1;
	printf("   --> using 2nd-order central difference approximation (stencil width = %i)\n",model->pad);
	
	for(i=0;i<model->npts;i++)
	  for(j=0;j<model->npts;j++)
	    {
	      index = j+(i*model->npts);
	      
	      model->A[index][index] = -4.0/h_squared;
	      
	      if(j > 0)		    /* (i,j-1) */
		model->A[index][index-1] = 1.0/h_squared;
	      
	      if(j < model->n-1 )   /* (i,j+1) */
		model->A[index][index+1] = 1.0/h_squared;
	      
	      if(i > 0)		    /* (i-1,j) */
		{
		  index2 = j+(i-1)*model->npts;
		  model->A[index][index2] = 1.0/h_squared;
		}
	      
	      if(i < model->npts-1) /* (i+1,j) */
		{
		  index2 = j+(i+1)*model->npts;
		  model->A[index][index2] = 1.0/h_squared;
		}
	    }
	break;
      }	/* end central_2nd_order */

    case(central_4th_order):
      {
	model->pad = 2;
	printf("   --> using 4th-order central difference approximation (stencil width = %i)\n",model->pad);
	
	for(i=0;i<model->npts;i++)
	  for(j=0;j<model->npts;j++)
	    {
	      index  = j+(i*model->npts);
	      
	      model->A[index][index] = -60.0/(12.0*h_squared);
	      
	      if(j > 0)	             /* (i,j-1) */
		model->A[index][index-1] = 16.0/(12.0*h_squared);
	      
	      if(j > 1)	             /* (i,j-2) */
		model->A[index][index-2] = -1.0/(12.0*h_squared);
	      
	      if(j < model->n-1)     /* (i,j+1) */
		model->A[index][index+1] = 16.0/(12.0*h_squared);
	      
	      if(j < model->n-2)     /* (i,j+2) */
		model->A[index][index+2] = -1.0/(12.0*h_squared);
	      
	      if(i > 0)		     /* (i-1,j) */
		{
		  index2 = j+(i-1)*model->npts;
		  model->A[index][index2] = 16.0/(12.0*h_squared);
		}
	      
	      if(i > 1)		     /* (i-2,j) */
		{
		  index2 = j+(i-2)*model->npts;
		  model->A[index][index2] = -1.0/(12.0*h_squared);
		}
	      
	      if(i < model->n-1 )    /* (i+1,j) */
		{
		  index2 = j+(i+1)*model->npts;
		  model->A[index][index2] = 16.0/(12.0*h_squared);
		}
	      
	      if(i < model->n-2 )    /* (i+2,j) */
		{
		  index2 = j+(i+2)*model->npts;
		  model->A[index][index2] = -1.0/(12.0*h_squared);
		}
	    }
	break;
      }	/* end central_4th_order */
      

    default:
      printf("\n** Error: unknown finite-difference method requested\n\n");
      exit(1);
	    
    }  /* end switch statement over FD methods */

  printf("   --> assembly complete\n");
  return;
}

void apply_bcs(pstruct *model)
{
  int i,j;
  int index;
  double soln;
  const int n    = model->n;
  const int npts = model->npts;

  double xval,yval;

  assert(model->pad >= 1);

    {
      /* BCs for north boundaries */

      for(i=0;i<model->pad;i++)
	for(j=0;j<model->npts;j++)      
	  {
	    xval = (i)*model->h;
	    yval = (j)*model->h;
	    soln = masa_eval_2d_exact_phi(xval,yval);

	    enforce_dirichlet_bc(i,j,soln,model);
	  }

      /* BCs for south boundaries */

      for(i=model->npts-model->pad;i<=model->npts-1;i++)
	for(j=0;j<model->npts;j++)      
	  {
	    xval = (i)*model->h;
	    yval = (j)*model->h;
	    soln = masa_eval_2d_exact_phi(xval,yval);
	    
	    enforce_dirichlet_bc(i,j,soln,model);
	  }

      /* BCs for west boundaries */

      for(i=0;i<model->npts;i++)      
	for(j=0;j<model->pad;j++)
	  {
	    xval = (i)*model->h;
	    yval = (j)*model->h;
	    soln = masa_eval_2d_exact_phi(xval,yval);
	    
	    enforce_dirichlet_bc(i,j,soln,model);
	  }

      /* BCs for east boundaries */

      for(i=0;i<model->npts;i++)      
	for(j=model->npts-model->pad;j<=model->npts-1;j++)
	  {
	    xval = (i)*model->h;
	    yval = (j)*model->h;
	    soln = masa_eval_2d_exact_phi(xval,yval);
	    
	    enforce_dirichlet_bc(i,j,soln,model);
	  }
    }

  return;
}

void problem_initialize(const int npts, const double length, pstruct *model)
{

/*! \fn           problem_initialize
 *  \brief        Allocates memory for desired problem size.
 *  \param npts    Desired # points in one direction (resulting solution vector is npts*npts)
 *  \param length Dimension of unit square
 *  \param model  Laplacian data structure.
 */

  int i;

  assert(npts   > 1 );
  assert(length > 0.);

  printf("\n** Initializing problem\n");

  model->npts   = npts;
  model->h      = length/(npts-1.0);
  model->n      = npts*npts;

  printf("   --> mesh size           = %-12.3f\n",model->h);

  /* Perform memory allocation */

  model->rhs = (double *) calloc(sizeof(double  ),     model->n);
  model->phi = (double *) calloc(sizeof(double  ),     model->n);
  model->A   = (double **)calloc(sizeof(double *),     model->n);

  if(model->rhs == NULL || model->phi == NULL || model->A  == NULL )
    {
      printf("ERROR: Unable to allocate memory\n");
      exit(0);
    }

  for(i=0;i<model->n;i++)
    {
      model->A[i] = calloc(sizeof(double),model->n);
      if(model->A[i] == NULL)
	{
	  printf("ERROR: Unable to allocate memory for matrix A\n");
	  exit(0);
	}
    }    
  
  printf("   --> memory initialization complete\n");
  return;
}


void solve_gauss(pstruct *model)
{
  int it=0;
  int i,j, itmp;
  double *phi_old;
  int   **sparse_index;
  int    *sparse_count;
  int col_index;
  double diff            = 0.0;
			
  const double  ZERO_TOL = 1.0e-20;
  const double CONVG_TOL = 1.0e-9;
  const int    ITER_MAX  = 30000;
  const int           n  = model->n;

  printf("\n** Solving system using Gauss-Seidel\n");

  phi_old      =        calloc(sizeof(double),model->n);
  sparse_index = (int **)calloc(sizeof(int *),model->n);
  sparse_count = ( int *)calloc(sizeof(int)  ,model->n);

  assert(phi_old != NULL && sparse_index != NULL && sparse_count != NULL);

  printf("   --> Building sparse matrix map\n");

  /* Build sparse-matrix mapping */

  for(i=0;i<n;i++)
    {

      /* Count number of non-zero entries on this row */

      sparse_count[i] = 0;
      for(j=0;j<n;j++)
	if(fabs(model->A[i][j]) > ZERO_TOL)
	  sparse_count[i]++;

      /* Allocate space and store non-zero entries for this row */

      sparse_index[i] = (int *)calloc(sizeof(int), sparse_count[i]);

      sparse_count[i] = 0;
      for(j=0;j<n;j++)
	if(fabs(model->A[i][j]) > ZERO_TOL)
	  sparse_index[i][sparse_count[i]++] = j;
    }

  /* iterate until convergencd */

  printf("   --> Solving linear system\n");

  while(it < ITER_MAX)
    {
      /* Save last iterate */

      for(i=0;i<model->n;i++)
	phi_old[i] = model->phi[i];

      for(i=0;i<n;i++)		/* loop-over rows  */
	{
	  model->phi[i] = model->rhs[i];
	  
	  for(j=0;j<sparse_count[i];j++)
	    {
	      col_index      = sparse_index[i][j];
	      if(col_index != i)
		model->phi[i] -= model->A[i][col_index]*model->phi[col_index];
	    }
	  
	  model->phi[i] =  model->phi[i] / model->A[i][i];
	}

      it++;

      if(converged(model->phi,phi_old,CONVG_TOL,model->n,&diff)) break;

    } /* end main iterative loop */

  free(phi_old);
  for(i=0;i<model->n;i++)
    free(sparse_index[i]);
  free(sparse_index);
  
  free(sparse_count);

  printf("   --> Converged in %i iters: tolerance =  %15.7g\n", it, diff);
}

int converged(double *a, double *b, double eps, int n, double *diff)
{
  int i;
  double sum=0.;
  for (i=0; i<n; i++)
    sum += (a[i]-b[i])*(a[i]-b[i]);

  sum /= ((double)(n));

  *diff=sqrt(sum);

  if(*diff < eps)
    return(1);
  else
    return(0);
}

void solve_cg(pstruct *model)
{
  int i,j;
  double res  = 0;
  int    iter = 0;
  int n = model->n;
  double mat[n][n];

  for(i=0;i<n;i++)
      for(j=0;j<n;j++)
	mat[i][j]=model->A[i][j];
  
  printf("\n** Solving system using Conjugate Gradient\n");
  printf("   --> Solving linear system\n");
  //cg(n,*model->A,model->rhs,model->phi,&res,&iter);
  cg(n,*mat,model->rhs,model->phi,&res,&iter);
  printf("   --> Converged in %i iterations. Residual =  %15.7g\n", iter, res);

}

double compute_l2_error(pstruct *model)
{
  double l2_error = 0.0;
  double xval,yval;
  double diff;
  int i,j,index;

  for(i=0;i<model->npts;i++)
    for(j=0;j<model->npts;j++)
      {
	index = j+(i*model->npts);
	
	xval = (i)*model->h;
	yval = (j)*model->h;

	diff = masa_eval_2d_exact_phi(xval,yval)-model->phi[index];

	l2_error += diff*diff;
      }

  l2_error = sqrt(l2_error / model->n );

  return(l2_error);
}

void compute_error(pstruct *model)
{
  double error = 0.0;
  double xval, yval;
  int i,j;
  int index;

  for(i=0;i<model->npts;i++)
    {
      for(j=0;j<model->npts;j++)
	{
	  index = j+(i*model->npts);
	  
	  xval = (i)*model->h;
	  yval = (j)*model->h;
	  
	  printf("%3i: Num solution = %12.5e, Analytic solution = %12.5e (diff = %f)\n",index,
		 model->phi[index], masa_eval_2d_exact_phi(xval,yval),
		 fabs(model->phi[index]-masa_eval_2d_exact_phi(xval,yval)));
	}
      printf("\n");
    }
  return;
}
