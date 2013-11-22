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
// $Id: 
//
// cg.c: conjugate gradient solver
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <math.h>
#include <stdio.h>

// takes two vectors, returns scalar product
double scalar(double *x1,double *x2, int n)
{
  int i;
  double scal=0;
  for(i=0;i<n;i++)
    scal += x1[i]*x2[i];
  
  return scal;
}

// return L infinity norm
double linf(double *x1,int n)
{
  int i;
  double scal=0;
  for(i=0;i<n;i++)
    {    
      //cout << "hit " << x1[i] << endl;
	  
      if(x1[i]>scal)
	{
	  scal=x1[i];
	}
    }
  return scal;
}

// return L2 norm
double l2(double *x1,int n)
{
  int i;
  double scal=0;
  for(i=0;i<n;i++)
    {    
      scal += x1[i]*x1[i];
    }
  //scal=sqrt(scal);
  return scal;
}

// j is column
// i is row
// arr[i][j] = arr[i*ncolumn + j]
//  x[i]     = x  [i]
int cg(int n,double* A,double* b,double* x,double* res,int* iter)
{
  int i,ii,j,jj;

  double r [n]; // residual
  double ro[n]; // old residual
  double p[n];  // search vector
  double alpha; // step length
  double beta;  // improvement
  
  double thresh = 1e-10;
  int it = 0;
  double temp[n];  // temp vector

  // // print solution
  // printf("\nA:\n");
  // printf("     j =");
  // for(j=0;j<n;j++)
  //   {
  //     printf("%6i ",j);
  //   }
  // printf("\n\n");

  // for(i=0;i<n;i++)
  //   {
  //     printf("i = %2i: ",i);

  //     for(j=0;j<n;j++)
  // 	{      
  // 	  printf("%6.2f ",A[i*n+j]);
  // 	}
  //     printf("\n");
  //   }

  // printf("\nb:\n");
  // for(i=0;i<n;i++)
  //   {
  //     printf("%g\n",b[i]);
  //   }
  
  for(i=0;i<n;i++)
    {
      // initialize residual to RHS
      r [i] = b[i];
      ro[i] = b[i];
      // init search vector to RHS
      p[i] = b[i];

      // init solution vector to zero
      x[i] = 0;

      temp[i] = 0;
    }

  //begin iteration
  while(linf(&r[0],n)>thresh)
    {
      // calc A p_{n-1}
      for(ii=0;ii<n;ii++)
	{
	  temp[ii]=0;
	  for(jj=0;jj<n;jj++)
	    {
	      temp[ii] += A[ii*n+jj]*p[jj];
	    }
	}
      
      // update step length
      alpha = scalar(&ro[0],&ro[0],n)/scalar(&p[0],&temp[0],n);

      // update approx solution & residual 
      for(ii=0;ii<n;ii++)
	{
	  x[ii] +=  alpha*p[ii];	  
	  r[ii] = ro[ii]-alpha*temp[ii];
	}
	
      // calc improvement this step
      beta = scalar(&r[0],&r[0],n)/scalar(&ro[0],&ro[0],n);

      // update search direction
      for(ii=0;ii<n;ii++)
	{
	  p[ii]  = r[ii] + beta * p[ii];
	  ro[ii] = r[ii]; // also update old residual
	}
      
      it++;

    }//done with iteration

  *res  = linf(&r[0],n);
  *iter = it;


  // display solution(s)
  printf("\nx:\n");
  for(i=0;i<n;i++)
    {
      printf("%g\n",x[i]);
    }

    
  return 0;
}
