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
// $Author$
// $Id$
//
// cg.cpp: conjugate gradient solver
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <iostream>
#include <math.h>

using namespace std;

// takes two vectors, returns scalar product
double scalar(double *x1,double *x2, int n)
{
  double scal=0;
  for(int i=0;i<n;i++)
    scal += x1[i]*x2[i];
  
  return scal;
}

// return L infinity norm
double linf(double *x1,int n)
{
  double scal=0;
  for(int i=0;i<n;i++)
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
  double scal=0;
  for(int i=0;i<n;i++)
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
int cg(int n,double* A,double* b,double* x)
{
  double r [n]; // residual
  double ro[n]; // old residual
  double p[n];  // search vector
  double alpha; // step length
  double beta;  // improvement
  
  double thresh = 1e-15;
  int it = 0;
  int count = 20;
  double temp[n];  // temp vector
  
  for(int i=0;i<n;i++)
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
  //while(linf(&r[0],n)>thresh)
  while(n > it)
    {
      // calc A p_{n-1}
      for(int ii=0;ii<n;ii++)
	{
	  temp[ii]=0;
	  for(int jj=0;jj<n;jj++)
	    {
	      temp[ii] += A[ii*n+jj]*p[jj];
	    }
	}
      
      // update step length
      alpha = scalar(&ro[0],&ro[0],n)/scalar(&p[0],&temp[0],n);
      //cout << "alpha is: " << alpha << endl;

      // update approx solution & residual 
      for(int ii=0;ii<n;ii++)
	{
	  x[ii] +=  alpha*p[ii];	  
	  r[ii] = ro[ii]-alpha*temp[ii];
	}
	
      // calc improvement this step
      beta = scalar(&r[0],&r[0],n)/scalar(&ro[0],&ro[0],n);

      // update search direction
      for(int ii=0;ii<n;ii++)
	{
	  p[ii]  = r[ii] + beta * p[ii];
	  ro[ii] = r[ii]; // also update old residual
	}
      
      it++;
      //cout << r[0] << endl;
      //cout << "current thresh exceeded by " << l2(&r[0],n) << endl;

    }//done with iteration

  cout << "Performed " << it << " Iteration(s)" << endl;
  return 0;
}
