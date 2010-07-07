// this is just a test routine to check the conjugate gradient solver
//
#include <iostream>
#include <math.h>
#include <masa.h>

using namespace MASA;
using namespace std;

int cg(int,double*,double*,double*);

int main()
{
  int n=3;        // size of everything
  double A[n][n]; // matrix
  double b[n];    // RHS
  double x[n];    // solution vector
   
  // init n:
  for(int i=0;i<n;i++)
    {
      /*
      for(int j=0;j<n;j++)
	{
	  A[i][j]=i*n+j+1;	  
	}
      */

      A[0][0]=1;
      A[0][1]=2;
      A[0][2]=3;

      A[1][0]=2;
      A[1][1]=1;
      A[1][2]=5;

      A[2][0]=3;
      A[2][1]=5;
      A[2][2]=1;

      b[0]=1;
      b[1]=2;
      b[2]=3;
      //x[i]=0;
    }

  // assume (n x n) (n x 1) and (n x 1)
  cg(n,&A[0][0],&b[0],&x[0]);
  
  cout << "A is: \n";
  for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
	{
	  // cout << "[ " << i << " " << j <<  " " << x[i][j] << " ] ";
	  cout << "[ " << A[i][j] << " ] ";
	}
      cout << endl;
    }

  cout << "\n\nb is: \n";
  for(int i=0;i<n;i++)
    {
      cout << "[ " << b[i] << " ]\n";
    }

  cout << "\n\nx is: \n";
  for(int i=0;i<n;i++)
    {
      cout << "[ " << x[i] << " ]\n";
    }

}// end program
