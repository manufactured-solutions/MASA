// $License$
// $Author$
// $Id$ 

//program that takes the 2nd derivative numerically and compares to the analytical solution

#include <iostream>

using namespace std;

double a = .2;
double b = .4;
double c = .71;
double d = .1791;
double e = 1.1;
double h= .01;

double poly(double x)
{
  // analytical solution
  double pol = a*x*x + b*x + c;
  return pol;
}

double source(double x)
{
  // analytical solution
  double pol = (a/12)*x*x*x*x + (b/6)*x*x*x + (c/2)*x*x + d*x + e;
  return pol;
}

double second(double x)
{
  double pol = (source(x-h)-2*source(x) + source(x+h))/(h*h);
  return pol;
}

int main()
{
  int nx = 100;
  double x;

  for(int dx=0;dx<nx;dx++)
    {
      x=dx*h;
      cout << x << " " << second(x) << " " << poly(x) << endl;      
    }

}
