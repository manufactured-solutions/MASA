 /*--------------------------------------------------------------------------
  *--------------------------------------------------------------------------
  *
  * Copyright (C) 2010 The PECOS Development Team
  *
  * Please see http://pecos.ices.utexas.edu for more information.
  *
  * This file is part of MASA.
  *
  * MASA is free software: you can redistribute it and/or modify it under
  * the terms of the GNU Lesser General Public License as published by the Free
  * Software Foundation, either version 3 of the License, or (at your option)
  * any later version.
  *
  * MASA is distributed in the hope that it will be useful, but WITHOUT ANY
  * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
  * details.
  *
  * You should have received a copy of the GNU Lesser General Public License along 
  * with MASA.  If not, see <http://www.gnu.org/licenses/>.
  *
  *--------------------------------------------------------------------------
  
  MASA -- Manufactured Analytical Solutions Abstraction Library

  A software interface that provides access to all manufactured solutions to 
  be used by various models throughout the center.
  
  *--------------------------------------------------------------------------
  */  

//
//   These are the MASA class member functions and constructors
// 

#include <masa_internal.h>
#include <assert.h>

using namespace MASA;

/* ------------------------------------------------
 *
 *         Manufactured Solution Class 
 *
 * -----------------------------------------------
 */ 
MASA::manufactured_solution::manufactured_solution()
{
  MASA_VAR_DEFAULT = -12345.67; // default --initialize each var to 'crazy' value
}
void MASA::manufactured_solution::get_var(string var, double* sol)
{
  int selector;
  selector = varmap[var];    // find location in pointer array

    if(selector==0)
    {
      cout << "\nMASA ERROR: No such variable exists\n";
    }
    else 
      {
	*sol = *vararr[selector];   // set to value 
      } 
    
}// done with get_var function

void MASA::manufactured_solution::display_var()
{
  int selector;
  cout << "\n Solution has " << varmap.size() << " variables.\n\n";
  for(map<string,int>::const_iterator it = varmap.begin(); it != varmap.end(); ++it)
    {      
      cout << it->first <<" is set to: " << *vararr[it->second] << '\n';
    }    
  
} // done with display all variable names

void MASA::manufactured_solution::set_var(string var, double val)
{
  int selector;
  selector = varmap[var];    // find location in pointer array

  if(selector==0)
    {
      cout << "\nMASA ERROR: No such variable to be set\n";
      // exit(1);not really a fatal error
    }
  else 
    {
      *vararr[selector] = val;   // set variable to new value    
    } 
}// done with set_var function

void MASA::manufactured_solution::sanity_check()
{
  for(map<string,int>::const_iterator it = varmap.begin(); it != varmap.end(); ++it)
    {      
      if(*vararr[it->second] == MASA_VAR_DEFAULT)
	{
	  cout << "\nMASA WARNING: " << it->first << " is not initialized!\n";
	}
    }    

}// done with set_var function

/* ------------------------------------------------
 *
 *         Polynomial Class
 *
 * -----------------------------------------------
 */ 

void Polynomial::set_coeffs( const std::vector<double> &coeffs_in )
{
  int num_coeffs = coeffs_in.size();

  coeffs.resize( num_coeffs );

  coeffs = coeffs_in;

  return;
}

double Polynomial::operator()( const double &x ) const
{
  int num_coeffs = coeffs.size();

  int n = num_coeffs-1;

  double y;

  // We use Horner's method here. 
  y = coeffs[n];

  for( int i = n-1; i >= 0; --i )
    {
      y = coeffs[i] + y*x;
    }

  return y;

}

void Polynomial::eval_derivs( const double &x, const int & k, std::vector<double> & derivs ) const
{

  // Zero out the vector first.
  for( int i = 0; i <= k; ++i) derivs[i] = 0.0;

  int num_coeffs = coeffs.size();

  int n = num_coeffs-1;

  // We use Horner's method here.
  // Can use proof by induction to prove recursion formula for derivatives.
  derivs[0] = coeffs[n];

  for( int i = n-1; i >= 0; --i )
    {

      for( int j = k; j > 0; --j)
	{
	  derivs[j] = j*derivs[j-1] + x*derivs[j]; 
	}

      derivs[0] = coeffs[i] + x*derivs[0];
      
    }
  
  return;
}

double Polynomial::get_coeffs( const int &coeff_index ) const
{
  assert( coeff_index >= 0 );
  assert( coeff_index <= (coeffs.size()-1) );
  return coeffs[coeff_index];
}

/* ------------------------------------------------
 *
 *         Test Problem
 *
 * -----------------------------------------------
 */ 

MASA::MASA_Test::MASA_Test()
{
  // here, we load up the map so we can key to specific variables
  // using input
  mmsname = "MASA_test_function";
  dimension = 1;
  
  //first variable "axp" -- load map and array
  varmap["axp"]=1;
  axp=MASA_VAR_DEFAULT;
  vararr.push_back(&axp);
  vararr.push_back(&axp);

  varmap["demo_var_2"]=2;
  demo_var_2=MASA_VAR_DEFAULT;
  vararr.push_back(&demo_var_2);

  varmap["demo_var_3"]=3;
  demo_var_3=MASA_VAR_DEFAULT;
  vararr.push_back(&demo_var_3);
  
}//done with constructor

double MASA::MASA_Test::eval_q_u(double x)
{
  double qt = demo_var_2 + demo_var_3 + axp;
  return qt;

}//done with constructor


/* ------------------------------------------------
 *
 *         Compressible Navier-Stokes EQUATIONs
 *
 *
 *
 * -----------------------------------------------
 */ 

MASA::ns_compress_2d::ns_compress_2d()
{
  mmsname = "navierstokes_compressible_2d";
    dimension=2;

}//done with constructor

MASA::ns_compress_3d::ns_compress_3d()
{
  mmsname = "navierstokes_compressible_3d";
    dimension=3;

}//done with constructor
