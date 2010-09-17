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
//
// masa_class.cpp: MASA class member functions and constructors
//                 Common to all MASA Manufactured Class Objects
// 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
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
  num_vars=0;                   // default -- will ++ for each registered variable
  dummy=0;
  vararr.push_back(&dummy);   // dummy used to start index at correct location
}

// define PI and other constants
const double MASA::manufactured_solution::PI = acos(-1);
const double MASA::manufactured_solution::MASA_VAR_DEFAULT = -12345.67; // default -- initialize each var to 'crazy' value

double MASA::manufactured_solution::get_var(string var)
{
  int selector=1;
  int error=1;
  
  // lets run though the list to check the variable does exist
  for(map<string,int>::const_iterator it = varmap.begin(); it != varmap.end(); ++it)
    {          
      error=var.rfind(it->first); // one value must be 0, as in equal, to exit with success
      if (error!=string::npos) // found a value
      {
	selector=0; // set flag for variable existance 
      }
    }    
  
  if(selector==0) // no error - the variable exists
    {
      selector = varmap[var];    // find location in pointer array
      return *vararr[selector];   // set to value 
      
    }
  else 
    {
      cout << "\nMASA ERROR: No such variable  (" << var << ") exists\n";
      return -20;
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
  int selector=1;
  int error=1;

  // lets run though the list to check the variable does exist
  for(map<string,int>::const_iterator it = varmap.begin(); it != varmap.end(); ++it)
    {          
      error=var.rfind(it->first); // one value must be 0, as in equal, to exit with success
      if (error!=string::npos) // found a value
      {
	selector=0; // set flag for variable existance 
      }
    }    
  
  if(selector==0)
    {      
      selector = varmap[var];    // find location in pointer array      
      *vararr[selector] = val;   // set variable to new value    
    }
  else 
    {      
      cout << "\nMASA ERROR: No such variable (" << var << ") exists to be set\n";
      return;
    } 

}// done with set_var function

void MASA::manufactured_solution::sanity_check()
{
  for(map<string,int>::const_iterator it = varmap.begin(); it != varmap.end(); ++it)
    {      
      if(*vararr[it->second] == MASA_VAR_DEFAULT)
	{
	  cout << "\nMASA WARNING: " << it->first << " has not been initialized!\n";
	}
    }    
  
  if(varmap.size() != num_vars)
    {
      cout << "\n MASA FATAL ERROR: mismatch in number of variables registered.\n"; 
      cout << "Are you calling the method manufactured_solution.register_var? This could be causing the error.\n"; 
      cout << "varmap.size() = " << varmap.size() << "; num_vars = " << num_vars << endl << endl;
      exit(1);
    }
}// done with set_var function

void MASA::manufactured_solution::register_var(string in,double* var)
{
  // first, check to ensure that no such variable has already been mapped
  if(varmap[in] <= 0) // 0 implies variable has not been registered
    {  
      // if variable has not been registered, register the variable
      num_vars++;           // this is not a bug-- we want to step num_vars up by one ONLY when adding a new variable.
      varmap[in]=num_vars;
      *var=MASA_VAR_DEFAULT;
      vararr.push_back(var);
    }
  else  // variable already registered! no unique identifier can exist!
    {
      cout << "\n MASA FATAL ERROR: \n"; 
      cout << "\n User has attempted to register two variables of the same name.\n"; 
      string error;
      return_name(&error);
      cout << " Info: error occured while constructing " << error << endl << endl;
      exit(1);      
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

MASA::masa_test::masa_test()
{
  // here, we load up the map so we can key to specific variables
  // using input
  mmsname = "masa_test_function";
  dimension = 1;

  register_var("dummy",&dummy);
  register_var("demo_var_2",&demo_var_2);
  register_var("demo_var_3",&demo_var_3);
 
}//done with constructor

// regression test blatantly stolen from paul bauman in the name of science
int MASA::manufactured_solution::poly_test()
{

  const double double_tol = 1.0e-15;

  int return_flag = 0;

  Polynomial poly;

  // a0 + a1*x + a2*x^2 + a3*x^3
  const double a0 = 1.0;
  const double a1 = 2.0;
  const double a2 = 3.0;
  const double a3 = 4.0;

  std::vector<double> a(4);
  a[0] = a0;
  a[1] = a1;
  a[2] = a2;
  a[3] = a3;

  poly.set_coeffs( a );

  // Check to make sure we get back what we set
  if( fabs( a0 - poly.get_coeffs( 0 ) ) > double_tol ) return_flag = 1;
  if( fabs( a1 - poly.get_coeffs( 1 ) ) > double_tol ) return_flag = 1;
  if( fabs( a2 - poly.get_coeffs( 2 ) ) > double_tol ) return_flag = 1;
  if( fabs( a3 - poly.get_coeffs( 3 ) ) > double_tol ) return_flag = 1;

  // Check polynomial evaluation
  const double x = 2.0;
  const double exact_value = 49.0;
  double computed_value = poly( x );
  if( fabs( exact_value - computed_value ) > double_tol ) return_flag = 1;

  // Check derivatives
  const double dx = 62;
  const double d2x = 54;
  const double d3x = 24;
  const double d4x = 0;

  std::vector<double> derivs(4);
  poly.eval_derivs( x, 4, derivs );
  
  if( fabs( exact_value - derivs[0] ) > double_tol ) return_flag = 1;
  if( fabs( dx - derivs[1] ) > double_tol ) return_flag = 1;
  if( fabs( d2x - derivs[2] ) > double_tol ) return_flag = 1;
  if( fabs( d3x - derivs[3] ) > double_tol ) return_flag = 1;
  if( fabs( d4x - derivs[4] ) > double_tol ) return_flag = 1;

  return return_flag;
}
