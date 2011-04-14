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
#include <limits>
#include <assert.h>

using namespace MASA;

/* ------------------------------------------------
 *
 *         Manufactured Solution Class 
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
Scalar MASA::manufactured_solution<Scalar>::pow(long double x, double y)
{
  Scalar tx = Scalar(x);
  Scalar ty = Scalar(y);

  return std::pow(tx,ty);
}

/* template <typename Scalar>
Scalar MASA::manufactured_solution<Scalar>::pow(double x, long double y)
{
  Scalar tx = Scalar(x);
  Scalar ty = Scalar(y);

  return std::pow(tx,ty);
}

template <typename Scalar>
Scalar MASA::manufactured_solution<Scalar>::pow(long double x, long double y)
{
  Scalar tx = Scalar(x);
  Scalar ty = Scalar(y);

  return std::pow(tx,ty);
}

template <typename Scalar>
Scalar MASA::manufactured_solution<Scalar>::pow(double x, double y)
{
  Scalar tx = Scalar(x);
  Scalar ty = Scalar(y);

  return std::pow(tx,ty);
} */

template <typename Scalar>
MASA::manufactured_solution<Scalar>::manufactured_solution()
{  
  std::vector<Scalar> dumvec;

  num_vars=0;                   // default -- will ++ for each registered variable
  dummy=0;
  dumvec.resize(2);
  vararr.push_back(&dummy);   // dummy used to start index at correct location
  vecarr.push_back(&dumvec);   // dummy used to start index at correct location
  }

// define PI and other constants
template <typename Scalar>
const Scalar MASA::manufactured_solution<Scalar>::PI = acos(Scalar(-1));

template <typename Scalar>
const Scalar MASA::manufactured_solution<Scalar>::pi = acos(Scalar(-1));

template <typename Scalar>
const Scalar MASA::manufactured_solution<Scalar>::MASA_VAR_DEFAULT = -12345.67; // default init each var to 'crazy' val

template <typename Scalar>
int MASA::manufactured_solution<Scalar>::register_vec(std::string in,std::vector<Scalar>& vec)
{

  // working fine here
  //std::cout << (*vec).size() << std::endl;

  // first, check to ensure that no such vec has already been mapped
  if(vecmap[in] <= 0) // 0 implies variable has not been registered
    {  
      // if vec has not been registered, register the vec
      num_vec++;           // we want to step num_vars up by one ONLY when adding a new vec.
      vecmap[in]=num_vec;
      vecarr.push_back(&vec);
    }
  else  // variable already registered! no unique identifier can exist!
    {
      std::cout << "\n MASA FATAL ERROR:: \n"; 
      std::cout << "\n Attempted to register two vectors of the same name.\n"; 
      std::string error;
      return_name(&error);
      std::cout << " Info: error occured while constructing " << error << std::endl << std::endl;
      return 1;
    }

  return 0; // smooth sailing
  
}// done with register_var function

template <typename Scalar>
int MASA::manufactured_solution<Scalar>::get_vec(std::string name,std::vector<Scalar>& vec)
{
  std::map<std::string,int>::const_iterator selector;
  
  // find variable
  selector = vecmap.find(name);
  
  // error handling
  if(selector == vecmap.end())
    {
      std::cout << "\nMASA ERROR!!!:: No such variable  (" << name << ") exists\n";
      return 1;
    }
  
  //std::cout << "this is: " << vec.size();
  vec = *vecarr[selector->second];   // set to value 
  //std::cout << "this is: " << vecarr[(*selector).second]->size();
  //std::cout << "this is: " << vec.size();
  return 0;

}

template <typename Scalar>
int MASA::manufactured_solution<Scalar>::set_vec(std::string var,std::vector<Scalar>& vec)
{
  //std::cout << "vec is: " << (*vec)[0] << std::endl;
  //std::cout << "vecarr is: " << vecarr[selector->second] << std::endl;

  std::map<std::string,int>::const_iterator selector;

  // find variable
  selector = vecmap.find(var);
  
  // error handling
 if(selector == vecmap.end())
    {
      std::cout << "\nMASA ERROR!!!:: No such array  (" << var << ") exists to be set\n";
      return 1;
    }
 
 // fix vector to same size and values as new guy
 *vecarr[selector->second] = vec;
  return 0; // exit with no error
 
}// done with set_vec function

template <typename Scalar>
int MASA::manufactured_solution<Scalar>::display_vec()
{
  std::vector<Scalar>* vec;
  
  std::cout << "\nMASA :: Solution has " << vecmap.size() << " vector(s).\n";
  std::cout << "*-------------------------------------*\n" ;

  for(std::map<std::string,int>::const_iterator it = vecmap.begin(); it != vecmap.end(); ++it)
    {      
      vec = vecarr[it->second];
      std::cout << it->first <<" is size: " << vec->size() << '\n';      
    }

  std::cout << "*-------------------------------------*\n" ;

  return 0;

}

template <typename Scalar>
Scalar MASA::manufactured_solution<Scalar>::get_var(std::string var)
{
  std::map<std::string,int>::const_iterator selector;
  
  // find variable
  selector = varmap.find(var);
  
  // error handling
  if(selector == varmap.end())
    {
      std::cout << "\nMASA ERROR!!!:: No such variable  (" << var << ") exists\n";
      return -20;
    }
  
  return *vararr[(*selector).second];   // set to value 
  
}// done with get_var function


template <typename Scalar>
int MASA::manufactured_solution<Scalar>::display_var()
{
  Scalar threshold = 5 * std::numeric_limits<Scalar>::epsilon();

  std::cout << "\nMASA :: Solution has " << varmap.size() << " variables.\n";
  std::cout << "*-------------------------------------*\n" ;

  for(std::map<std::string,int>::const_iterator it = varmap.begin(); it != varmap.end(); ++it)
    {   
      // adding conditional to avoid confusing our users about uninitalized variables
      // this is because the default is a bit odd... -12345.7 might appear 'set'

      if((*vararr[it->second] - MASA_VAR_DEFAULT) <= threshold)
	{
	  std::cout << it->first <<" is set to: Uninitialized\n";
	}
      else //value has been set
	{
	  std::cout.precision(16);
	  std::cout << it->first <<" is set to: " << *vararr[it->second] << '\n';
	}
      
    }    

  std::cout << "*-------------------------------------*\n" ;
  
  return 0;

} // done with display all variable names

template <typename Scalar>
int MASA::manufactured_solution<Scalar>::set_var(std::string var, Scalar val)
{
  std::map<std::string,int>::const_iterator selector;
  
  // find variable
  selector = varmap.find(var);
  
  // error handling
  if(selector == varmap.end())
    {
      std::cout << "\nMASA ERROR!!!:: No such variable  (" << var << ") exists to be set\n";
      return 1;
    }
  
  // set new value
  *vararr[(*selector).second] = val;
  return 0; // exit with no error

}// done with set_var function

template <typename Scalar>
int MASA::manufactured_solution<Scalar>::purge_var()
{
  // MASA_VAR_DEFAULT
  for(std::map<std::string,int>::const_iterator it = varmap.begin(); it != varmap.end(); ++it)
    {      
      *vararr[it->second]=MASA_VAR_DEFAULT;      
    }
  return 0;
}// done with purge_var function


template <typename Scalar>
Scalar MASA::manufactured_solution<Scalar>::pass_function(Scalar (*in_func)(Scalar),Scalar a)
{

  // just want to evaluate function here
  Scalar out = in_func(a);
  return out;
}

template <typename Scalar>
int MASA::manufactured_solution<Scalar>::sanity_check()
{
  int flag=0;
  Scalar thresh = 1.0e-10;

  // check all scalar values
  for(std::map<std::string,int>::const_iterator it = varmap.begin(); it != varmap.end(); ++it)
    {      
      if(fabs((*vararr[it->second] - MASA_VAR_DEFAULT)/MASA_VAR_DEFAULT) < thresh)
	{
	  std::cout << "\nMASA WARNING:: " << it->first << " has not been initialized!\n";
	  std::cout << "Current value is: " << fabs((*vararr[it->second] - MASA_VAR_DEFAULT)/MASA_VAR_DEFAULT) << std::endl;
	  flag += 1;
	}
    }    
  
  if(int(varmap.size()) != num_vars)
    {
      std::cout << "\n MASA FATAL ERROR:: mismatch in number of variables registered.\n"; 
      std::cout << "Are you calling the method manufactured_solution.register_var? This could be causing the error.\n"; 
      std::cout << "varmap.size() = " << varmap.size() << "; num_vars = " << num_vars << std::endl << std::endl;
      masa_exit(1);
    }

  
  // check all vector values
  for(std::map<std::string,int>::const_iterator it = vecmap.begin(); it != vecmap.end(); ++it)
    {    
      std::vector<Scalar>* vec;
      vec = vecarr[it->second];
	
      if((*vec).size() == 0)
	{
	  std::cout << "\nMASA WARNING:: vector " << it->first << " has not been initialized!\n";
	  flag += 1;
	}
    }    
  
  if(flag != 0)
    {
      return 1; // not all values init
    }
  else
    {
      return 0; // all values init: smooth sailing
    }
  

}// done with ssanity check

template <typename Scalar>
int MASA::manufactured_solution<Scalar>::register_var(std::string in,Scalar* var)
{
  // first, check to ensure that no such variable has already been mapped
  if(varmap[in] <= 0) // 0 implies variable has not been registered
    {  
      // if variable has not been registered, register the variable
      num_vars++;           // we want to step num_vars up by one ONLY when adding a new variable.
      varmap[in]=num_vars;
      *var=MASA_VAR_DEFAULT;
      vararr.push_back(var);
    }
  else  // variable already registered! no unique identifier can exist!
    {
      std::cout << "\n MASA FATAL ERROR:: \n"; 
      std::cout << "\n Attempted to register two variables of the same name.\n"; 
      std::string error;
      return_name(&error);
      std::cout << " Info: error occured while constructing " << error << std::endl << std::endl;
      return 1;
    }

  return 0; // smooth sailing
  
}// done with register_var function

/* ------------------------------------------------
 *
 *         Polynomial Class
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
void Polynomial<Scalar>::set_coeffs( const std::vector<Scalar> &coeffs_in )
{
  int num_coeffs = coeffs_in.size();

  coeffs.resize( num_coeffs );

  coeffs = coeffs_in;

  return;
}

template <typename Scalar>
Scalar Polynomial<Scalar>::operator()( const Scalar &x, int *err) const
{  
  int num_coeffs = coeffs.size();
  if(num_coeffs == 0)
    {
      *err=1;
      return 0;
    }

  int n = num_coeffs-1;
  Scalar y;

  // We use Horner's method here. 
  y = coeffs[n];

  for( int i = n-1; i >= 0; --i )
    {
      y = coeffs[i] + y*x;
    }

  *err = 0; // no errors
  return y;

}

template <typename Scalar>
void Polynomial<Scalar>::eval_derivs( const Scalar x, const int k, std::vector<Scalar> & derivs ) const
{

  // Zero out the vector first.
  for( int i = 0; i < k; ++i) derivs[i] = 0.0;

  int num_coeffs = coeffs.size();

  int n = num_coeffs-1;

  // We use Horner's method here.
  // Can use proof by induction to prove recursion formula for derivatives.
  derivs[0] = coeffs[n];

  for( int i = n-1; i >= 0; --i )
    {

      for( int j = k-1; j > 0; --j)
	{
	  derivs[j] = j*derivs[j-1] + x*derivs[j]; 
	}

      derivs[0] = coeffs[i] + x*derivs[0];
      
    }
  
  return;
}

template <typename Scalar>
Scalar Polynomial<Scalar>::get_coeffs( const int &coeff_index ) const
{
  assert( coeff_index >= 0 );
  assert( coeff_index <= (int(coeffs.size())-1) );
  return coeffs[coeff_index];
}

/* ------------------------------------------------
 *
 *         Test Problem
 *
 * -----------------------------------------------
 */ 

template <typename Scalar>
MASA::masa_test_function<Scalar>::masa_test_function()
{
  // here, we load up the map so we can key to specific variables
  // using input
  this->mmsname = "masa_test_function";
  this->dimension = 1;

  // to do 
  // WARNING: this is designed to fail! 
  // This function tests the MASA error handling for: 
  // registering two variables of the same name
  register_var("demo_var_2",&demo_var_2);
  register_var("demo_var_3",&demo_var_3);
 
}//done with constructor

// regression test blatantly stolen from paul bauman in the name of science
template <typename Scalar>
int MASA::manufactured_solution<Scalar>::poly_test()
{
  Scalar thresh = 5 * std::numeric_limits<Scalar>::epsilon();
  
  int return_flag = 0;
  int ierr = -1;

  Polynomial<Scalar> poly;

  // a0 + a1*x + a2*x^2 + a3*x^3
  const Scalar a0 = 1.0;
  const Scalar a1 = 2.0;
  const Scalar a2 = 3.0;
  const Scalar a3 = 4.0;

  std::vector<Scalar> a(4);
  a[0] = a0;
  a[1] = a1;
  a[2] = a2;
  a[3] = a3;

  // check poly will return failure (1) if no coeff set
  Scalar computed_value = poly( 0 , &ierr);
  if(ierr != 1) return_flag=1;

  poly.set_coeffs( a );

  // Check to make sure we get back what we set
  if( fabs( a0 - poly.get_coeffs( 0 ) ) > thresh ) return_flag = 1;
  if( fabs( a1 - poly.get_coeffs( 1 ) ) > thresh ) return_flag = 1;
  if( fabs( a2 - poly.get_coeffs( 2 ) ) > thresh ) return_flag = 1;
  if( fabs( a3 - poly.get_coeffs( 3 ) ) > thresh ) return_flag = 1;

  // Check polynomial evaluation
  const Scalar x = 2.0;
  const Scalar exact_value = 49.0;
  
  // evaluate and check for 'good' return value (0)
  computed_value = poly( x , &ierr);
  if(ierr != 0) return_flag=1;

  if( fabs( exact_value - computed_value ) > thresh ) return_flag = 1;

  // Check derivatives
  const Scalar dx = 62;
  const Scalar d2x = 54;
  const Scalar d3x = 24;

  std::vector<Scalar> derivs(4);

  poly.eval_derivs( x, 4, derivs );
  
  if( fabs( exact_value - derivs[0] ) > thresh ) return_flag = 1;
  if( fabs( dx  - derivs[1] ) > thresh ) return_flag = 1;
  if( fabs( d2x - derivs[2] ) > thresh ) return_flag = 1;
  if( fabs( d3x - derivs[3] ) > thresh ) return_flag = 1;

  return return_flag;
}

template <typename Scalar>
int MASA::masa_test_function<Scalar>::init_var()
{
  int err = 0;

  // designed to fail -- 2nd var does not exist
  // 3rd var has already been registered
  err += this->set_var("demo_var_2",1);   
  err += this->set_var("demo_var_12",1);   
  err += this->register_var("demo_var_3",&demo_var_3);

  // now for a really epic fail: user sets var 
  // array instead of calling register_var method
  //vararr.push_back(&demo_var_3);
  this->num_vars++;

  return err;

}

template <typename Scalar>
MASA::masa_uninit<Scalar>::masa_uninit()
{
  // nothing to see here
  this->mmsname = "masa_uninit";
  this->dimension = 1;
}

template <typename Scalar>
int MASA::masa_uninit<Scalar>::init_var()
{

  return 0;

}

// ----------------------------------------
//   Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::manufactured_solution);
MASA_INSTANTIATE_ALL(MASA::masa_test_function);
MASA_INSTANTIATE_ALL(MASA::masa_uninit);
MASA_INSTANTIATE_ALL(MASA::Polynomial);
