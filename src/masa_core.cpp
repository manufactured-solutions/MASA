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
// $Author$
// $Id$
//
// masa_core.cpp: this is the core set of routines -- the only functions that 
//                should be called by users -- all other files are internal
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <config.h>        // for MASA_EXCEPTIONS conditional
#include <smasa.h>
#include <map>

using namespace MASA;
using namespace std;

// Anonymous namespace for local helper class/functions
namespace {

using namespace MASA;

template <typename Scalar>
class MasterMS
{
public:
  MasterMS () : _master_pointer(NULL), _master_map() {}

  ~MasterMS () {
    if (!_master_map.empty()) // workaround for icpc 12.1.6 map bug
      for(typename map<std::string,manufactured_solution<Scalar>*>::iterator iter = this->_master_map.begin(); iter != this->_master_map.end(); iter++)
        delete iter->second;

    // workaround for icpc 12.1.6 "double destruct globals" bug
    _master_map.clear();
  }

  const manufactured_solution<Scalar>& get_ms() const {
    verify_pointer_sanity();
    return *_master_pointer;
  }

  manufactured_solution<Scalar>& get_ms() {
    verify_pointer_sanity();
    return *_master_pointer;
  }

  void select_mms(const std::string& my_name);

  void init_mms(const std::string& my_name, const std::string& masa_name);

  void list_mms() const;

  unsigned int size() const { return _master_map.size(); }

private:
  //
  //  this function checks the user has an active mms
  //
  void verify_pointer_sanity() const;

  manufactured_solution<Scalar>*             _master_pointer; // pointer to currently selected manufactured solution
  std::map<std::string, manufactured_solution<Scalar> *> _master_map; // global map b/t unique name and manuf class
};


template <typename Scalar>
int get_list_mms(std::vector<manufactured_solution<Scalar>*>& anim)
{
  // Build a temporary vector of MMS objects, then sort them into our map by name
  anim.push_back(new masa_test_function<Scalar>());   // test function
  anim.push_back(new masa_uninit<Scalar>()); // another test function
  
  //  ** register solutions here - lets keep this alphabetical ** 

  // axisymmetric solutions
  anim.push_back(new axi_cns<Scalar>());
  anim.push_back(new axi_euler<Scalar>());

  // SMASA::
  anim.push_back(new cp_normal<Scalar>());

  // euler 
  anim.push_back(new euler_1d<Scalar>());
  anim.push_back(new euler_2d<Scalar>());
  anim.push_back(new euler_3d<Scalar>());
  anim.push_back(new euler_transient_1d<Scalar>());
  anim.push_back(new euler_transient_2d<Scalar>());
  anim.push_back(new euler_transient_3d<Scalar>());

  anim.push_back(new euler_chem_1d<Scalar>());

  // favre averaged navier stokes
  anim.push_back(new fans_sa_steady_wall_bounded<Scalar>());
  anim.push_back(new fans_sa_transient_free_shear<Scalar>());

  // heat equation
  anim.push_back(new heateq_1d_steady_const<Scalar>());
  anim.push_back(new heateq_2d_steady_const<Scalar>());
  anim.push_back(new heateq_3d_steady_const<Scalar>());

  anim.push_back(new heateq_1d_steady_var<Scalar>());
  anim.push_back(new heateq_2d_steady_var<Scalar>());
  anim.push_back(new heateq_3d_steady_var<Scalar>());

  anim.push_back(new heateq_1d_unsteady_const<Scalar>());
  anim.push_back(new heateq_2d_unsteady_const<Scalar>());
  anim.push_back(new heateq_3d_unsteady_const<Scalar>());

  anim.push_back(new heateq_1d_unsteady_var<Scalar>());
  anim.push_back(new heateq_2d_unsteady_var<Scalar>());
  anim.push_back(new heateq_3d_unsteady_var<Scalar>());

  // laplacian
  anim.push_back(new laplace_2d<Scalar>());

  // navier stokes
  anim.push_back(new navierstokes_2d_compressible<Scalar>());
  anim.push_back(new navierstokes_3d_compressible<Scalar>());
  anim.push_back(new navierstokes_4d_compressible_powerlaw<Scalar>());
  anim.push_back(new navierstokes_ablation_1d_steady<Scalar>());

  // radiation
  anim.push_back(new radiation_integrated_intensity<Scalar>());

  // reynolds averaged navier stokes
  anim.push_back(new rans_sa<Scalar>());

  // sod shock tube
  anim.push_back(new sod_1d<Scalar>());

  // MetaPhysicL-based solutions
#ifdef HAVE_METAPHYSICL
  anim.push_back(new ad_cns_2d_crossterms<Scalar>());
  anim.push_back(new ad_cns_3d_crossterms<Scalar>());
  anim.push_back(new convdiff_steady_nosource_1d<Scalar>());
  anim.push_back(new navierstokes_3d_incompressible<Scalar>());
  anim.push_back(new navierstokes_3d_incompressible_homogeneous<Scalar>());
  anim.push_back(new navierstokes_3d_transient_sutherland<Scalar>());
#endif // HAVE_METAPHYSICL

  // automatically generated MMS:

  anim.push_back(new burgers_equation<Scalar>());
  anim.push_back(new axi_euler_transient<Scalar>());
  anim.push_back(new axi_cns_transient<Scalar>());


  // --l33t-- DO NOT EDIT THIS LINE OR ANY BELOW IT


  return 0;
}

// Instantiations for every precision

MasterMS<double>      masa_master_double;
MasterMS<long double> masa_master_longdouble;

// Function to return a MasterMS by precision
template <typename Scalar>
MasterMS<Scalar>&      masa_master() { return masa_master_double; }
template <>
MasterMS<long double>& masa_master() { return masa_master_longdouble; }

}

template <typename Scalar>
void MasterMS<Scalar>::verify_pointer_sanity() const
{
  if(_master_pointer == 0)
    {    
      std::cout << "MASA FATAL ERROR:: No initialized Manufactured Solution!" << std::endl;
      std::cout << "Have you called masa_init?" << std::endl;
      masa_exit(1);
    }  
}

//
//  limited masa exception handling
//
void MASA::masa_exit(int ex)
{

#ifdef MASA_EXCEPTIONS
  std::cout << "MASA:: caught exception " << ex << std::endl;
  throw(ex);
#else
  std::cout << "MASA:: ABORTING\n";
  exit(ex);
#endif
  
}

template <typename Scalar>
int MASA::masa_purge_default_param()
{
  return masa_master<Scalar>().get_ms().purge_var();
}

template <typename Scalar>
Scalar MASA::pass_func(Scalar (*in_func)(Scalar),Scalar a)
{
  return masa_master<Scalar>().get_ms().pass_function(in_func,a);
}

//
//  this function selects an already initialized manufactured class
//
template <typename Scalar>
void MasterMS<Scalar>::select_mms(const std::string& my_name)
{
  // check that the class does exist
  typename std::map<std::string, manufactured_solution<Scalar> *>::iterator it=_master_map.find(my_name);
  if(it != _master_map.end()) // found a name
    { 
      std::cout << "MASA :: selected " << my_name << std::endl;
      _master_pointer=it->second; // set pointer to currently selected solution      
    }      
  else 
    {
      std::cout << "\nMASA FATAL ERROR:: No such manufactured solution (" << my_name << ") has been initialized.\n";
      this->list_mms();
      masa_exit(1);
    } 
}


//
//  this function tests if a value is set to default
//
template <typename Scalar>
int MASA::masa_test_default(Scalar input)
{
  Scalar MASA_VAR_DEFAULT = -12345.67;
  Scalar uninit = -1.33;
  Scalar thresh = 5 * std::numeric_limits<Scalar>::epsilon();

  if( fabs((input - MASA_VAR_DEFAULT)/MASA_VAR_DEFAULT) < thresh)
    {
      exit(1);
    }

  if(fabs((input - uninit)/uninit) < thresh)
    {
      exit(1);
    }
  
  exit(0);
}

template <typename Scalar>
int MASA::masa_select_mms(std::string name)
{
  masa_master<Scalar>().select_mms(name);
  return 0;
}

//
//  this function will initiate a masa manufactured class
//
template <typename Scalar>
int MASA::masa_init(std::string unique_name, std::string str)
{
  masa_master<Scalar>().init_mms(unique_name, str);
  
  return 0; // steady as she goes
}


template <typename Scalar>
void MasterMS<Scalar>::init_mms(const std::string& my_name,
                                const std::string& masa_name)
{
  std::vector<manufactured_solution<Scalar>*> anim;
  get_list_mms<Scalar>(anim); //construct maps of MMS objects

  std::string mapped_name = masa_name;
  MASA::masa_map(&mapped_name);

  for (unsigned int i=0; i != anim.size(); ++i)
    {
      std::string name;
      anim[i]->return_name(&name);
      if (name.empty())
        {
          std::cout << "MASA FATAL ERROR:: manufactured solution has no name!\n";
          masa_exit(1);
        }
      if (name == mapped_name)
        {
          _master_map[my_name] = _master_pointer = anim[i];
          return;
        }
      else
        delete anim[i];
    }

  std::cout << "MASA FATAL ERROR:: no manufactured solution named " << masa_name << " found!\n";
  masa_exit(1);
}


template <typename Scalar>
void MasterMS<Scalar>::list_mms() const
{
  std::string str;

  // output the size of the map
  std::cout << "Number of initialized solutions: " << this->size() << std::endl;
  for(typename map<std::string,manufactured_solution<Scalar>*>::const_iterator iter = this->_master_map.begin(); iter != this->_master_map.end(); iter++)
    {
      (iter->second)->return_name(&str);
      std::cout << iter->first << " : " << str << std::endl;
    }
}



template <typename Scalar>
int MASA::masa_list_mms()
{
  masa_master<Scalar>().list_mms();
  return 0;
}

//
// function that prints all registered masa solutions
//
template <typename Scalar>
int MASA::masa_printid()
{
  std::vector<manufactured_solution<Scalar>*> anim;

  get_list_mms(anim); //construct list 

  std::cout << std::endl;
  std::cout << "\nMASA :: Available Solutions:\n";
  std::cout << "*-------------------------------------*" ;

  for (typename std::vector<manufactured_solution<Scalar>*>::const_iterator it = anim.begin(); it != anim.end(); ++it) 
    {
      std::string name;
      (*it)->return_name(&name); // get name
      std::cout << std::endl << name;
      delete *it;
    } // done with for loop 

  std::cout << "\n*-------------------------------------*\n" ;

  return 0; // steady as she goes
}// done with masa print id

template <typename Scalar>
void MASA::masa_set_param(std::string param,Scalar paramval)
{
  masa_master<Scalar>().get_ms().set_var(param,paramval);
}

template <typename Scalar>
void MASA::masa_set_vec(std::string vector_name,std::vector<Scalar>& new_vector)
{
  masa_master<Scalar>().get_ms().set_vec(vector_name,new_vector);
}

//
// Set all parameters to default values
//
template <typename Scalar>
int MASA::masa_init_param()
{
  return masa_master<Scalar>().get_ms().init_var();
}

//
// Function that returns value of parameter selected by string
// 

template <typename Scalar>
Scalar MASA::masa_get_param(std::string param)
{
  return masa_master<Scalar>().get_ms().get_var(param);
}

//
// Function that returns vector -- selected by string
// 

template <typename Scalar>
int MASA::masa_get_vec(std::string vector_name,std::vector<Scalar>& vector)
{
  return masa_master<Scalar>().get_ms().get_vec(vector_name,vector);
}


template <typename Scalar>
int MASA::masa_display_param()
{
  return masa_master<Scalar>().get_ms().display_var();
}

template <typename Scalar>
int MASA::masa_display_vec()
{
  return masa_master<Scalar>().get_ms().display_vec();
}

/* ------------------------------------------------
 *
 *         1D functions
 *
 * -----------------------------------------------
 */ 


  // --------------------------------
  // source terms
  // --------------------------------

template <typename Scalar>
Scalar MASA::masa_eval_source_t(Scalar x) //x 
{
  return masa_master<Scalar>().get_ms().eval_q_t(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_t(Scalar x,Scalar t) //x,t
{
  return masa_master<Scalar>().get_ms().eval_q_t(x,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_u(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_u(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_v(Scalar x)  // for SA model
{
  return masa_master<Scalar>().get_ms().eval_q_v(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_w(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_w(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_rho(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_u(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_u(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_v(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_v(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_w(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_w(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_e(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_e(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_boundary(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_u_boundary(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_N(Scalar x,Scalar (*in_func)(Scalar))
{
  return masa_master<Scalar>().get_ms().eval_q_rho_N(x,in_func);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_N2(Scalar x,Scalar (*in_func)(Scalar))
{
  return masa_master<Scalar>().get_ms().eval_q_rho_N2(x,in_func);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_C(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_C(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_C3(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_C3(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_C(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_C(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_C3(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_C3(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_e(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_q_e(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_e(Scalar x,Scalar (*in_func)(Scalar))
{
  return masa_master<Scalar>().get_ms().eval_q_e(x,in_func);
}

  // --------------------------------
  // analytical terms
  // --------------------------------

template <typename Scalar>
Scalar MASA::masa_eval_exact_t(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_exact_t(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_u(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_exact_u(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_v(Scalar x) // for SA model
{
  return masa_master<Scalar>().get_ms().eval_exact_v(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_w(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_exact_w(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_p(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_exact_p(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho_N(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho_N(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho_N2(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho_N2(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho_C(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho_C(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho_C3(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho_C3(x);
}

// --------------------------------
// smasa: 1D
// --------------------------------

template <typename Scalar>
Scalar MASA::masa_eval_likelyhood(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_likelyhood(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_loglikelyhood(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_loglikelyhood(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_prior(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_prior(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_posterior(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_posterior(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_central_moment(int x)
{
  return masa_master<Scalar>().get_ms().eval_cen_mom(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_posterior_mean()
{
  return masa_master<Scalar>().get_ms().eval_post_mean();
}

template <typename Scalar>
Scalar MASA::masa_eval_posterior_variance()
{
  return masa_master<Scalar>().get_ms().eval_post_var();
}

// --------------------------------
// gradients: 1D
// --------------------------------

template <typename Scalar>
Scalar MASA::masa_eval_grad_t(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_g_u(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_u(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_g_u(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_v(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_g_v(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_w(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_g_w(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_p(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_g_p(x);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_rho(Scalar x)
{
  return masa_master<Scalar>().get_ms().eval_g_rho(x);
}


/* ------------------------------------------------
 *
 *         2D functions
 *
 * -----------------------------------------------
 */ 

  // --------------------------------
  // source terms
  // --------------------------------

template <typename Scalar>
Scalar MASA::masa_eval_source_t(Scalar x,Scalar y,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_q_t(x,y,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_f(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_q_f(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_u(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_q_u(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_v(Scalar x,Scalar y)
{  
  return masa_master<Scalar>().get_ms().eval_q_v(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_w(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_q_w(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_q_rho(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_e(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_q_e(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_u(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_u(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_v(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_v(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_w(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_w(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_e(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_e(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_nu(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_q_nu(x,y);
}


  // --------------------------------
  // analytical terms
  // --------------------------------

template <typename Scalar>
Scalar MASA::masa_eval_exact_t(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_exact_t(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_u(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_exact_u(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_phi(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_exact_phi(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_v(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_exact_v(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_w(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_exact_w(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_p(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_exact_p(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_nu(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_exact_nu(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho_C(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho_C(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho_C3(Scalar x,Scalar y)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho_C3(x,y);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_t(Scalar x,Scalar y,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_t(x,y,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_u(Scalar x,Scalar y,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_u(x,y,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_v(Scalar x,Scalar y,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_v(x,y,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_w(Scalar x,Scalar y,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_w(x,y,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_p(Scalar x,Scalar y,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_p(x,y,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_rho(Scalar x,Scalar y,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_rho(x,y,i);
}

/* ------------------------------------------------
 *
 *         3D functions
 *
 * -----------------------------------------------
 */ 

  // --------------------------------
  // source terms
  // --------------------------------

template <typename Scalar>
Scalar MASA::masa_eval_source_t(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_q_t(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_u(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_q_u(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_u(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_q_u(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_v(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_q_v(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_v(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_q_v(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_w(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_q_w(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_w(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_q_w(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho(Scalar x,Scalar y, Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_q_rho(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho(Scalar x,Scalar y, Scalar z, Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_q_rho(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_e(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_q_e(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_e(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_q_e(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_u(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_u(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_u(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_u(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_v(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_v(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_v(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_v(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_w(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_w(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_w(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_w(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_e(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_e(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_rho_e(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_q_rho_e(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_source_nu(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_q_nu(x,y,z);
}

  // --------------------------------
  // analytical terms
  // --------------------------------

template <typename Scalar>
Scalar MASA::masa_eval_exact_t(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_exact_t(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_t(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_exact_t(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_u(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_exact_u(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_u(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_exact_u(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_v(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_exact_v(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_v(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_exact_v(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_w(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_exact_w(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_w(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_exact_w(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_p(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_exact_p(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_p(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_exact_p(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho(Scalar x,Scalar y,Scalar z,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho(x,y,z,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_nu(Scalar x,Scalar y,Scalar t)
{
  return masa_master<Scalar>().get_ms().eval_exact_nu(x,y,t);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho_C(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho_C(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_exact_rho_C3(Scalar x,Scalar y,Scalar z)
{
  return masa_master<Scalar>().get_ms().eval_exact_rho_C3(x,y,z);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_t(Scalar x,Scalar y,Scalar z,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_t(x,y,z,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_t(Scalar x,Scalar y,Scalar z,Scalar t,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_t(x,y,z,t,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_u(Scalar x,Scalar y,Scalar z,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_u(x,y,z,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_u(Scalar x,Scalar y,Scalar z,Scalar t,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_u(x,y,z,t,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_v(Scalar x,Scalar y,Scalar z,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_v(x,y,z,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_v(Scalar x,Scalar y,Scalar z,Scalar t,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_v(x,y,z,t,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_w(Scalar x,Scalar y,Scalar z,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_w(x,y,z,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_w(Scalar x,Scalar y,Scalar z,Scalar t,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_w(x,y,z,t,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_p(Scalar x,Scalar y,Scalar z,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_p(x,y,z,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_p(Scalar x,Scalar y,Scalar z,Scalar t,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_p(x,y,z,t,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_rho(Scalar x,Scalar y,Scalar z,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_rho(x,y,z,i);
}

template <typename Scalar>
Scalar MASA::masa_eval_grad_rho(Scalar x,Scalar y,Scalar z,Scalar t,int i)
{
  return masa_master<Scalar>().get_ms().eval_g_rho(x,y,z,t,i);
}

/* ------------------------------------------------
 *
 *         utility functions
 *
 * -----------------------------------------------
 */ 


template <typename Scalar>
int MASA::masa_get_name(std::string* name)
{
  masa_master<Scalar>().get_ms().return_name(name); // set string to name
  return 0;
}


template <typename Scalar>
int MASA::masa_get_dimension(int* dim)
{
  masa_master<Scalar>().get_ms().return_dim(dim); // set string to name
  return 0;
}

template <typename Scalar>
int MASA::masa_test_poly()
{
  return masa_master<Scalar>().get_ms().poly_test(); // return error condition
}

template <typename Scalar>
int MASA::masa_sanity_check()
{
  return masa_master<Scalar>().get_ms().sanity_check(); // set string to name
}

int MASA::masa_version_stdout()
{
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::cout << "MASA Library: Version = " << MASA_LIB_VERSION;
  std::cout << " (" << MASA::masa_get_numeric_version() << ")" << std::endl << std::endl;

  std::cout << MASA_LIB_RELEASE << std::endl << std::endl;

  std::cout << "Build Date   = " << MASA_BUILD_DATE     << std::endl;
  std::cout << "Build Host   = " << MASA_BUILD_HOST     << std::endl;
  std::cout << "Build User   = " << MASA_BUILD_USER     << std::endl;
  std::cout << "Build Arch   = " << MASA_BUILD_ARCH     << std::endl;
  std::cout << "Build Rev    = " << MASA_BUILD_VERSION  << std::endl << std::endl;

  std::cout << "C++ Config   = " << MASA_CXX << " "     << MASA_CXXFLAGS << std::endl;
  std::cout << "F90 Config   = " << MASA_FC MASA_FCFLAGS << std::endl << std::endl;
  std::cout << "Optional Features:" << std::endl;
  std::cout << "   Python support enabled = ";
#ifdef SWIG_INTERFACES
  std::cout << "yes" << std::endl;
#else
  std::cout << "no" << std::endl;
#endif
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  return 0;
}

int MASA::masa_get_numeric_version()
{
  // Note: return format follows the versioning convention xx.yy.zz where
  //
  // xx = major version number
  // yy = minor version number
  // zz = micro version number
  //
  // For example:
  // v.   0.23  -> 002300 = 2300
  // v   0.23.1 -> 002301 = 2301
  // v. 10.23.2 -> 102302

  int major_version = 0;
  int minor_version = 0;
  int micro_version = 0;

  #ifdef MASA_MAJOR_VERSION
  major_version = MASA_MAJOR_VERSION;
  #endif

  #ifdef MASA_MINOR_VERSION
  minor_version = MASA_MINOR_VERSION;
  #endif

  #ifdef MASA_MICRO_VERSION
  micro_version = MASA_MICRO_VERSION;
  #endif

  return(major_version*10000 + minor_version*100 + micro_version);

}


// Instantiations

#define INSTANTIATE_ALL_FUNCTIONS(Scalar) \
  template int masa_init      <Scalar>(std::string, std::string); \
  template int masa_test_default <Scalar>(Scalar);		  \
  template int masa_select_mms<Scalar>(std::string); \
  template int masa_list_mms  <Scalar>(); \
  template int masa_purge_default_param <Scalar>(); \
  template Scalar pass_func             <Scalar>(Scalar (*)(Scalar),Scalar); \
  template int    masa_init_param<Scalar>(); \
  template void   masa_set_param<Scalar>(std::string,Scalar); \
  template Scalar masa_get_param<Scalar>(std::string); \
  template void   masa_set_vec<Scalar>(std::string,std::vector<Scalar>&); \
  template int masa_get_vec<Scalar>(std::string,std::vector<Scalar>&); \
  template Scalar masa_eval_source_t  <Scalar>(Scalar);         \
  template Scalar masa_eval_source_t  <Scalar>(Scalar,Scalar);  \
  template Scalar masa_eval_source_f  <Scalar>(Scalar,Scalar);  \
  template Scalar masa_eval_source_u  <Scalar>(Scalar); \
  template Scalar masa_eval_source_v  <Scalar>(Scalar);         \
  template Scalar masa_eval_source_w  <Scalar>(Scalar); \
  template Scalar masa_eval_source_e  <Scalar>(Scalar); \
  template Scalar masa_eval_source_e<Scalar>(Scalar,Scalar (*)(Scalar));  \
  template Scalar masa_eval_source_rho<Scalar>(Scalar); \
  template Scalar masa_eval_source_rho_u<Scalar>(Scalar);  \
  template Scalar masa_eval_source_rho_v<Scalar>(Scalar);  \
  template Scalar masa_eval_source_rho_w<Scalar>(Scalar);  \
  template Scalar masa_eval_source_rho_e<Scalar>(Scalar);  \
  template Scalar masa_eval_source_rho_N<Scalar>(Scalar,Scalar (*)(Scalar));  \
  template Scalar masa_eval_source_rho_N2<Scalar>(Scalar,Scalar (*)(Scalar));       \
  template Scalar masa_eval_source_boundary<Scalar>(Scalar);    \
  template Scalar masa_eval_source_C<Scalar>(Scalar);  \
  template Scalar masa_eval_source_C3<Scalar>(Scalar);      \
  template Scalar masa_eval_source_rho_C<Scalar>(Scalar);  \
  template Scalar masa_eval_source_rho_C3<Scalar>(Scalar);      \
  template Scalar masa_eval_exact_t      <Scalar>(Scalar);         \
  template Scalar masa_eval_exact_t      <Scalar>(Scalar,Scalar);  \
  template Scalar masa_eval_exact_u      <Scalar>(Scalar); \
  template Scalar masa_eval_exact_v      <Scalar>(Scalar);         \
  template Scalar masa_eval_exact_w      <Scalar>(Scalar); \
  template Scalar masa_eval_exact_p      <Scalar>(Scalar); \
  template Scalar masa_eval_exact_rho    <Scalar>(Scalar); \
  template Scalar masa_eval_exact_rho_N   <Scalar>(Scalar); \
  template Scalar masa_eval_exact_rho_N2  <Scalar>(Scalar); \
  template Scalar masa_eval_exact_rho_C   <Scalar>(Scalar); \
  template Scalar masa_eval_exact_rho_C3  <Scalar>(Scalar); \
  template Scalar masa_eval_posterior  <Scalar>(Scalar); \
  template Scalar masa_eval_prior      <Scalar>(Scalar); \
  template Scalar masa_eval_central_moment <Scalar>(int); \
  template Scalar masa_eval_posterior_mean <Scalar>(); \
  template Scalar masa_eval_posterior_variance <Scalar>(); \
  template Scalar masa_eval_likelyhood <Scalar>(Scalar); \
  template Scalar masa_eval_loglikelyhood <Scalar>(Scalar); \
  template Scalar masa_eval_exact_rho_C   <Scalar>(Scalar,Scalar);  \
  template Scalar masa_eval_exact_rho_C3  <Scalar>(Scalar,Scalar);      \
  template Scalar masa_eval_exact_rho_C   <Scalar>(Scalar,Scalar,Scalar);   \
  template Scalar masa_eval_exact_rho_C3  <Scalar>(Scalar,Scalar,Scalar);   \
  template Scalar masa_eval_source_t  <Scalar>(Scalar,Scalar,Scalar);  \
  template Scalar masa_eval_source_u  <Scalar>(Scalar,Scalar); \
  template Scalar masa_eval_source_v  <Scalar>(Scalar,Scalar); \
  template Scalar masa_eval_source_w  <Scalar>(Scalar,Scalar);  \
  template Scalar masa_eval_source_e  <Scalar>(Scalar,Scalar); \
  template Scalar masa_eval_source_rho_u<Scalar>(Scalar,Scalar);  \
  template Scalar masa_eval_source_rho_v<Scalar>(Scalar,Scalar);   \
  template Scalar masa_eval_source_rho_w<Scalar>(Scalar,Scalar);  \
  template Scalar masa_eval_source_rho_e<Scalar>(Scalar,Scalar);  \
  template Scalar masa_eval_source_rho<Scalar>(Scalar,Scalar); \
  template Scalar masa_eval_source_nu <Scalar>(Scalar,Scalar); \
  template Scalar masa_eval_exact_t      <Scalar>(Scalar,Scalar,Scalar);  \
  template Scalar masa_eval_exact_u      <Scalar>(Scalar,Scalar); \
  template Scalar masa_eval_exact_v      <Scalar>(Scalar,Scalar); \
  template Scalar masa_eval_exact_w      <Scalar>(Scalar,Scalar);  \
  template Scalar masa_eval_exact_p      <Scalar>(Scalar,Scalar); \
  template Scalar masa_eval_exact_rho    <Scalar>(Scalar,Scalar); \
  template Scalar masa_eval_exact_phi    <Scalar>(Scalar,Scalar); \
  template Scalar masa_eval_exact_nu     <Scalar>(Scalar,Scalar); \
  template Scalar masa_eval_source_t  <Scalar>(Scalar,Scalar,Scalar,Scalar);  \
  template Scalar masa_eval_source_u  <Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_source_v  <Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_source_w  <Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_source_e  <Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_source_rho_u<Scalar>(Scalar,Scalar,Scalar);     \
  template Scalar masa_eval_source_rho_v<Scalar>(Scalar,Scalar,Scalar);     \
  template Scalar masa_eval_source_rho_w<Scalar>(Scalar,Scalar,Scalar);     \
  template Scalar masa_eval_source_rho_e<Scalar>(Scalar,Scalar,Scalar);     \
  template Scalar masa_eval_source_rho<Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_source_nu <Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_exact_u      <Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_exact_v      <Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_exact_w      <Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_exact_p      <Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_exact_rho    <Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_exact_nu     <Scalar>(Scalar,Scalar,Scalar); \
  template Scalar masa_eval_exact_t      <Scalar>(Scalar,Scalar,Scalar,Scalar);  \
  template Scalar masa_eval_exact_u      <Scalar>(Scalar,Scalar,Scalar,Scalar); \
  template Scalar masa_eval_exact_v      <Scalar>(Scalar,Scalar,Scalar,Scalar); \
  template Scalar masa_eval_exact_w      <Scalar>(Scalar,Scalar,Scalar,Scalar);  \
  template Scalar masa_eval_exact_p      <Scalar>(Scalar,Scalar,Scalar,Scalar); \
  template Scalar masa_eval_exact_rho    <Scalar>(Scalar,Scalar,Scalar,Scalar); \
  template Scalar masa_eval_source_u  <Scalar>(Scalar,Scalar,Scalar,Scalar); \
  template Scalar masa_eval_source_v  <Scalar>(Scalar,Scalar,Scalar,Scalar); \
  template Scalar masa_eval_source_w  <Scalar>(Scalar,Scalar,Scalar,Scalar); \
  template Scalar masa_eval_source_e  <Scalar>(Scalar,Scalar,Scalar,Scalar); \
  template Scalar masa_eval_source_rho_u<Scalar>(Scalar,Scalar,Scalar,Scalar);     \
  template Scalar masa_eval_source_rho_v<Scalar>(Scalar,Scalar,Scalar,Scalar);     \
  template Scalar masa_eval_source_rho_w<Scalar>(Scalar,Scalar,Scalar,Scalar);     \
  template Scalar masa_eval_source_rho_e<Scalar>(Scalar,Scalar,Scalar,Scalar);     \
  template Scalar masa_eval_source_rho<Scalar>(Scalar,Scalar,Scalar,Scalar); \
  template Scalar masa_eval_grad_t  <Scalar>(Scalar); \
  template Scalar masa_eval_grad_t  <Scalar>(Scalar,Scalar,int); \
  template Scalar masa_eval_grad_t  <Scalar>(Scalar,Scalar,Scalar,int); \
  template Scalar masa_eval_grad_t  <Scalar>(Scalar,Scalar,Scalar,Scalar,int); \
  template Scalar masa_eval_grad_u  <Scalar>(Scalar); \
  template Scalar masa_eval_grad_u  <Scalar>(Scalar,Scalar,int); \
  template Scalar masa_eval_grad_u  <Scalar>(Scalar,Scalar,Scalar,int); \
  template Scalar masa_eval_grad_u  <Scalar>(Scalar,Scalar,Scalar,Scalar,int); \
  template Scalar masa_eval_grad_v  <Scalar>(Scalar); \
  template Scalar masa_eval_grad_v  <Scalar>(Scalar,Scalar,int); \
  template Scalar masa_eval_grad_v  <Scalar>(Scalar,Scalar,Scalar,int); \
  template Scalar masa_eval_grad_v  <Scalar>(Scalar,Scalar,Scalar,Scalar,int); \
  template Scalar masa_eval_grad_w  <Scalar>(Scalar); \
  template Scalar masa_eval_grad_w  <Scalar>(Scalar,Scalar,int); \
  template Scalar masa_eval_grad_w  <Scalar>(Scalar,Scalar,Scalar,int); \
  template Scalar masa_eval_grad_w  <Scalar>(Scalar,Scalar,Scalar,Scalar,int); \
  template Scalar masa_eval_grad_p  <Scalar>(Scalar); \
  template Scalar masa_eval_grad_p  <Scalar>(Scalar,Scalar,int); \
  template Scalar masa_eval_grad_p  <Scalar>(Scalar,Scalar,Scalar,int); \
  template Scalar masa_eval_grad_p  <Scalar>(Scalar,Scalar,Scalar,Scalar,int); \
  template Scalar masa_eval_grad_rho<Scalar>(Scalar); \
  template Scalar masa_eval_grad_rho<Scalar>(Scalar,Scalar,int); \
  template Scalar masa_eval_grad_rho<Scalar>(Scalar,Scalar,Scalar,int); \
  template Scalar masa_eval_grad_rho<Scalar>(Scalar,Scalar,Scalar,Scalar,int); \
  template int masa_test_poly<Scalar>();                             \
  template int masa_printid<Scalar>(); \
  template int masa_display_param<Scalar>(); \
  template int masa_display_vec<Scalar>(); \
  template int masa_get_name<Scalar>(std::string*); \
  template int masa_get_dimension<Scalar>(int*); \
  template int masa_sanity_check<Scalar>()

namespace MASA {

INSTANTIATE_ALL_FUNCTIONS(double);
INSTANTIATE_ALL_FUNCTIONS(long double);

}
