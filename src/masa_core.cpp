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
// this is the core set of routines -- the only functions that should be called by users
// all other files are internal
//

#include <masa_internal.h>

using namespace MASA;

map<string, manufactured_solution*> masa_master_list; // global map between unique name and manufactured class
manufactured_solution* masa_master_pointer;           // pointer to currently selected manufactured solution

//
//  this function selects an already initialized manufactured class
//
int MASA::masa_select_mms(string name)
{
  string nametemp;
  int selector;
  // lets run though the list to check the variable does exist
  map<string,manufactured_solution*>::iterator it;
  it=masa_master_list.find(name);
  if(it != masa_master_list.end()) // found a name
    { 
      // cout << "selected " << name << endl;
      masa_master_list[name]=masa_master_pointer; // set pointer to currently selected solution      
    }      
  else 
    {
      cout << "\nMASA ERROR: No such manufactured solution (" << name << ") has been initialized.\n";
      exit(1);
    } 

  return 0;
  
}// done with masa_select

// helper function (deprecated!)
int masa_v2o(void* obid, manufactured_solution** manfac)
{ 
  string name;
  manufactured_solution* temp;
  temp=static_cast<manufactured_solution*>(obid); // cast the void ptr back as an obj
  *manfac=temp;
  return 0;

}

// get list mms
//
// this function returns a vector of pointers to 
// manufactured solutions available
//
// adding function because it is called in several places

int get_list_mms(vector<manufactured_solution*>* anim)
{
  anim->push_back(new MASA_Test()); // test function
  
  // register solutions here
  anim->push_back(new heateq_1d_steady_const());
  anim->push_back(new heateq_2d_steady_const());
  anim->push_back(new heateq_3d_steady_const());

  anim->push_back(new heateq_1d_unsteady_const());
  anim->push_back(new heateq_2d_unsteady_const());
  anim->push_back(new heateq_3d_unsteady_const());

  anim->push_back(new heateq_1d_unsteady_var());
  anim->push_back(new heateq_2d_unsteady_var());
  anim->push_back(new heateq_3d_unsteady_var());

  anim->push_back(new heateq_1d_steady_var());
  anim->push_back(new heateq_2d_steady_var());
  anim->push_back(new heateq_3d_steady_var());

  anim->push_back(new euler_1d());
  anim->push_back(new euler_2d());
  anim->push_back(new euler_3d());

  anim->push_back(new navierstokes_2d_compressible());
  anim->push_back(new navierstokes_3d_compressible());

  return 0;

}

//
//  this function will initiate a masa manufactured class
//
int MASA::masa_init(string unique_name, string str)
{
  int flag=0;
  string name,temp;
  int error=1;
  vector<manufactured_solution*> anim;

  get_list_mms(&anim); //construct list 
  
  temp=str;
  masa_map(&temp); // remove lowercase, etc.

  for (vector<manufactured_solution*>::const_iterator it = anim.begin(); it != anim.end(); ++it) 
    {
      (*it)->return_name(&name); // get name
      error=temp.rfind(name);   // look for name -- must be identical to full name, after masa_map edits
      if (error!=string::npos) // found a value
	{
	  masa_master_list[unique_name]=*it; // write down name 
	  masa_master_pointer=*it; // set as currently active manufactured solution
	  flag=1;                 // set exit flag
	}
      else // strings not identical
	{
	  delete *it; // this calls the deconstructor
	}
    }// done with for loop 
  
  
  if(flag==1)
    {      
      //cout << "\nMASA got it\n";
    }
  else
    {
      cout << "\nMASA FATAL ERROR: No Manufactured Solutions of that Type\n";
      exit(1); // error code, terminate
    }
  
  return 0; // steady as she goes

}

//
//  this function will initiate a masa manufactured class
//
int MASA::masa_curr_mms(string* str)
{

  // lets run though the list to check the variable does exist
  masa_master_pointer->return_name(str);
  // cout << masa_master_list[masa_master_pointer] << endl;    

  return 0;
}

//
//  this function will initiate a masa manufactured class
//
int MASA::masa_list_mms()
{
  string str;

  // output the size of the map
  cout << "Number of initialized solutions: " << masa_master_list.size() << endl;
  for(map<string,manufactured_solution*>::iterator iter = masa_master_list.begin(); iter != masa_master_list.end(); iter++)
    {
      (iter->second)->return_name(&str);
      cout << iter->first << " : " << str << endl;
    }
  return 0;
}

//
// function that searches all registered masa solutions
// for a selected manufactured solution
// (deprecated)
int MASA::masa_getid(void** objid,string str)
{
  int flag=0;
  string name;
  vector<manufactured_solution*> anim;

  get_list_mms(&anim); //construct list 
 
  // masa_map();  

  for (vector<manufactured_solution*>::const_iterator it = anim.begin(); it != anim.end(); ++it) {
    (*it)->return_name(&name); // get name
    if(str.compare(name) == 0) // strings are identical
      {
	*objid=(*it); // cast object pointer as void
	flag=1;      // set flag
      }
    else // strings not identical
      delete *it; // this calls the deconstructor
    
  }// done with for loop 

  if(flag==1)
    {      
      //cout << "\nMASA got it\n";
    }
  else
    {
      cout << "\nMASA FATAL ERROR: No Manufactured Solutions of that Type\n";
      exit(1); // error code, terminate
    }

  return 0; // steady as she goes
}// done with masa getid

//
// function that prints all registered masa solutions
//
int MASA::masa_printid()
{
  vector<manufactured_solution*> anim;
  string name;

  get_list_mms(&anim); //construct list 

  cout << endl;
  // masa_map();
  for (vector<manufactured_solution*>::const_iterator it = anim.begin(); it != anim.end(); ++it) {
    (*it)->return_name(&name); // get name
    cout << endl << name;
    delete *it; // this calls the deconstructor
  }// done with for loop 
  
  return 0; // steady as she goes
}// done with masa print id

int MASA::masa_set_param(string param,double paramval)
{
  masa_master_pointer->set_var(param,paramval);
  return 0;
}

//
// Set all parameters to default values
//
int MASA::masa_init_param()
{
  masa_master_pointer->init_var();
  return 0;
}

//
// Function that returns value of parameter selected by string
// 

int MASA::masa_get_param(string param,double *paramval)
{

  masa_master_pointer->get_var(param,paramval);

  return 0;
}

int MASA::masa_display_param()
{

  masa_master_pointer->display_var();

  return 0;
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

int MASA::masa_eval_t_source(double x,double* field) //x 
{

  *field=masa_master_pointer->eval_q_t(x);

  return 0;
}

int MASA::masa_eval_t_source(double x,double t,double* field) //x,t
{

  *field=masa_master_pointer->eval_q_t(x,t);

  return 0;
}

int MASA::masa_eval_u_source(double x,double* field)
{
  *field=masa_master_pointer->eval_q_u(x);
  return 0;
}

int MASA::masa_eval_rho_source(double x,double* field)
{
  *field=masa_master_pointer->eval_q_rho(x);
  return 0;
}

int MASA::masa_eval_e_source(double x,double* field)
{
  *field=masa_master_pointer->eval_q_e(x);
  return 0;
}

  // --------------------------------
  // analytical terms
  // --------------------------------

int MASA::masa_eval_t_an(double x,double* field)
{
  *field=masa_master_pointer->eval_an_t(x);
  return 0;
}

int MASA::masa_eval_u_an(double x,double* field)
{
  *field=masa_master_pointer->eval_an_u(x);
  return 0;
}

int MASA::masa_eval_p_an(double x,double* field)
{
  *field=masa_master_pointer->eval_an_p(x);
  return 0;
}

int MASA::masa_eval_rho_an(double x,double* field)
{
  *field=masa_master_pointer->eval_an_rho(x);
  return 0;
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

int MASA::masa_eval_t_source(double x,double y,double t,double* field)
{
  *field=masa_master_pointer->eval_q_t(x,y,t);
  return 0;
}

int MASA::masa_eval_u_source(double x,double y,double* field)
{

  *field=masa_master_pointer->eval_q_u(x,y);
  return 0;
}

int MASA::masa_eval_v_source(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_q_v(x,y);
  return 0;
}

int MASA::masa_eval_rho_source(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_q_rho(x,y);
  return 0;
}

int MASA::masa_eval_e_source(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_q_e(x,y);
  return 0;
}

  // --------------------------------
  // analytical terms
  // --------------------------------

int MASA::masa_eval_t_an(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_an_t(x,y);
  return 0;
}

int MASA::masa_eval_u_an(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_an_u(x,y);
  return 0;
}

int MASA::masa_eval_v_an(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_an_v(x,y);
  return 0;
}

int MASA::masa_eval_p_an(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_an_p(x,y);
  return 0;
}

int MASA::masa_eval_rho_an(double x,double y,double* field)
{
  *field=masa_master_pointer->eval_an_rho(x,y);
  return 0;
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

int MASA::masa_eval_t_source(double x,double y,double z,double t,double* field)
{
  *field=masa_master_pointer->eval_q_t(x,y,z,t);
  return 0;
}

int MASA::masa_eval_u_source(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_q_u(x,y,z);
  return 0;
}

int MASA::masa_eval_v_source(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_q_v(x,y,z);
  return 0;
}

int MASA::masa_eval_w_source(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_q_w(x,y,z);
  return 0;
}

int MASA::masa_eval_rho_source(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_q_rho(x,y,z);
  return 0;
}

int MASA::masa_eval_e_source(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_q_e(x,y,z);
  return 0;
}

  // --------------------------------
  // analytical terms
  // --------------------------------

int MASA::masa_eval_t_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_t(x,y,z);
  return 0;
}

int MASA::masa_eval_u_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_u(x,y,z);
  return 0;
}

int MASA::masa_eval_v_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_v(x,y,z);
  return 0;
}

int MASA::masa_eval_w_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_w(x,y,z);
  return 0;
}

int MASA::masa_eval_p_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_p(x,y,z);
  return 0;
}

int MASA::masa_eval_rho_an(double x,double y,double z,double* field)
{
  *field=masa_master_pointer->eval_an_rho(x,y,z);
  return 0;
}

/* ------------------------------------------------
 *
 *         utility functions
 *
 * -----------------------------------------------
 */ 

int MASA::masa_get_name(string* name)
{
  masa_master_pointer->return_name(name); // set string to name
  return 0;
}

int MASA::masa_get_dimension(int* dim)
{
  masa_master_pointer->return_dim(dim); // set string to name
  return 0;
}

int MASA::masa_sanity_check()
{
  masa_master_pointer->sanity_check(); // set string to name
  return 0;
}
