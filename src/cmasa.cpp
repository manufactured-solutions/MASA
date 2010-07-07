#include <cmasa.h>
#include <masa.h>
#include <string>

using namespace std;
using namespace MASA;

int cmasa_init(const char* specificname,const char* functionname)
{
  string sn(specificname);
  string fn(functionname);

  masa_init(sn,fn);

  return 0;
}

int cmasa_select_mms(const char* function_user_wants)
{
  string fuw(function_user_wants);
  masa_select_mms(fuw);
  return 0;
}

/*
int MASA::cmasa_curr_mms(char* function_name)
{
  string fn(function_name);
  masa_curr_mms(fn);
  return 0;
}
*/

int cmasa_list_mms()
{
  masa_list_mms();
  return 0;
}

int cmasa_init_param()
{
  masa_init_param();
  return 0;
}

int cmasa_sanity_check()
{
  masa_sanity_check();
  return 0;
}

int cmasa_display_param()
{
  masa_display_param();
  return 0;
}

// --------------------------------
//
//    Source and Analytical Terms
// 
// --------------------------------

// --------------------------------
// source term(s) -- 1D
// --------------------------------
int cmasa_eval_t_source  (double x ,double t,double* sol){masa_eval_t_source  (x ,t,sol);return 0;}
int cmasa_eval_t_source  (double x,double* sol){masa_eval_t_source  (x,sol); return 0;}
int cmasa_eval_u_source  (double x,double* sol){masa_eval_u_source  (x,sol); return 0;}
int cmasa_eval_e_source  (double x,double* sol){masa_eval_e_source  (x,sol); return 0;}
int cmasa_eval_rho_source(double x,double* sol){masa_eval_rho_source(x,sol); return 0;}

/*int cmasa_eval_t_an      (double,double*);        // x
int cmasa_eval_t_an      (double,double,double*); // x,t
int cmasa_eval_u_an      (double,double*);
int cmasa_eval_p_an      (double,double*);
int cmasa_eval_rho_an    (double,double*);

// --------------------------------
// source term(s) -- 2D
// --------------------------------
int cmasa_eval_t_source  (double,double,double,double*); //x,y,t
int cmasa_eval_u_source  (double,double,double*);
int cmasa_eval_v_source  (double,double,double*);
int cmasa_eval_e_source  (double,double,double*);
int cmasa_eval_rho_source(double,double,double*);

int cmasa_eval_t_an      (double,double,double,double*); //x,y,t
int cmasa_eval_u_an      (double,double,double*);
int cmasa_eval_v_an      (double,double,double*);
int cmasa_eval_p_an      (double,double,double*);
int cmasa_eval_rho_an    (double,double,double*);

// --------------------------------
// source term(s) -- 3D
// --------------------------------
int cmasa_eval_t_source  (double,double,double,double,double*); // x,y,z,t
int cmasa_eval_u_source  (double,double,double,double*);
int cmasa_eval_v_source  (double,double,double,double*);
int cmasa_eval_w_source  (double,double,double,double*);
int cmasa_eval_e_source  (double,double,double,double*);
int cmasa_eval_rho_source(double,double,double,double*);

int cmasa_eval_t_an      (double,double,double,double,double*); // x,y,z,t
int cmasa_eval_u_an      (double,double,double,double*);
int cmasa_eval_v_an      (double,double,double,double*);
int cmasa_eval_w_an      (double,double,double,double*);
int cmasa_eval_p_an      (double,double,double,double*);
int cmasa_eval_rho_an    (double,double,double,double*);
*/
