// $License$
// $Author$
// $Copyright$
// $Id:$

#include <cmasa.h>
#include <masa.h>
#include <string>

using namespace std;
using namespace MASA;

extern "C" int cmasa_init(const char* specificname,const char* functionname)
{
  string sn(specificname);
  string fn(functionname);

  masa_init(sn,fn);

  return 0;
}

extern "C" int cmasa_select_mms(const char* function_user_wants)
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

extern "C" int cmasa_list_mms()
{
  masa_list_mms();
  return 0;
}

extern "C" int cmasa_init_param()
{
  masa_init_param();
  return 0;
}

extern "C" int cmasa_sanity_check()
{
  masa_sanity_check();
  return 0;
}

extern "C" int cmasa_display_param()
{
  masa_display_param();
  return 0;
}

extern "C" int cmasa_set_param(const char* param,double val)
{
  return(masa_set_param(param,val));
}

extern "C" int cmasa_get_param(const char* param, double* val)
{
  return(masa_get_param(param,val));
}

// --------------------------------
//
//    Source and Analytical Terms
// 
// --------------------------------

// --------------------------------
// source term(s) -- 1D
// --------------------------------

extern "C" int cmasa_eval_1d_t_source  (double x,double* sol){return(masa_eval_t_source(x,sol));};
extern "C" int cmasa_eval_1d_u_source  (double x,double* sol){return(masa_eval_u_source(x,sol));};
extern "C" int cmasa_eval_1d_e_source  (double x,double* sol){return(masa_eval_e_source(x,sol));};
extern "C" int cmasa_eval_1d_rho_source(double x,double* sol){return(masa_eval_rho_source(x,sol));};

extern "C" int cmasa_eval_1d_t_an      (double x,double* sol){return(masa_eval_t_an(x,sol));};
extern "C" int cmasa_eval_1d_u_an      (double x,double* sol){return(masa_eval_u_an(x,sol));};
extern "C" int cmasa_eval_1d_p_an      (double x,double* sol){return(masa_eval_p_an(x,sol));};
extern "C" int cmasa_eval_1d_rho_an    (double x,double* sol){return(masa_eval_rho_an(x,sol));};

// --------------------------------
// source term(s) -- 2D
// --------------------------------

extern "C" int cmasa_eval_2d_t_source  (double x,double y,double* sol){return(masa_eval_t_source(x,y,sol));};
extern "C" int cmasa_eval_2d_u_source  (double x,double y,double* sol){return(masa_eval_u_source(x,y,sol));};
extern "C" int cmasa_eval_2d_v_source  (double x,double y,double* sol){return(masa_eval_v_source(x,y,sol));};
extern "C" int cmasa_eval_2d_e_source  (double x,double y,double* sol){return(masa_eval_e_source(x,y,sol));};
extern "C" int cmasa_eval_2d_rho_source(double x,double y,double* sol){return(masa_eval_rho_source(x,y,sol));};

extern "C" int cmasa_eval_2d_t_an      (double x,double y,double* sol){return(masa_eval_t_an(x,y,sol));};
extern "C" int cmasa_eval_2d_u_an      (double x,double y,double* sol){return(masa_eval_u_an(x,y,sol));};
extern "C" int cmasa_eval_2d_v_an      (double x,double y,double* sol){return(masa_eval_v_an(x,y,sol));};
extern "C" int cmasa_eval_2d_p_an      (double x,double y,double* sol){return(masa_eval_p_an(x,y,sol));};
extern "C" int cmasa_eval_2d_rho_an    (double x,double y,double* sol){return(masa_eval_rho_an(x,y,sol));};

// --------------------------------
// source term(s) -- 3D
// --------------------------------

extern "C" int cmasa_eval_3d_t_source  (double x,double y,double z,double* sol){return(masa_eval_t_source  (x,y,z,sol));};
extern "C" int cmasa_eval_3d_u_source  (double x,double y,double z,double* sol){return(masa_eval_u_source  (x,y,z,sol));};
extern "C" int cmasa_eval_3d_v_source  (double x,double y,double z,double* sol){return(masa_eval_v_source  (x,y,z,sol));};
extern "C" int cmasa_eval_3d_w_source  (double x,double y,double z,double* sol){return(masa_eval_w_source  (x,y,z,sol));};
extern "C" int cmasa_eval_3d_e_source  (double x,double y,double z,double* sol){return(masa_eval_e_source  (x,y,z,sol));};
extern "C" int cmasa_eval_3d_rho_source(double x,double y,double z,double* sol){return(masa_eval_rho_source(x,y,z,sol));};

extern "C" int cmasa_eval_3d_t_an      (double x,double y,double z,double* sol){return(masa_eval_t_an  (x,y,z,sol));};
extern "C" int cmasa_eval_3d_u_an      (double x,double y,double z,double* sol){return(masa_eval_u_an  (x,y,z,sol));};
extern "C" int cmasa_eval_3d_v_an      (double x,double y,double z,double* sol){return(masa_eval_v_an  (x,y,z,sol));};
extern "C" int cmasa_eval_3d_w_an      (double x,double y,double z,double* sol){return(masa_eval_w_an  (x,y,z,sol));};
extern "C" int cmasa_eval_3d_p_an      (double x,double y,double z,double* sol){return(masa_eval_p_an  (x,y,z,sol));};
extern "C" int cmasa_eval_3d_rho_an    (double x,double y,double z,double* sol){return(masa_eval_rho_an(x,y,z,sol));};
