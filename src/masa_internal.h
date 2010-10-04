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
// $Author: nick $
// $Id: masa_core.cpp 12639 2010-08-24 23:33:29Z nick $
//
//
// masa_internal.h: internal MASA class definitions
//                  users are NOT recommended to muck around!
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa.h>
#include <math.h>
#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>

using namespace std;

namespace MASA 
{

  // masa map functions here
  // probably want to hide this from the user eventually
  
  int masa_map_solution  (string, string);
  int masa_map_temporal  (string, string);
  int masa_map_coeff     (string, string); 
  int masa_map           (string, string);
  int masa_map_dimension (string, string);
 
  // masa_shell
  //void masa_shell_choose_solution();

  /* -------------------------------------------------------------------------------------------   
   * manufactured_solution Base Class
   *
   * this is an abstract class
   * DO NOT EDIT THIS CLASS
   *
   * All manufactured solutions will inherit this class and it's methods 
   *
   * -------------------------------------------------------------------------------------------   
   */ 
  class manufactured_solution 
  {
    
  private: 
    double Tan;                          // analytical solution 
    double Q;                            // source term
    double gradT;                        // gradient 
    
  protected:
    static const double PI;              // 3.1415... defined in constructor
    static const double MASA_VAR_DEFAULT;   
    double dummy;
    int num_vars;

    map<string,int> varmap;               // map to each variable
    vector<double*>  vararr;              // arr of pointers to each variable
    string mmsname;                       // the name of the manufactured solution  
    int dimension;                        // dimension of the solution

  public: 

    // functions to override
    virtual ~manufactured_solution(){};       // destructor
    //virtual int init_var(){cout << "MASA ERROR:: NO DEFAULT VALUES AVAILABLE\n"; return 1;};                                                                  // inits all variables to selected values
    virtual int init_var() = 0;                                                                  // inits all variables to selected values

    // analytical solution(s)
    virtual double eval_an_t(double)                {cout << "MASA ERROR:: Analytical Solution (T) is unavailable or not properly loaded.\n"; return -1.33;}; // returns value of analytical solution
    virtual double eval_an_t(double,double)         {cout << "MASA ERROR:: Analytical Solution (T) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 2d problems
    virtual double eval_an_t(double,double,double)  {cout << "MASA ERROR:: Analytical Solution (T) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems
    virtual double eval_an_t(double,double,double,double)  {cout << "MASA ERROR:: Analytical Solution (T) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 4d problems

    virtual double eval_an_u(double)                {cout << "MASA ERROR:: Analytical Solution (u) is unavailable or not properly loaded.\n"; return -1.33;}; // returns value of analytical solution
    virtual double eval_an_u(double,double)         {cout << "MASA ERROR:: Analytical Solution (u) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 2d problems
    virtual double eval_an_u(double,double,double)  {cout << "MASA ERROR:: Analytical Solution (u) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems

    virtual double eval_an_v(double)                {cout << "MASA ERROR:: Analytical Solution (v) is unavailable for 1D problems.\n"; return -1.33;};        // returns value of analytical solution
    virtual double eval_an_v(double,double)         {cout << "MASA ERROR:: Analytical Solution (v) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 2d problems
    virtual double eval_an_v(double,double,double)  {cout << "MASA ERROR:: Analytical Solution (v) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems

    virtual double eval_an_w(double)                {cout << "MASA ERROR:: Analytical Solution (w) is unavailable for 1d problems.\n"; return -1.33;};        // returns value of analytical solution
    virtual double eval_an_w(double,double)         {cout << "MASA ERROR:: Analytical Solution (w) is unavailable for 2d problems.\n"; return -1.33;};        // overloaded for 2d problems
    virtual double eval_an_w(double,double,double)  {cout << "MASA ERROR:: Analytical Solution (w) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems
    
    virtual double eval_an_p(double)                {cout << "MASA ERROR:: Analytical Solution (e) is unavailable or not properly loaded.\n"; return -1.33;}; // returns value of analytical solution
    virtual double eval_an_p(double,double)         {cout << "MASA ERROR:: Analytical Solution (e) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 2d problems
    virtual double eval_an_p(double,double,double)  {cout << "MASA ERROR:: Analytical Solution (e) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems

    virtual double eval_an_rho(double)              {cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;}; // returns value of analytical solution
    virtual double eval_an_rho(double,double)       {cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 2d problems
    virtual double eval_an_rho(double,double,double){cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems
   
    // source terms
    virtual double eval_q_t(double)                {cout << "MASA ERROR:: Solution has not been properly loaded.\n"; return -1.33;};                    // returns value of source term (temp)
    virtual double eval_q_t(double, double)        {cout << "MASA ERROR:: Solution has not been properly loaded.\n"; return -1.33;};                    // returns value of source term (temp)
    virtual double eval_q_t(double,double,double)  {cout << "MASA ERROR:: Solution has not been properly loaded.\n"; return -1.33;};                    // returns value of source term (temp)
    virtual double eval_q_t(double,double,double,double)  {cout << "MASA ERROR:: Solution has not been properly loaded.\n"; return -1.33;};             // returns value of source term (x,y,z,t)

    virtual double eval_q_u(double)                {cout << "MASA ERROR:: Solution has not been properly loaded.\n"; return -1.33;};                    // returns value of source term (u)
    virtual double eval_q_u(double,double)         {cout << "MASA ERROR:: Solution has not been properly loaded.\n"; return -1.33;};                    // overloaded for 2d problems
    virtual double eval_q_u(double,double,double)  {cout << "MASA ERROR:: Solution has not been properly loaded.\n"; return -1.33;};                    // overloaded for 3d problems

    virtual double eval_q_v(double)                {cout << "MASA ERROR:: Source term (v) is unavailable for 1d problems -- eval_q_v has too few arguments.\n"; return -1.33;}; // returns value of source term (v)
    virtual double eval_q_v(double,double)         {cout << "MASA ERROR:: Source term (v) is unavailable or not properly loaded.\n"; return -1.33;};                            // overloaded for 2d problems
    virtual double eval_q_v(double,double,double)  {cout << "MASA ERROR:: Source term (v) is unavailable or not properly loaded.\n"; return -1.33;};                            // overloaded for 3d problems 

    virtual double eval_q_w(double)                {cout << "MASA ERROR:: Source term (w) is unavailable for 1d problems -- eval_q_w has too few arguments.\n"; return -1.33;};  // returns value of source term (w)
    virtual double eval_q_w(double,double)         {cout << "MASA ERROR:: Source term (w) is unavailable for 2d problems -- eval_q_w has too few arguments.\n"; return -1.33;};  // overloaded for 2d problems
    virtual double eval_q_w(double,double,double)  {cout << "MASA ERROR:: Source Term (w) is unavailable or not properly loaded.\n"; return -1.33;};                             // overloaded for 3d problems

    virtual double eval_q_e(double)                {cout << "MASA ERROR:: Source Term (e) is unavailable or not properly loaded.\n"; return -1.33;};    // returns value of source term (energy)
    virtual double eval_q_e(double,double)         {cout << "MASA ERROR:: Source Term (e) is unavailable or not properly loaded.\n"; return -1.33;};    // returns value of source term (energy)
    virtual double eval_q_e(double,double,double)  {cout << "MASA ERROR:: Source Term (e) is unavailable or not properly loaded.\n"; return -1.33;};    // returns value of source term (energy)

    virtual double eval_q_rho(double)              {cout << "MASA ERROR:: Source Term (rho) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (density)
    virtual double eval_q_rho(double,double)       {cout << "MASA ERROR:: Source Term (rho) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (density)
    virtual double eval_q_rho_u(double,double)     {cout << "MASA ERROR:: Source Term (rho) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (density) -- 1d Sod
    virtual double eval_q_rho(double,double,double){cout << "MASA ERROR:: Source Term (rho) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (density)

    // gradient 
    virtual double eval_1d_g(int,double)               {cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};         // returns value of 1d gradient
    virtual double eval_2d_g(int,double,double)        {cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};         // returns value of 2d gradient
    virtual double eval_3d_g(int,double,double,double) {cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};         // returns value of 3d gradient

    // member functions solution classes will inherit
    manufactured_solution();                                     // constructor
    int set_var(string,double);                                  // sets variable value    
    int register_var(string, double*);                           // this registers a variable
    int sanity_check();                                          // checks that all variables to the class have been initalized
    int poly_test();                                             // regression method for poly class (see below)
    double get_var(string);                                      // returns variable value
    void display_var();                                          // print all variable names and values
    void return_name(string* inname){inname->assign(mmsname);};  // method: returns name
    void return_dim (int* indim)    {*indim=dimension;};         // method: returns dimension of solution
    
  }; // done with MMS base class

  /* -------------------------------------------------------------------------------------------   
   * Polynomial Base Class
   *
   * Blatantly stealing paul bauman's polynomial class definitions in the name of science
   * 
   * In addition to inheriting manufactured class, must inherit polynomial 
   *
   * -------------------------------------------------------------------------------------------   
   */ 
  class Polynomial
  {
  public:
    
    void set_coeffs( const vector<double> & );
    
    // Evaluates polynomial.
    double operator()( const double &, int *) const;
    
    // Evaluates polynomial and deriviatives up to order specified by user.
    void eval_derivs( const double, const int, vector<double> & ) const;
    
    double get_coeffs( const int & ) const;    
    
  protected:
    
    // We assume that the coefficents are ordered as follows:
    // y(x) = a0 + a1*x + a2*x^2 + ... + an*x^n
    // so that there should be n+1 coeff so that
    // coeffs[0] = a0
    // coeffs[1] = a1
    // and so on.
    vector<double> coeffs;
    
  };  

  // ------------------------------------------------------
  // ---------- all other mms classes ------------
  // ------------------------------------------------------

  // just a demo class
  class masa_test_function : public manufactured_solution 
  {
  private:
    double dummy;
    double demo_var_2;
    double demo_var_3;
  public:
    masa_test_function(); // constructor
    int init_var();        // default problem values
    int poly_test();
    //double eval_q_t(double);
  }; // done with masa_test


  // class with no source/analytical terms, to test virtual function defaults
  class masa_uninit : public manufactured_solution 
  {
  private:

  public:
    masa_uninit();  // constructor
    int init_var(); // default problem values
  }; 

  // ------------------------------------------------------
  // ---------- heat equation /steady / constant ------------
  // ------------------------------------------------------

  class heateq_1d_steady_const : public manufactured_solution 
  {
  private:
    double A_x;
    double k_0;    
  public:
    heateq_1d_steady_const(); // constructor
    int init_var();        // default problem values
    double eval_q_t (double);  // source term evaluator
    double eval_an_t(double);   //analytical solution

  };

  class heateq_2d_steady_const : public manufactured_solution 
  {
  private:
    double A_x;
    double k_0;
    double B_y;

  public:
    heateq_2d_steady_const();       // constructor
    int init_var();        // default problem values
    double eval_q_t (double,double); // source term evaluator
    double eval_an_t(double,double); // analytical term evaluator
  };
  
  class heateq_3d_steady_const : public manufactured_solution 
  {
  private:
    double A_x;
    double k_0;
    double B_y;
    double C_z;

  public:
    heateq_3d_steady_const(); // constructor
    int init_var();        // default problem values
    double eval_q_t (double,double,double); //evaluate source term
    double eval_an_t(double,double,double); // analytical term evaluator
  };
  // ------------------------------------------------------
  // ---------- heat equation / unsteady / constant -------
  // ------------------------------------------------------

  class heateq_1d_unsteady_const : public manufactured_solution 
  {
  private:
    double A_x;
    double A_t;
    double D_t;
    double k_0;
    double cp_0;
    double rho;
    
  public:
    heateq_1d_unsteady_const(); // constructor
    int init_var();        // default problem values
    double eval_q_t(double,double); // needs x,t
  };

  class heateq_2d_unsteady_const : public manufactured_solution 
  {
  private:
    double A_x;
    double A_t;
    double B_y;
    double B_t;
    double D_t;
    double rho;
    double k_0;
    double cp_0;

  public:
    heateq_2d_unsteady_const(); // constructor
    int init_var();        // default problem values
    double eval_q_t(double,double,double); // needs x,y,t
  };
  
  class heateq_3d_unsteady_const : public manufactured_solution 
  {
  private:
    double A_x;
    double A_t;
    double B_y;
    double B_t;
    double C_z;
    double C_t;
    double D_t;
    double k_0;
    double cp_0;
    double rho;
    
  public:
    heateq_3d_unsteady_const(); // constructor
    int init_var();        // default problem values
    double eval_q_t(double,double,double,double); // needs x,y,z,t    
  };

  // ------------------------------------------------------
  // ---------- heat equation / unsteady / var ------------
  // ------------------------------------------------------

  class heateq_1d_unsteady_var : public manufactured_solution 
  {
  private:
    double A_x;
    double A_t;
    double D_t;
    double rho;
    double k_0;
    double k_1;
    double k_2;
    double cp_0;
    double cp_1;
    double cp_2;
    
  public:
    heateq_1d_unsteady_var(); // constructor
    int init_var();        // default problem values
    double eval_q_t(double,double); // needs x,t
  };

  class heateq_2d_unsteady_var : public manufactured_solution 
  {
  private:
    double A_x;
    double A_t;
    double B_y;
    double B_t;
    double D_t;
    double rho;
    double k_0;
    double k_1;
    double k_2;
    double cp_0;
    double cp_1;
    double cp_2;

  public:
    heateq_2d_unsteady_var(); // constructor
    int init_var();        // default problem values
    double eval_q_t(double,double,double); // needs x,y,t    
  };
  
  class heateq_3d_unsteady_var : public manufactured_solution 
  {
  private:
    double A_x;
    double A_t;
    double B_y;
    double B_t;
    double C_z;
    double C_t;
    double D_t;
    double rho;
    double k_0;
    double k_1;
    double k_2;
    double cp_0;
    double cp_1;
    double cp_2;

  public:
    heateq_3d_unsteady_var(); // constructor
    int init_var();        // default problem values
    double eval_q_t(double,double,double,double); // needs x,y,z,t    
  };

  // ------------------------------------------------------
  // ---------- heat equation / steady / var ------------
  // ------------------------------------------------------

  class heateq_1d_steady_var : public manufactured_solution 
  {
  private:
    double A_x;
    double k_0;
    double k_1;
    double k_2;

  public:
    heateq_1d_steady_var(); // constructor
    int init_var();        // default problem values
    double eval_q_t(double); // needs x
  };

  class heateq_2d_steady_var : public manufactured_solution 
  {
  private:
    double A_x;
    double k_0;
    double k_1;
    double k_2;
    double B_y;

  public:
    heateq_2d_steady_var(); // constructor
    int init_var();        // default problem values
    double eval_q_t(double,double); // needs x,y
  };
  
  class heateq_3d_steady_var : public manufactured_solution 
  {
  private:
    double A_x;
    double k_0;
    double k_1;
    double k_2;
    double B_y;
    double C_z;

  public:
    heateq_3d_steady_var(); // constructor
    int init_var();        // default problem values
    double eval_q_t(double,double,double); // needs x,y,z
    
  };

  // ------------------------------------------------------
  // ---------- euler ------------
  // ------------------------------------------------------
  class euler_1d : public manufactured_solution
  {

    double R;                             // (or is this the ideal gas constant?) ratio of specific heat capacities, defined in constructor
    double k;                             // Boltzmanns constant

    double u_0;
    double u_x;
    double rho_0;
    double rho_x;
    double p_0;
    double p_x;
    double a_px;
    double a_rhox;
    double a_ux;
    double Gamma;
    double mu;
    double L;
    

  public:
    euler_1d(); // constructor    
    int init_var();          // default problem values

    double eval_q_u   (double); // source terms
    double eval_q_e   (double);
    double eval_q_rho (double);
    
    double eval_an_u  (double); // analytical
    double eval_an_p  (double);
    double eval_an_rho(double);

  };

  class euler_2d : public manufactured_solution
  {

    double R;     // (or is this the ideal gas constant?) ratio of specific heat capacities, defined in constructor
    double k;     // Boltzmanns constant

    double u_0;
    double u_x;
    double u_y;
    double v_0;
    double v_x;
    double v_y;
    double rho_0;
    double rho_x;
    double rho_y;
    double p_0;
    double p_x;
    double p_y;
    double a_px;
    double a_py;
    double a_rhox;
    double a_rhoy;
    double a_ux;
    double a_uy;
    double a_vx;
    double a_vy;
    double Gamma;
    double mu;
    double L;
    
  public:
    euler_2d(); // constructor    
    int init_var();        // default problem values

    double eval_q_u  (double,double);
    double eval_q_v  (double,double);
    double eval_q_e  (double,double);
    double eval_q_rho(double,double);

    double eval_an_u  (double,double); // analytical
    double eval_an_v  (double,double);
    double eval_an_p  (double,double);
    double eval_an_rho(double,double);

  };

  class euler_3d : public manufactured_solution
  {
    double R;     // (or is this the ideal gas constant?) ratio of specific heat capacities, defined in constructor
    double k;     // Boltzmanns constant

    double u_0;
    double u_x;
    double u_y;
    double u_z;
    double v_0;
    double v_x;
    double v_y;
    double v_z;
    double w_0;
    double w_x;
    double w_y;
    double w_z;
    double rho_0;
    double rho_x;
    double rho_y;
    double rho_z;
    double p_0;
    double p_x;
    double p_y;
    double p_z;
    double a_px;
    double a_py;
    double a_pz;
    double a_rhox;
    double a_rhoy;
    double a_rhoz;
    double a_ux;
    double a_uy;
    double a_uz;
    double a_vx;
    double a_vy;
    double a_vz;
    double a_wx;
    double a_wy;
    double a_wz;
    double mu;
    double Gamma;
    double L;    
    
  public:
    euler_3d(); // constructor
    int init_var();        // default problem values

    double eval_q_u  (double,double,double); // source terms
    double eval_q_v  (double,double,double);
    double eval_q_w  (double,double,double);
    double eval_q_e  (double,double,double);
    double eval_q_rho(double,double,double);

    double eval_an_u  (double,double,double); // analytical
    double eval_an_v  (double,double,double);
    double eval_an_w  (double,double,double);
    double eval_an_q  (double,double,double);
    double eval_an_p  (double,double,double);
    double eval_an_rho(double,double,double);

  };

  // ------------------------------------------------------
  // ---------- axisymmetric euler ------------
  // ------------------------------------------------------
  class axi_euler : public manufactured_solution
  {

    double R;                             // (or is this the ideal gas constant?) ratio of specific heat capacities, defined in constructor

    double p_0;
    double p_1;
    double rho_0;
    double rho_1;
    double u_1;
    double w_0;
    double w_1;
    double a_pr;
    double a_pz;
    double a_rhor;
    double a_rhoz;
    double a_ur;
    double a_uz;
    double a_wr;
    double a_wz;
    double L;
    double mu;
    double Gamma;    

  public:
    axi_euler(); // constructor    
    int init_var();          // default problem values

    double eval_q_u   (double,double); // radial velocity 
    double eval_q_w   (double,double); // axial 
    double eval_q_e   (double,double);
    double eval_q_rho (double,double);
    
    double eval_an_u  (double,double); // analytical
    double eval_an_w  (double,double);
    double eval_an_p  (double,double);
    double eval_an_rho(double,double);

  };

  // ------------------------------------------------------
  // ---------- axisymmetric compressible navier stokes ---
  // ------------------------------------------------------
  class axi_cns : public manufactured_solution
  {

    double R;                             // (or is this the ideal gas constant?) ratio of specific heat capacities, defined in constructor
    double k;                             // Boltzmanns constant

    double p_0;
    double p_1;
    double rho_0;
    double rho_1;
    double u_1;
    double w_0;
    double w_1;
    double a_pr;
    double a_pz;
    double a_rhor;
    double a_rhoz;
    double a_ur;
    double a_uz;
    double a_wr;
    double a_wz;
    double L;
    double mu;
    double Gamma;
   
  public:
    axi_cns(); // constructor    
    int init_var();          // default problem values

    double eval_q_u   (double,double); // radial velocity 
    double eval_q_w   (double,double); // axial 
    double eval_q_e   (double,double);
    double eval_q_rho (double,double);
    
    double eval_an_u  (double,double); // analytical
    double eval_an_w  (double,double);
    double eval_an_p  (double,double);
    double eval_an_rho(double,double);

  };

  // ------------------------------------------------------
  // ---------- sod 1d ------------
  // ------------------------------------------------------
  class sod_1d : public manufactured_solution
  {
    double Gamma;
    double mu;

  public:
    sod_1d(); // constructor    
    int init_var();          // default problem values

    double eval_q_rho   (double,double);
    double eval_q_rho_u (double,double);
    double func         (double);  // helper function
    double rtbis        (double,double,double,int);
    double eval_q_t     (double); // regression test function
    double eval_q_t     (double,double); // regression test function
  };

  // ------------------------------------------------------
  // ---------- RANS: Spelart Alamaras Channel  ------------
  // ------------------------------------------------------
  class rans_sa : public manufactured_solution
  {
    
    double cb1;
    double cb2;
    double cv1;
    double cw2;
    double cw3;
    double sigma;
    double kappa;  
    double re_tau;

    // adding for modified SA
    double cv2;
    double cv3;

    // parameters not for users
    double etam;
    double a1;
    double b1;

  public:
    rans_sa(); // constructor    
    int init_var();          // default problem values
    
    double eval_q_u (double); // velocity term
    double eval_q_v (double); // eddy viscosity term
    double eval_an_u(double); // analytical
    double eval_an_v(double);

    // member functions not exposed by API
    double   u(double);
    double  du(double);
    double d2u(double);

    double   nu(double);
    double  dnu(double);
    double d2nu(double);

    double production(double);
    double destruction(double);
    double transport(double);
    
    double cw1();
    double chi(double);
    double fv1(double);
    double fv2(double);
    double   s(double);
    double  sb(double);
    double   r(double);
    double   g(double);
    double  fw(double);
    double dvt(double);
    double  vt(double);

  };

  // ------------------------------------------------------
  // ---------- compressible navier stokes  ------------
  // ------------------------------------------------------

  class navierstokes_2d_compressible : public manufactured_solution
  {    
    double R;     // (or is this the ideal gas constant?) ratio of specific heat capacities, defined in constructor
    double k;     // Boltzmanns constant

    double u_0;
    double u_x;
    double u_y;
    double v_0;
    double v_x;
    double v_y;
    double rho_0;
    double rho_x;
    double rho_y;
    double p_0;
    double p_x;
    double p_y;
    double a_px;
    double a_py;
    double a_rhox;
    double a_rhoy;
    double a_ux;
    double a_uy;
    double a_vx;
    double a_vy;
    double Gamma;
    double mu;
    double L;
    
  public:
    navierstokes_2d_compressible(); // constructor
    int init_var();        // default problem values

    double eval_q_u  (double,double);
    double eval_q_v  (double,double);
    double eval_q_e  (double,double);
    double eval_q_rho(double,double);

    double eval_an_u  (double,double); // analytical
    double eval_an_v  (double,double);
    double eval_an_p  (double,double);
    double eval_an_rho(double,double);

  };
  
  class navierstokes_3d_compressible : public manufactured_solution
  {
    double R;     // (or is this the ideal gas constant?) ratio of specific heat capacities, defined in constructor
    double k;     // Boltzmanns constant

    double u_0;
    double u_x;
    double u_y;
    double u_z;
    double v_0;
    double v_x;
    double v_y;
    double v_z;
    double w_0;
    double w_x;
    double w_y;
    double w_z;
    double rho_0;
    double rho_x;
    double rho_y;
    double rho_z;
    double p_0;
    double p_x;
    double p_y;
    double p_z;
    double a_px;
    double a_py;
    double a_pz;
    double a_rhox;
    double a_rhoy;
    double a_rhoz;
    double a_ux;
    double a_uy;
    double a_uz;
    double a_vx;
    double a_vy;
    double a_vz;
    double a_wx;
    double a_wy;
    double a_wz;
    double mu;
    double Gamma;
    double L;    
    
  public:
    navierstokes_3d_compressible(); //constructor
    int init_var();        // default problem values

    double eval_q_u  (double,double,double); // source terms
    double eval_q_v  (double,double,double);
    double eval_q_w  (double,double,double);
    double eval_q_e  (double,double,double);
    double eval_q_rho(double,double,double);

    double eval_an_u  (double,double,double); // analytical
    double eval_an_v  (double,double,double);
    double eval_an_w  (double,double,double);
    double eval_an_q  (double,double,double);
    double eval_an_p  (double,double,double);
    double eval_an_rho(double,double,double);

  }; // done with navier stokes 3d class
  
} // end MASA namespace
