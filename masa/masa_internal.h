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
//   These are the internal MASA class definitions
//   Users are NOT recommended to muck around with this
//

#include <masa.h>
#include <math.h>

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
    double Tan;                             // analytical solution 
    double Q;                               // source term
    double gradT;                           // gradient 
    
  protected:
    double MASA_VAR_DEFAULT;
    map<string,int> varmap;               // map to each variable
    vector<double*>  vararr;              // arr of pointers to each variable
    string mmsname;                       // the name of the manufactured solution  
    int dimension;                        // dimension of the solution

  public: 

    // functions to override
    virtual ~manufactured_solution(){};       // destructor

    // analytical solution
    virtual double eval_an()         {cout << "MASA ERROR: Analytical Solution is unavailable or not properly loaded.";  return -1.33;}; // returns value of analytical solution

    // source terms
    virtual double eval_q_u(double)  {cout << "MASA ERROR: Solution has not been properly loaded."; return -1.33;};                      // returns value of source term (u)
    virtual double eval_q_v(double)  {cout << "MASA ERROR: Source term (v) is unavailable or not properly loaded."; return -1.33;};      // returns value of source term (v)
    virtual double eval_q_w(double)  {cout << "MASA ERROR: Source Term (w) is unavailable or not properly loaded."; return -1.33;};      // returns value of source term (w)
    virtual double eval_q_e(double)  {cout << "MASA ERROR: Source Term (e) is unavailable or not properly loaded."; return -1.33;};      // returns value of source term (energy)
    virtual double eval_q_rho(double){cout << "MASA ERROR: Source Term (rho) is unavailable or not properly loaded."; return -1.33;};    // returns value of source term (density)

    // gradient 
    virtual double eval_1d_g(double) {cout << "MASA ERROR: gradient is unavailable or not properly loaded.";   return -1.33;};           // returns value of gradient
    virtual double eval_2d_g(double) {cout << "MASA ERROR: gradient is unavailable or not properly loaded.";   return -1.33;};           // returns value of gradient
    virtual double eval_3d_g(double) {cout << "MASA ERROR: gradient is unavailable or not properly loaded.";   return -1.33;};           // returns value of gradient

    // member functions solution classes will inherit
    manufactured_solution();                                     // constructor
    void get_var(string,double*);                                // returns variable value
    void set_var(string,double);                                 // sets variable value    
    void display_var();                                          // print all variable names and values
    void sanity_check();                                         // checks that all variables to the class have been initalized
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
    double operator()( const double & ) const;
    
    // Evaluates polynomial and deriviatives up to order specified by user.
    void eval_derivs( const double &, const int &, vector<double> & ) const;
    
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
  class MASA_Test : public manufactured_solution 
  {
  private:
    double axp;
    double demo_var_2;
    double demo_var_3;
  public:
    MASA_Test(); // constructor
    double eval_q_u(double);
  }; // done with heat_eq_1d

  // ------------------------------------------------------
  // ---------- heat equation /steady / constant ------------
  // ------------------------------------------------------

  class heateq_1d_steady_const : public manufactured_solution 
  {
  private:
    double axp;

    double ax;
    double k0;    
  public:
    heateq_1d_steady_const(); // constructor
    double eval_q_u(double)      ; // source term evaluator
  };

  class heateq_2d_steady_const : public manufactured_solution 
  {
  private:
    double axp;

    double ax;
    double k0;
    double by;

  public:
    heateq_2d_steady_const(); // constructor
    double eval_q_u(double)      ; // source term evaluator
  };
  
  class heateq_3d_steady_const : public manufactured_solution 
  {
  private:
    double axp;

    double ax;
    double k0;
    double by;
    double cz;

  public:
    heateq_3d_steady_const(); // constructor
    
  };
  // ------------------------------------------------------
  // ---------- heat equation / unsteady / constant -------
  // ------------------------------------------------------

  class heateq_1d_unsteady_const : public manufactured_solution 
  {
  private:
    double axp;

    double ax;
    double k0;
    double dt;
    double cp0;
    double at;
    double rho;
    
  public:
    heateq_1d_unsteady_const(); // constructor
    
  };

  class heateq_2d_unsteady_const : public manufactured_solution 
  {
  private:
    double axp;

    double ax;
    double k0;
    double dt;
    double cp0;
    double at;
    double rho;
    double by;
    double bt;

  public:
    heateq_2d_unsteady_const(); // constructor
    
  };
  
  class heateq_3d_unsteady_const : public manufactured_solution 
  {
  private:
    double axp;

    double ax;
    double k0;
    double dt;
    double cp0;
    double at;
    double rho;
    double by;
    double bt;
    double cz;
    double ct;
    
  public:
    heateq_3d_unsteady_const(); // constructor
    
  };

  // ------------------------------------------------------
  // ---------- heat equation / unsteady / var ------------
  // ------------------------------------------------------

  class heateq_1d_unsteady_var : public manufactured_solution 
  {
  private:
    double axp;

    double ax;
    double k0;
    double dt;
    double cp0;
    double at;
    double rho;
    double cp1;
    double cp2;
    double k1;
    double k2;
    
  public:
    heateq_1d_unsteady_var(); // constructor
    
  };

  class heateq_2d_unsteady_var : public manufactured_solution 
  {
  private:
    double axp;

    double ax;
    double k0;
    double dt;
    double cp0;
    double at;
    double rho;
    double cp1;
    double cp2;
    double k1;
    double k2;
    double by;
    double bt;

  public:
    heateq_2d_unsteady_var(); // constructor
    
  };
  
  class heateq_3d_unsteady_var : public manufactured_solution 
  {
  private:
    double axp;

    double ax;
    double k0;
    double dt;
    double cp0;
    double at;
    double rho;
    double cp1;
    double cp2;
    double k1;
    double k2;
    double by;
    double bt;
    double cz;
    double ct;

  public:
    heateq_3d_unsteady_var(); // constructor
    
  };

  // ------------------------------------------------------
  // ---------- heat equation / steady / var ------------
  // ------------------------------------------------------

  class heateq_1d_steady_var : public manufactured_solution 
  {
  private:
    double axp;

    double ax;
    double k0;
    double k1;
    double k2;

  public:
    heateq_1d_steady_var(); // constructor
    
  };

  class heateq_2d_steady_var : public manufactured_solution 
  {
  private:
    double axp;
    
    double ax;
    double k0;
    double k1;
    double k2;
    double by;

  public:
    heateq_2d_steady_var(); // constructor
    
  };
  
  class heateq_3d_steady_var : public manufactured_solution 
  {
  private:
    double axp;

    double ax;
    double k0;
    double k1;
    double k2;
    double by;
    double cz;

  public:
    heateq_3d_steady_var(); // constructor
    
  };

  // ------------------------------------------------------
  // ---------- euler ------------
  // ------------------------------------------------------

  class euler_2d : public manufactured_solution
  {
    double axp;

    double u0;
    double ux;
    double uy;
    double v0;
    double vx;
    double vy;
    double rho0;
    double rhox;
    double rhoy;
    double p0;
    double px;
    double py;
    double apx;
    double apy;
    double arhox;
    double arhoy;
    double aux;
    double auy;
    double avx;
    double avy;
    double gamma;
    double mu;
    double l;
    
  public:
    euler_2d(); // constructor    

  };

  class euler_3d : public manufactured_solution
  {
    double axp;
    
    double u0;
    double ux;
    double uy;
    double uz;
    double v0;
    double vx;
    double vy;
    double vz;
    double w0;
    double wx;
    double wy;
    double wz;
    double rho0;
    double rhox;
    double rhoy;
    double rhoz;
    double p0;
    double px;
    double py;
    double pz;
    double apx;
    double apy;
    double apz;
    double arhox;
    double arhoy;
    double arhoz;
    double aux;
    double auy;
    double auz;
    double avx;
    double avy;
    double avz;
    double awx;
    double awy;
    double awz;
    double mu;
    double gamma;
    double l;
    
  public:
    euler_3d(); // constructor
    
  };

  // ------------------------------------------------------
  // ---------- compressible navier stokes  ------------
  // ------------------------------------------------------

  class navierstokes_compressible_2d : public manufactured_solution
  {    
    double axp;

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
    navierstokes_compressible_2d(); // constructor
    
  };
  
  class navierstokes_compressible_3d : public manufactured_solution
  {
    double axp;

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
    navierstokes_compressible_3d(); //constructor
    
  }; // done with navier stokes 3d class
  
} // end MASA namespace
