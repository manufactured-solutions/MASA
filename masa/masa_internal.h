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
  
  // this is an abstract class
  // DO NOT EDIT THIS CLASS
  class manufactured_solution 
  {
    
  private: 
    double Tan;                             // analytical solution 
    double Q;                               // source term
    double gradT;                           // gradient 
    
  protected:
    string mmsname;                       // the name of the manufactured solution  
    map<string,int> varmap;              // map to each variable
    vector<double*>  vararr;              // arr of pointers to each variable

  public: 

    // functions to override
    virtual ~manufactured_solution(){};       // destructor
    virtual double eval_an(){cout << "MASA ERROR";  return -1.33;}; // returns value of solution
    virtual double eval_q_u(double){cout << "MASA ERROR"; return -1.33;}; // returns value of source term (u)
    virtual double eval_q_v(double){cout << "MASA ERROR"; return -1.33;}; // returns value of source term (v)
    virtual double eval_q_w(double){cout << "MASA ERROR"; return -1.33;}; // returns value of source term (w)
    virtual double eval_g(double){cout << "MASA ERROR";   return -1.33;}; // returns value of gradient

    // functions to inherit
    void get_var(string,double*);     // returns variable value
    void set_var(string,double);     // sets variable value    
    void return_name(string* inname){inname->assign(mmsname);};  // method: returns name
    
  }; // done with MMS base class

  // ------------------------------------------------------
  // ---------- all other mms classes ------------
  // ------------------------------------------------------
  // just a demo class
  class MASA_Test : public manufactured_solution 
  {
  private:
    double axp;
  public:
    MASA_Test(); // constructor
    
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
    int axp;
  public:
    heateq_2d_steady_const(); // constructor
    
  };
  
  class heateq_3d_steady_const : public manufactured_solution 
  {
  private:
    int axp;
  public:
    heateq_3d_steady_const(); // constructor
    
  };
  // ------------------------------------------------------
  // ---------- heat equation / unsteady / constant -------
  // ------------------------------------------------------

  class heateq_1d_unsteady_const : public manufactured_solution 
  {
  private:
    int axp;
  public:
    heateq_1d_unsteady_const(); // constructor
    
  };

  class heateq_2d_unsteady_const : public manufactured_solution 
  {
  private:
    int axp;
  public:
    heateq_2d_unsteady_const(); // constructor
    
  };
  
  class heateq_3d_unsteady_const : public manufactured_solution 
  {
  private:
    int axp;
  public:
    heateq_3d_unsteady_const(); // constructor
    
  };

  // ------------------------------------------------------
  // ---------- heat equation / unsteady / var ------------
  // ------------------------------------------------------

  class heateq_1d_unsteady_var : public manufactured_solution 
  {
  private:
    int axp;
  public:
    heateq_1d_unsteady_var(); // constructor
    
  };

  class heateq_2d_unsteady_var : public manufactured_solution 
  {
  private:
    int axp;
  public:
    heateq_2d_unsteady_var(); // constructor
    
  };
  
  class heateq_3d_unsteady_var : public manufactured_solution 
  {
  private:
    int axp;
  public:
    heateq_3d_unsteady_var(); // constructor
    
  };

  // ------------------------------------------------------
  // ---------- heat equation / steady / var ------------
  // ------------------------------------------------------

  class heateq_1d_steady_var : public manufactured_solution 
  {
  private:
    int axp;
  public:
    heateq_1d_steady_var(); // constructor
    
  };

  class heateq_2d_steady_var : public manufactured_solution 
  {
  private:
    int axp;
  public:
    heateq_2d_steady_var(); // constructor
    
  };
  
  class heateq_3d_steady_var : public manufactured_solution 
  {
  private:
    int axp;
  public:
    heateq_3d_steady_var(); // constructor
    
  };

  // ------------------------------------------------------
  // ---------- euler ------------
  // ------------------------------------------------------

  class euler_2d : public manufactured_solution
  {
    
  public:
    euler_2d(); // constructor
    
  };

  class euler_3d : public manufactured_solution
  {
    
  public:
    euler_3d(); // constructor
    
  };

  // ------------------------------------------------------
  // ---------- compressible navier stokes  ------------
  // ------------------------------------------------------

  class ns_compress_2d : public manufactured_solution
  {
    
  public:
    ns_compress_2d(); // constructor
    
  };

  class ns_compress_3d : public manufactured_solution
  {
    
  public:
    ns_compress_3d(); // constructor
    
  };

} // end MASA namespace
