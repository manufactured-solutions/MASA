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
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <vector>
#include <stdlib.h>

#include <config.h>

using std::cos;
using std::sin;
using std::pow;

// Macro for declaring MASA classes with all supported Scalar types
#define MASA_INSTANTIATE_ALL(my_class) template class my_class<double>; \
                                       template class my_class<long double>


// overload for pgi compilers
#ifdef portland_compiler
  long double pow(long double, double);
#endif

// overload for intel compilers
#ifdef intel_compiler
  long double pow(long double, double);
#endif

namespace MASA
{

  // masa map functions here
  // probably want to hide this from the user eventually
  int masa_map_solution  (std::string, std::string);
  int masa_map_temporal  (std::string, std::string);
  int masa_map_coeff     (std::string, std::string);
  int masa_map           (std::string, std::string);
  int masa_map_dimension (std::string, std::string);

  // masa map internal functions
  void uptolow(std::string&);
  void remove_line(std::string&);
  void remove_whitespace(std::string&);

  /*
   * -------------------------------------------------------------------------------------------
   *
   * manufactured_solution Base Class
   *
   * this is an abstract class
   * DO NOT EDIT THIS CLASS
   *
   * All manufactured solutions will inherit this class and it's methods
   *
   * -------------------------------------------------------------------------------------------
   */

  template <typename Scalar>
  class manufactured_solution
  {

  private:
    Scalar Tan;                          // analytical solution
    Scalar Q;                            // source term
    Scalar gradT;                        // gradient

  protected:
    static const Scalar MASA_VAR_DEFAULT;
    Scalar dummy;
    int num_vars;
    int num_vec;

    std::map<std::string,int> varmap;    // map to each variable
    std::vector<Scalar*>  vararr;        // arr of pointers to each variable

    std::map<std::string,int> vecmap;           // map to each vector
    std::vector<std::vector<Scalar>* >  vecarr; // vector of pointers to each vector

    std::string mmsname;                 // the name of the manufactured solution
    int dimension;                       // dimension of the solution

  public:
    static const Scalar pi;
    static const Scalar PI;

    // functions to override
    virtual ~manufactured_solution(){};   // destructor
    virtual int init_var() = 0;           // inits all variables to selected values

  /*
   * -------------------------------------------------------------------------------------------
   *
   * analytical solution(s)
   *
   * -------------------------------------------------------------------------------------------
   */

    virtual Scalar eval_exact_t(Scalar)                       {std::cout << "MASA ERROR:: Analytical Solution (T) is unavailable or not properly loaded.\n"; return -1.33;}; // returns value of analytical solution
    virtual Scalar eval_exact_t(Scalar,Scalar)                {std::cout << "MASA ERROR:: Analytical Solution (T) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 2d problems
    virtual Scalar eval_exact_t(Scalar,Scalar,Scalar)         {std::cout << "MASA ERROR:: Analytical Solution (T) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems
    virtual Scalar eval_exact_t(Scalar,Scalar,Scalar,Scalar)  {std::cout << "MASA ERROR:: Analytical Solution (T) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 4d problems

    virtual Scalar eval_exact_u(Scalar)                       {std::cout << "MASA ERROR:: Analytical Solution (u) is unavailable or not properly loaded.\n"; return -1.33;}; // returns value of analytical solution
    virtual Scalar eval_exact_u(Scalar,Scalar)                {std::cout << "MASA ERROR:: Analytical Solution (u) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 2d problems
    virtual Scalar eval_exact_u(Scalar,Scalar,Scalar)         {std::cout << "MASA ERROR:: Analytical Solution (u) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems
    virtual Scalar eval_exact_u(Scalar,Scalar,Scalar,Scalar)  {std::cout << "MASA ERROR:: Analytical Solution (u) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 4d problems

    virtual Scalar eval_exact_v(Scalar)                       {std::cout << "MASA ERROR:: Analytical Solution (v) is unavailable for 1D problems.\n"; return -1.33;};        // returns value of analytical solution
    virtual Scalar eval_exact_v(Scalar,Scalar)                {std::cout << "MASA ERROR:: Analytical Solution (v) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 2d problems
    virtual Scalar eval_exact_v(Scalar,Scalar,Scalar)         {std::cout << "MASA ERROR:: Analytical Solution (v) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems
    virtual Scalar eval_exact_v(Scalar,Scalar,Scalar,Scalar)  {std::cout << "MASA ERROR:: Analytical Solution (v) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 4d problems

    virtual Scalar eval_exact_w(Scalar)                       {std::cout << "MASA ERROR:: Analytical Solution (w) is unavailable for 1d problems.\n"; return -1.33;};        // returns value of analytical solution
    virtual Scalar eval_exact_w(Scalar,Scalar)                {std::cout << "MASA ERROR:: Analytical Solution (w) is unavailable for 2d problems.\n"; return -1.33;};        // overloaded for 2d problems
    virtual Scalar eval_exact_w(Scalar,Scalar,Scalar)         {std::cout << "MASA ERROR:: Analytical Solution (w) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems
    virtual Scalar eval_exact_w(Scalar,Scalar,Scalar,Scalar)  {std::cout << "MASA ERROR:: Analytical Solution (w) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 4d problems

    virtual Scalar eval_exact_p(Scalar)                       {std::cout << "MASA ERROR:: Analytical Solution (e) is unavailable or not properly loaded.\n"; return -1.33;}; // returns value of analytical solution
    virtual Scalar eval_exact_p(Scalar,Scalar)                {std::cout << "MASA ERROR:: Analytical Solution (e) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 2d problems
    virtual Scalar eval_exact_p(Scalar,Scalar,Scalar)         {std::cout << "MASA ERROR:: Analytical Solution (e) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems
    virtual Scalar eval_exact_p(Scalar,Scalar,Scalar,Scalar)  {std::cout << "MASA ERROR:: Analytical Solution (e) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 4d problems

    virtual Scalar eval_exact_phi(Scalar,Scalar)              {std::cout << "MASA ERROR:: Analytical Solution (phi) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 2d problems

    virtual Scalar eval_exact_rho(Scalar)                      {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;}; // returns value of analytical solution
    virtual Scalar eval_exact_rho(Scalar,Scalar)               {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 2d problems
    virtual Scalar eval_exact_rho(Scalar,Scalar,Scalar)        {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 3d problems
    virtual Scalar eval_exact_rho(Scalar,Scalar,Scalar,Scalar) {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;}; // overloaded for 4d problems

    virtual Scalar eval_exact_nu (Scalar,Scalar)               {std::cout << "MASA ERROR:: Analytical Solution (nu) is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_exact_nu (Scalar,Scalar,Scalar)        {std::cout << "MASA ERROR:: Analytical Solution (nu) is unavailable or not properly loaded.\n"; return -1.33;};

    virtual Scalar eval_exact_rho_N(Scalar)            {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_exact_rho_N2(Scalar)           {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_exact_rho_C(Scalar)            {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_exact_rho_C3(Scalar)           {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_exact_rho_C(Scalar,Scalar)     {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_exact_rho_C3(Scalar,Scalar)    {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_exact_rho_C(Scalar,Scalar,Scalar)  {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_exact_rho_C3(Scalar,Scalar,Scalar) {std::cout << "MASA ERROR:: Analytical Solution (rho) is unavailable or not properly loaded.\n"; return -1.33;};

    virtual Scalar eval_exact_u_boundary(Scalar)            {std::cout << "MASA ERROR:: Analytical Solution (boundary) is unavailable or not properly loaded.\n"; return -1.33;};

  /*
   * -------------------------------------------------------------------------------------------
   *
   * source terms
   *
   * -------------------------------------------------------------------------------------------
   */

    virtual Scalar eval_q_f(Scalar, Scalar)               {std::cout << "MASA ERROR:: Solution (f) has not been properly loaded.\n"; return -1.33;};  // returns value of source term (f)

    virtual Scalar eval_q_t(Scalar)                       {std::cout << "MASA ERROR:: Solution has not been properly loaded.\n"; return -1.33;};  // returns value of source term (temp)
    virtual Scalar eval_q_t(Scalar, Scalar)               {std::cout << "MASA ERROR:: Solution has not been properly loaded.\n"; return -1.33;};  // returns value of source term (temp)
    virtual Scalar eval_q_t(Scalar,Scalar,Scalar)         {std::cout << "MASA ERROR:: Solution has not been properly loaded.\n"; return -1.33;};  // returns value of source term (temp)
    virtual Scalar eval_q_t(Scalar,Scalar,Scalar,Scalar)  {std::cout << "MASA ERROR:: Solution has not been properly loaded.\n"; return -1.33;};  // returns value of source term (x,y,z,t)

    virtual Scalar eval_q_u(Scalar)                       {std::cout << "MASA ERROR:: Source term (u) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (u)
    virtual Scalar eval_q_u(Scalar,Scalar)                {std::cout << "MASA ERROR:: Source term (u) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 2d problems
    virtual Scalar eval_q_u(Scalar,Scalar,Scalar)         {std::cout << "MASA ERROR:: Source term (u) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 3d problems
    virtual Scalar eval_q_u(Scalar,Scalar,Scalar,Scalar)  {std::cout << "MASA ERROR:: Source term (u) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 4d problems

    virtual Scalar eval_q_v(Scalar)                       {std::cout << "MASA ERROR:: Source term (v) is unavailable for 1d problems -- eval_q_v has too few arguments.\n"; return -1.33;}; // returns value of source term (v)
    virtual Scalar eval_q_v(Scalar,Scalar)                {std::cout << "MASA ERROR:: Source term (v) is unavailable or not properly loaded.\n"; return -1.33;};                            // overloaded for 2d problems
    virtual Scalar eval_q_v(Scalar,Scalar,Scalar)         {std::cout << "MASA ERROR:: Source term (v) is unavailable or not properly loaded.\n"; return -1.33;};                            // overloaded for 3d problems
    virtual Scalar eval_q_v(Scalar,Scalar,Scalar,Scalar)  {std::cout << "MASA ERROR:: Source term (v) is unavailable or not properly loaded.\n"; return -1.33;};                            // overloaded for 4d problems

    virtual Scalar eval_q_w(Scalar)                       {std::cout << "MASA ERROR:: Source term (w) is unavailable for 1d problems -- eval_q_w has too few arguments.\n"; return -1.33;};  // returns value of source term (w)
    virtual Scalar eval_q_w(Scalar,Scalar)                {std::cout << "MASA ERROR:: Source term (w) is unavailable for 2d problems -- eval_q_w has too few arguments.\n"; return -1.33;};  // overloaded for 2d problems
    virtual Scalar eval_q_w(Scalar,Scalar,Scalar)         {std::cout << "MASA ERROR:: Source Term (w) is unavailable or not properly loaded.\n"; return -1.33;};                             // overloaded for 3d problems
    virtual Scalar eval_q_w(Scalar,Scalar,Scalar,Scalar)  {std::cout << "MASA ERROR:: Source Term (w) is unavailable or not properly loaded.\n"; return -1.33;};                             // overloaded for 4d problems

    virtual Scalar eval_q_e(Scalar)                       {std::cout << "MASA ERROR:: Source Term (e) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (energy)
    virtual Scalar eval_q_e(Scalar,Scalar (*)(Scalar))    {std::cout << "MASA ERROR:: Source Term (e) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (energy)
    virtual Scalar eval_q_e(Scalar,Scalar)                {std::cout << "MASA ERROR:: Source Term (e) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 2d problems
    virtual Scalar eval_q_e(Scalar,Scalar,Scalar)         {std::cout << "MASA ERROR:: Source Term (e) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 3d problems
    virtual Scalar eval_q_e(Scalar,Scalar,Scalar,Scalar)  {std::cout << "MASA ERROR:: Source Term (e) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 4d problems

    virtual Scalar eval_q_rho(Scalar)                     {std::cout << "MASA ERROR:: Source Term (rho) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (density)
    virtual Scalar eval_q_rho(Scalar,Scalar)              {std::cout << "MASA ERROR:: Source Term (rho) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 2d problems
    virtual Scalar eval_q_rho(Scalar,Scalar,Scalar)       {std::cout << "MASA ERROR:: Source Term (rho) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 3d problems
    virtual Scalar eval_q_rho(Scalar,Scalar,Scalar,Scalar){std::cout << "MASA ERROR:: Source Term (rho) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 4d problems

    virtual Scalar eval_q_nu (Scalar,Scalar){std::cout << "MASA ERROR:: Source Term (nu) is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_q_nu (Scalar,Scalar,Scalar){std::cout << "MASA ERROR:: Source Term (nu) is unavailable or not properly loaded.\n"; return -1.33;};

    virtual Scalar eval_q_rho_u(Scalar)                     {std::cout << "MASA ERROR:: Source Term (rho*u) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (density*u)
    virtual Scalar eval_q_rho_u(Scalar,Scalar)              {std::cout << "MASA ERROR:: Source Term (rho*u) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 2d problems
    virtual Scalar eval_q_rho_u(Scalar,Scalar,Scalar)       {std::cout << "MASA ERROR:: Source Term (rho*u) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 3d problems
    virtual Scalar eval_q_rho_u(Scalar,Scalar,Scalar,Scalar){std::cout << "MASA ERROR:: Source Term (rho*u) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 4d problems

    virtual Scalar eval_q_rho_v(Scalar)                     {std::cout << "MASA ERROR:: Source Term (rho*v) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (density*v)
    virtual Scalar eval_q_rho_v(Scalar,Scalar)              {std::cout << "MASA ERROR:: Source Term (rho*v) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 2d problems
    virtual Scalar eval_q_rho_v(Scalar,Scalar,Scalar)       {std::cout << "MASA ERROR:: Source Term (rho*v) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 3d problems
    virtual Scalar eval_q_rho_v(Scalar,Scalar,Scalar,Scalar){std::cout << "MASA ERROR:: Source Term (rho*v) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 4d problems

    virtual Scalar eval_q_rho_w(Scalar)                     {std::cout << "MASA ERROR:: Source Term (rho*w) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (density*w)
    virtual Scalar eval_q_rho_w(Scalar,Scalar)              {std::cout << "MASA ERROR:: Source Term (rho*w) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 2d problems
    virtual Scalar eval_q_rho_w(Scalar,Scalar,Scalar)       {std::cout << "MASA ERROR:: Source Term (rho*w) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 3d problems
    virtual Scalar eval_q_rho_w(Scalar,Scalar,Scalar,Scalar){std::cout << "MASA ERROR:: Source Term (rho*w) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 4d problems

    virtual Scalar eval_q_u_boundary (Scalar)        {std::cout << "MASA ERROR:: Source Term (u_boundary) is unavailable or not properly loaded.\n"; return -1.33;};

    virtual Scalar eval_q_rho_e(Scalar)                      {std::cout << "MASA ERROR:: Source Term (rho*e) is unavailable or not properly loaded.\n"; return -1.33;};  // returns value of source term (density*e)
    virtual Scalar eval_q_rho_e(Scalar,Scalar)               {std::cout << "MASA ERROR:: Source Term (rho*e) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 2d problems
    virtual Scalar eval_q_rho_e(Scalar,Scalar,Scalar)        {std::cout << "MASA ERROR:: Source Term (rho*e) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 3d problems
    virtual Scalar eval_q_rho_e(Scalar,Scalar,Scalar,Scalar) {std::cout << "MASA ERROR:: Source Term (rho*e) is unavailable or not properly loaded.\n"; return -1.33;};  // overloaded for 4d problems

    virtual Scalar eval_q_rho_N (Scalar,Scalar (*)(Scalar))    {std::cout << "MASA ERROR:: Source Term (N )    is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_q_rho_N2(Scalar,Scalar (*)(Scalar))    {std::cout << "MASA ERROR:: Source Term (N2)    is unavailable or not properly loaded.\n"; return -1.33;};

    virtual Scalar eval_q_C (Scalar)    {std::cout << "MASA ERROR:: Source Term (N )    is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_q_C3(Scalar)    {std::cout << "MASA ERROR:: Source Term (N2)    is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_q_rho_C (Scalar)    {std::cout << "MASA ERROR:: Source Term (N )    is unavailable or not properly loaded.\n"; return -1.33;};
    virtual Scalar eval_q_rho_C3(Scalar)    {std::cout << "MASA ERROR:: Source Term (N2)    is unavailable or not properly loaded.\n"; return -1.33;};

  /*
   * -------------------------------------------------------------------------------------------
   *
   * Gradient Terms
   *
   * -------------------------------------------------------------------------------------------
   */

    virtual Scalar eval_g_t(Scalar)                          {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 1d gradient
    virtual Scalar eval_g_t(Scalar,Scalar,int)               {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 2d gradient
    virtual Scalar eval_g_t(Scalar,Scalar,Scalar,int)        {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d gradient
    virtual Scalar eval_g_t(Scalar,Scalar,Scalar,Scalar,int) {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d, time-varying gradient

    virtual Scalar eval_g_u(Scalar)                          {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 1d gradient
    virtual Scalar eval_g_u(Scalar,Scalar,int)               {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 2d gradient
    virtual Scalar eval_g_u(Scalar,Scalar,Scalar,int)        {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d gradient
    virtual Scalar eval_g_u(Scalar,Scalar,Scalar,Scalar,int) {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d, time-varying gradient

    virtual Scalar eval_g_v(Scalar)                          {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 1d gradient
    virtual Scalar eval_g_v(Scalar,Scalar,int)               {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 2d gradient
    virtual Scalar eval_g_v(Scalar,Scalar,Scalar,int)        {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d gradient
    virtual Scalar eval_g_v(Scalar,Scalar,Scalar,Scalar,int) {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d, time-varying gradient

    virtual Scalar eval_g_w(Scalar)                          {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 1d gradient
    virtual Scalar eval_g_w(Scalar,Scalar,int)               {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 2d gradient
    virtual Scalar eval_g_w(Scalar,Scalar,Scalar,int)        {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d gradient
    virtual Scalar eval_g_w(Scalar,Scalar,Scalar,Scalar,int) {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d, time-varying gradient

    virtual Scalar eval_g_p(Scalar)                          {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 1d gradient
    virtual Scalar eval_g_p(Scalar,Scalar,int)               {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 2d gradient
    virtual Scalar eval_g_p(Scalar,Scalar,Scalar,int)        {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d gradient
    virtual Scalar eval_g_p(Scalar,Scalar,Scalar,Scalar,int) {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d, time-varying gradient

    virtual Scalar eval_g_rho(Scalar)                          {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 1d gradient
    virtual Scalar eval_g_rho(Scalar,Scalar,int)               {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 2d gradient
    virtual Scalar eval_g_rho(Scalar,Scalar,Scalar,int)        {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d gradient
    virtual Scalar eval_g_rho(Scalar,Scalar,Scalar,Scalar,int) {std::cout << "MASA ERROR:: gradient is unavailable or not properly loaded.\n";   return -1.33;};  // returns value of 3d, time-varying gradient

  /*
   * -------------------------------------------------------------------------------------------
   *
   * member functions solution classes will inherit (not override)
   *
   * -------------------------------------------------------------------------------------------
   */
    manufactured_solution();                                     // constructor
    int purge_var();                                             // dump defaults
    Scalar pass_function(Scalar (*)(Scalar),Scalar);
    int set_var(std::string,Scalar);                             // sets variable value
    int set_vec(std::string,std::vector<Scalar>&);               // sets vector value
    int register_var(std::string, Scalar*);                      // this registers a variable
    int register_vec(std::string, std::vector<Scalar>& );        // this registers a vector

    int sanity_check();                                          // checks that all variables to the class have been initalized
    int poly_test();                                             // regression method for poly class (see below)
    int get_vec(std::string,std::vector<Scalar>&);               // returns vector values
    Scalar get_var(std::string);                                 // returns variable value
    int display_var();                                           // print all variable names and values
    int display_vec();                                           // print all variable names and values
    void return_name(std::string* inname){inname->assign(mmsname);};  // method: returns name
    void return_dim (int* indim)    {*indim=dimension;};         // method: returns dimension of solution

  /*
   * -------------------------------------------------------------------------------------------
   *
   * smasa: functions (stochastic masa)
   *
   * -------------------------------------------------------------------------------------------
   */

    virtual Scalar eval_post_mean () {std::cout << "SMASA ERROR:: posterior mean is unavailable.\n";return -1.33;};
    virtual Scalar eval_post_var  () {std::cout << "SMASA ERROR:: posterior variance is unavailable.\n";return -1.33;};

    virtual Scalar eval_cen_mom      (int /*x*/) {std::cout << "SMASA ERROR:: central moment is unavailable or not properly loaded.\n";return -1.33;};
    virtual Scalar eval_likelyhood   (Scalar){std::cout << "SMASA ERROR:: likelyhood is unavailable or not properly loaded.\n";return -1.33;};
    virtual Scalar eval_loglikelyhood(Scalar){std::cout << "SMASA ERROR:: loglikelyhood is unavailable or not properly loaded.\n";return -1.33;};
    virtual Scalar eval_prior        (Scalar){std::cout << "SMASA ERROR:: prior is unavailable or not properly loaded.\n";return -1.33;};
    virtual Scalar eval_posterior    (Scalar){std::cout << "SMASA ERROR:: posterior is unavailable or not properly loaded.\n";return -1.33;};

  /*
   * -------------------------------------------------------------------------------------------
   *
   * smasa: distance between probability distributions
   *
   * -------------------------------------------------------------------------------------------
   */

    virtual Scalar eval_hellinger (Scalar) {std::cout << "SMASA ERROR:: Hellinger distance is unavailable or not properly loaded.\n";return -1.33;};
    virtual Scalar eval_kolmogorov(Scalar) {std::cout << "SMASA ERROR:: Kolmogorov distance is unavailable or not properly loaded.\n";return -1.33;};

  }; // done with MMS base class

  /*
   * -------------------------------------------------------------------------------------------
   *
   * Polynomial Base Class
   *
   * Blatantly stealing paul bauman's polynomial class definitions in the name of science
   *
   * In addition to inheriting manufactured class, must inherit polynomial
   *
   * -------------------------------------------------------------------------------------------
   */
  template <typename Scalar>
  class Polynomial
  {
  public:

    void set_coeffs( const std::vector<Scalar> & );

    // Evaluates polynomial.
    Scalar operator()( const Scalar &, int *) const;

    // Evaluates polynomial and deriviatives up to order specified by user.
    void eval_derivs( const Scalar, const int, std::vector<Scalar> & ) const;
    Scalar get_coeffs( const int & ) const;

  protected:

    // We assume that the coefficents are ordered as follows:
    // y(x) = a0 + a1*x + a2*x^2 + ... + an*x^n
    // so that there should be n+1 coeff so that
    // coeffs[0] = a0
    // coeffs[1] = a1
    // and so on.
    std::vector<Scalar> coeffs;

  };

  // ------------------------------------------------------
  // ---------- all other mms classes ---------------------
  // ------------------------------------------------------

  // just a demo class
  template <typename Scalar>
  class masa_test_function : public manufactured_solution<Scalar>
  {
  private:
    Scalar demo_var_2;
    Scalar demo_var_3;
  public:
    masa_test_function(); // constructor
    int init_var();        // default problem values
    int poly_test();
    //Scalar eval_q_t(Scalar);
  }; // done with masa_test


  // class with no source/analytical terms, to test virtual function defaults
  template <typename Scalar>
  class masa_uninit : public manufactured_solution<Scalar>
  {
  private:

  public:
    masa_uninit();  // constructor
    int init_var(); // default problem values
  };

  // ------------------------------------------------------
  // ---------- heat equation /steady / constant ----------
  // ------------------------------------------------------

  template <typename Scalar>
  class heateq_1d_steady_const : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar k_0;
  public:
    heateq_1d_steady_const(); // constructor
    int init_var();        // default problem values
    Scalar eval_q_t (Scalar);  // source term evaluator
    Scalar eval_exact_t(Scalar);   //analytical solution

  };

  template <typename Scalar>
  class heateq_2d_steady_const : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar k_0;
    Scalar B_y;

  public:
    heateq_2d_steady_const();       // constructor
    int init_var();        // default problem values
    Scalar eval_q_t (Scalar,Scalar); // source term evaluator
    Scalar eval_exact_t(Scalar,Scalar); // analytical term evaluator
  };

  template <typename Scalar>
  class heateq_3d_steady_const : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar k_0;
    Scalar B_y;
    Scalar C_z;

  public:
    heateq_3d_steady_const(); // constructor
    int init_var();        // default problem values
    Scalar eval_q_t (Scalar,Scalar,Scalar); //evaluate source term
    Scalar eval_exact_t(Scalar,Scalar,Scalar); // analytical term evaluator
  };
  // ------------------------------------------------------
  // ---------- heat equation / unsteady / constant -------
  // ------------------------------------------------------

  template <typename Scalar>
  class heateq_1d_unsteady_const : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar A_t;
    Scalar D_t;
    Scalar k_0;
    Scalar cp_0;
    Scalar rho;

  public:
    heateq_1d_unsteady_const(); // constructor
    int init_var();        // default problem values
    Scalar eval_q_t(Scalar,Scalar); // needs x,t
  };

  template <typename Scalar>
  class heateq_2d_unsteady_const : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar A_t;
    Scalar B_y;
    Scalar B_t;
    Scalar D_t;
    Scalar rho;
    Scalar k_0;
    Scalar cp_0;

  public:
    heateq_2d_unsteady_const();                // constructor
    int init_var();                            // default problem values
    Scalar eval_q_t(Scalar,Scalar,Scalar);     // needs x,y,t
    Scalar eval_exact_t(Scalar,Scalar,Scalar); // analytical term evaluator
  };

  template <typename Scalar>
  class heateq_3d_unsteady_const : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar A_t;
    Scalar B_y;
    Scalar B_t;
    Scalar C_z;
    Scalar C_t;
    Scalar D_t;
    Scalar k_0;
    Scalar cp_0;
    Scalar rho;

  public:
    heateq_3d_unsteady_const(); // constructor
    int init_var();        // default problem values
    Scalar eval_q_t(Scalar,Scalar,Scalar,Scalar); // needs x,y,z,t
  };

  // ------------------------------------------------------
  // ---------- heat equation / unsteady / var ------------
  // ------------------------------------------------------

  template <typename Scalar>
  class heateq_1d_unsteady_var : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar A_t;
    Scalar D_t;
    Scalar rho;
    Scalar k_0;
    Scalar k_1;
    Scalar k_2;
    Scalar cp_0;
    Scalar cp_1;
    Scalar cp_2;

  public:
    heateq_1d_unsteady_var(); // constructor
    int init_var();        // default problem values
    Scalar eval_q_t(Scalar,Scalar); // needs x,t
  };

  template <typename Scalar>
  class heateq_2d_unsteady_var : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar A_t;
    Scalar B_y;
    Scalar B_t;
    Scalar D_t;
    Scalar rho;
    Scalar k_0;
    Scalar k_1;
    Scalar k_2;
    Scalar cp_0;
    Scalar cp_1;
    Scalar cp_2;

  public:
    heateq_2d_unsteady_var(); // constructor
    int init_var();        // default problem values
    Scalar eval_q_t(Scalar,Scalar,Scalar); // needs x,y,t
  };

  template <typename Scalar>
  class heateq_3d_unsteady_var : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar A_t;
    Scalar B_y;
    Scalar B_t;
    Scalar C_z;
    Scalar C_t;
    Scalar D_t;
    Scalar rho;
    Scalar k_0;
    Scalar k_1;
    Scalar k_2;
    Scalar cp_0;
    Scalar cp_1;
    Scalar cp_2;

  public:
    heateq_3d_unsteady_var(); // constructor
    int init_var();        // default problem values
    Scalar eval_q_t(Scalar,Scalar,Scalar,Scalar); // needs x,y,z,t
  };

  // ------------------------------------------------------
  // ---------- heat equation / steady / var --------------
  // ------------------------------------------------------

  template <typename Scalar>
  class heateq_1d_steady_var : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar k_0;
    Scalar k_1;
    Scalar k_2;

  public:
    heateq_1d_steady_var(); // constructor
    int init_var();        // default problem values
    Scalar eval_q_t(Scalar); // needs x
  };

  template <typename Scalar>
  class heateq_2d_steady_var : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar k_0;
    Scalar k_1;
    Scalar k_2;
    Scalar B_y;

  public:
    heateq_2d_steady_var(); // constructor
    int init_var();        // default problem values
    Scalar eval_q_t(Scalar,Scalar); // needs x,y
  };

  template <typename Scalar>
  class heateq_3d_steady_var : public manufactured_solution<Scalar>
  {
  private:
    Scalar A_x;
    Scalar k_0;
    Scalar k_1;
    Scalar k_2;
    Scalar B_y;
    Scalar C_z;

  public:
    heateq_3d_steady_var(); // constructor
    int init_var();        // default problem values
    Scalar eval_q_t(Scalar,Scalar,Scalar); // needs x,y,z

  };

  // ------------------------------------------------------
  // ------------------------- euler ----------------------
  // ------------------------------------------------------
  template <typename Scalar>
  class euler_1d : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar k;                             // Boltzmanns constant

    Scalar u_0;
    Scalar u_x;
    Scalar rho_0;
    Scalar rho_x;
    Scalar p_0;
    Scalar p_x;
    Scalar a_px;
    Scalar a_rhox;
    Scalar a_ux;
    Scalar Gamma;
    Scalar mu;
    Scalar L;


  public:
    euler_1d(); // constructor
    int init_var();          // default problem values

    Scalar eval_q_rho_u (Scalar); // source terms
    Scalar eval_q_rho_e (Scalar);
    Scalar eval_q_rho   (Scalar);

    Scalar eval_exact_u  (Scalar); // analytical
    Scalar eval_exact_p  (Scalar);
    Scalar eval_exact_rho(Scalar);

    Scalar eval_g_u  (Scalar);   // gradient of source term
    Scalar eval_g_p  (Scalar);
    Scalar eval_g_rho(Scalar);

  };

  template <typename Scalar>
  class euler_2d : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar k;     // Boltzmanns constant

    Scalar u_0;
    Scalar u_x;
    Scalar u_y;
    Scalar v_0;
    Scalar v_x;
    Scalar v_y;
    Scalar rho_0;
    Scalar rho_x;
    Scalar rho_y;
    Scalar p_0;
    Scalar p_x;
    Scalar p_y;
    Scalar a_px;
    Scalar a_py;
    Scalar a_rhox;
    Scalar a_rhoy;
    Scalar a_ux;
    Scalar a_uy;
    Scalar a_vx;
    Scalar a_vy;
    Scalar Gamma;
    Scalar mu;
    Scalar L;

  public:
    euler_2d(); // constructor
    int init_var();        // default problem values

    Scalar eval_q_rho_u (Scalar,Scalar);
    Scalar eval_q_rho_v (Scalar,Scalar);
    Scalar eval_q_rho_e (Scalar,Scalar);
    Scalar eval_q_rho   (Scalar,Scalar);

    Scalar eval_exact_u  (Scalar,Scalar); // analytical
    Scalar eval_exact_v  (Scalar,Scalar);
    Scalar eval_exact_p  (Scalar,Scalar);
    Scalar eval_exact_rho(Scalar,Scalar);

    Scalar eval_g_u(Scalar,Scalar,int);   // gradient of source term
    Scalar eval_g_v(Scalar,Scalar,int);
    Scalar eval_g_p(Scalar,Scalar,int);

    Scalar eval_g_rho(Scalar,Scalar,int);

  };

  template <typename Scalar>
  class euler_3d : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar k;     // Boltzmanns constant

    Scalar u_0;
    Scalar u_x;
    Scalar u_y;
    Scalar u_z;
    Scalar v_0;
    Scalar v_x;
    Scalar v_y;
    Scalar v_z;
    Scalar w_0;
    Scalar w_x;
    Scalar w_y;
    Scalar w_z;
    Scalar rho_0;
    Scalar rho_x;
    Scalar rho_y;
    Scalar rho_z;
    Scalar p_0;
    Scalar p_x;
    Scalar p_y;
    Scalar p_z;
    Scalar a_px;
    Scalar a_py;
    Scalar a_pz;
    Scalar a_rhox;
    Scalar a_rhoy;
    Scalar a_rhoz;
    Scalar a_ux;
    Scalar a_uy;
    Scalar a_uz;
    Scalar a_vx;
    Scalar a_vy;
    Scalar a_vz;
    Scalar a_wx;
    Scalar a_wy;
    Scalar a_wz;
    Scalar mu;
    Scalar Gamma;
    Scalar L;

  public:
    euler_3d(); // constructor
    int init_var();        // default problem values

    Scalar eval_q_rho_u  (Scalar,Scalar,Scalar); // source terms
    Scalar eval_q_rho_v  (Scalar,Scalar,Scalar);
    Scalar eval_q_rho_w  (Scalar,Scalar,Scalar);
    Scalar eval_q_rho_e  (Scalar,Scalar,Scalar);
    Scalar eval_q_rho    (Scalar,Scalar,Scalar);

    Scalar eval_exact_u  (Scalar,Scalar,Scalar); // analytical
    Scalar eval_exact_v  (Scalar,Scalar,Scalar);
    Scalar eval_exact_w  (Scalar,Scalar,Scalar);
    Scalar eval_exact_q  (Scalar,Scalar,Scalar);
    Scalar eval_exact_p  (Scalar,Scalar,Scalar);
    Scalar eval_exact_rho(Scalar,Scalar,Scalar);

    Scalar eval_g_u(Scalar,Scalar,Scalar,int);   // gradient of source term
    Scalar eval_g_v(Scalar,Scalar,Scalar,int);
    Scalar eval_g_w(Scalar,Scalar,Scalar,int);

    Scalar eval_g_p  (Scalar,Scalar,Scalar,int);
    Scalar eval_g_rho(Scalar,Scalar,Scalar,int);

  };
  // ------------------------------------------------------
  // ------------------------- euler - transient ----------
  // ------------------------------------------------------

  template <typename Scalar>
  class euler_transient_1d : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar k;                             // Boltzmanns constant

    Scalar u_0;
    Scalar u_x;
    Scalar u_t;
    Scalar a_ux;
    Scalar a_ut;

    Scalar rho_0;
    Scalar rho_x;
    Scalar rho_t;
    Scalar a_rhox;
    Scalar a_rhot;

    Scalar p_0;
    Scalar p_x;
    Scalar p_t;
    Scalar a_px;
    Scalar a_pt;

    Scalar Gamma;
    Scalar mu;
    Scalar L;

  public:
    euler_transient_1d(); // constructor
    int init_var();          // default problem values

    Scalar eval_q_rho_u (Scalar,Scalar);
    Scalar eval_q_rho_e (Scalar,Scalar);
    Scalar eval_q_rho   (Scalar,Scalar);

    Scalar eval_exact_u    (Scalar,Scalar);
    Scalar eval_exact_p    (Scalar,Scalar);
    Scalar eval_exact_rho  (Scalar,Scalar);

  };



  // ------------------------------------------------------
  // ----------------   euler + chem     ------------------
  // ------------------------------------------------------

  template <typename Scalar>
  class euler_chem_1d : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar R;
    Scalar Cf1_N;
    Scalar Cf1_N2;
    Scalar etaf1_N;
    Scalar etaf1_N2;
    Scalar Ea_N;
    Scalar Ea_N2;
    Scalar R_N;
    Scalar R_N2;
    Scalar theta_v_N2;
    Scalar M_N;
    Scalar h0_N;
    Scalar h0_N2;

    Scalar L;
    Scalar u_x;
    Scalar u_0;
    Scalar a_ux;

    Scalar rho_N_0;
    Scalar rho_N_x;
    Scalar a_rho_N_x;

    Scalar rho_N2_0;
    Scalar rho_N2_x;
    Scalar a_rho_N2_x;

    Scalar T_0;
    Scalar T_x;
    Scalar a_Tx;

  public:
    euler_chem_1d(); // constructor
    int init_var();          // default problem values

    Scalar eval_q_rho_u  (Scalar);
    Scalar eval_q_rho_e  (Scalar);
    Scalar eval_q_rho_N  (Scalar,Scalar (*)(Scalar));
    Scalar eval_q_rho_N2 (Scalar,Scalar (*)(Scalar));

    Scalar eval_exact_t      (Scalar);
    Scalar eval_exact_u      (Scalar);
    Scalar eval_exact_rho    (Scalar);
    Scalar eval_exact_rho_N  (Scalar);
    Scalar eval_exact_rho_N2 (Scalar);

  };

  // ------------------------------------------------------
  // ---------------- axisymmetric euler ------------------
  // ------------------------------------------------------
  template <typename Scalar>
  class axi_euler : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar R;                             // (or is this the ideal gas constant?) ratio of specific heat capacities, defined in constructor

    Scalar p_0;
    Scalar p_r;
    Scalar p_z;
    Scalar rho_0;
    Scalar rho_r;
    Scalar rho_z;
    Scalar u_r;
    Scalar u_z;
    Scalar w_0;
    Scalar w_r;
    Scalar w_z;
    Scalar a_pr;
    Scalar a_pz;
    Scalar a_rhor;
    Scalar a_rhoz;
    Scalar a_ur;
    Scalar a_uz;
    Scalar a_wr;
    Scalar a_wz;
    Scalar L;
    Scalar mu;
    Scalar Gamma;

  public:
    axi_euler(); // constructor
    int init_var();          // default problem values

    Scalar eval_q_rho_u (Scalar,Scalar); // radial velocity
    Scalar eval_q_rho_w (Scalar,Scalar); // axial
    Scalar eval_q_rho_e (Scalar,Scalar);
    Scalar eval_q_rho   (Scalar,Scalar);

    Scalar eval_exact_u  (Scalar,Scalar); // analytical
    Scalar eval_exact_w  (Scalar,Scalar);
    Scalar eval_exact_p  (Scalar,Scalar);
    Scalar eval_exact_rho(Scalar,Scalar);

  };

  // ------------------------------------------------------
  // ---------- axisymmetric compressible navier stokes ---
  // ------------------------------------------------------
  template <typename Scalar>
  class axi_cns : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar R;                             // (or is this the ideal gas constant?) ratio of specific heat capacities, defined in constructor
    Scalar k;                             // Boltzmanns constant

    Scalar p_0;
    Scalar p_1;
    Scalar rho_0;
    Scalar rho_1;
    Scalar u_1;
    Scalar w_0;
    Scalar w_1;
    Scalar a_pr;
    Scalar a_pz;
    Scalar a_rhor;
    Scalar a_rhoz;
    Scalar a_ur;
    Scalar a_uz;
    Scalar a_wr;
    Scalar a_wz;
    Scalar L;
    Scalar mu;
    Scalar Gamma;

  public:
    axi_cns(); // constructor
    int init_var();          // default problem values

    Scalar eval_q_rho_u (Scalar,Scalar); // radial velocity
    Scalar eval_q_rho_w (Scalar,Scalar); // axial
    Scalar eval_q_rho_e (Scalar,Scalar);
    Scalar eval_q_rho   (Scalar,Scalar);

    Scalar eval_exact_u  (Scalar,Scalar); // analytical
    Scalar eval_exact_w  (Scalar,Scalar);
    Scalar eval_exact_p  (Scalar,Scalar);
    Scalar eval_exact_rho(Scalar,Scalar);

  };

  // ------------------------------------------------------
  // --------------------- sod 1d -------------------------
  // ------------------------------------------------------
  template <typename Scalar>
  class sod_1d : public manufactured_solution<Scalar>
  {
    Scalar Gamma;
    Scalar mu;
    Scalar pl, pr, rhol, rhor, cl, cr;

  public:
    sod_1d(); // constructor
    int init_var();          // default problem values

    Scalar eval_q_rho   (Scalar,Scalar);
    Scalar eval_q_p     (Scalar,Scalar);
    Scalar eval_q_rho_u (Scalar,Scalar);
    Scalar func         (Scalar);  // helper function
    Scalar rtbis        (Scalar,Scalar,Scalar,int);
    Scalar eval_q_t     ();
    Scalar eval_q_t     (Scalar x);
  };

  // ------------------------------------------------------
  // ---------- RANS: Spelart Alamaras (Channel) ----------
  // ------------------------------------------------------
  template <typename Scalar>
  class   fans_sa_steady_wall_bounded : public manufactured_solution<Scalar>
  {

    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    // fluid properties
    Scalar mu;
    Scalar R;
    Scalar p_0;
    Scalar Pr;
    Scalar Pr_t;

    // SA calibration model
    Scalar eta1;
    Scalar eta_v;
    Scalar kappa;
    Scalar sigma;
    Scalar c_b1;
    Scalar c_b2;
    Scalar c_v1;
    Scalar c_v2;
    Scalar c_v3;
    Scalar c_w2;
    Scalar c_w3;

    // manufactured solutions
    Scalar T_inf;
    Scalar M_inf;
    Scalar r_T;
    Scalar Gamma;
    Scalar alpha;
    Scalar C_cf;
    Scalar C;
    Scalar b;

    // constants that are computed from parameters
    Scalar C1;
    Scalar u_inf;
    Scalar rho_inf;
    Scalar T_aw;
    Scalar rho_w;
    Scalar A;
    Scalar F_c;
    Scalar nu_w;
    Scalar Re_x;
    Scalar c_f;
    Scalar u_tau;
    Scalar u_eq_plus;
    Scalar y_plus;
    Scalar u_eq;
    Scalar U;
    Scalar V;
    Scalar T;
    Scalar RHO;
    Scalar NU_SA;
    Scalar chi;
    Scalar f_v1;
    Scalar mu_t;
    Scalar d_ueqplus_yplus;

    Scalar c_w1;
    Scalar d;
    Scalar Omega;
    Scalar Sm1;
    Scalar Sm;
    Scalar Sm2;
    Scalar Sm_orig;
    Scalar g;
    Scalar r;
    Scalar S_sa;
    Scalar cp;

    Scalar sigma_sa;
    Scalar f_w;
    Scalar f_v2;

    Scalar D2ueqDx2;
    Scalar D2ueqDy2;
    Scalar D2uDx2;
    Scalar D2uDy2;
    Scalar D2vDxy;
    
    Scalar D2vDx2;
    Scalar D2vDy2;
    Scalar D2TDx2;
    Scalar D2TDy2;
    Scalar D2uDxy;

  public:
    fans_sa_steady_wall_bounded(); // constructor
    int init_var();
    Scalar update(Scalar,Scalar);

    Scalar eval_q_rho_u (Scalar,Scalar);
    Scalar eval_q_rho_v (Scalar,Scalar);
    Scalar eval_q_rho_e (Scalar,Scalar);
    Scalar eval_q_rho   (Scalar,Scalar);
    Scalar eval_q_nu    (Scalar,Scalar);

    Scalar eval_exact_u  (Scalar,Scalar);
    Scalar eval_exact_v  (Scalar,Scalar);
    Scalar eval_exact_t  (Scalar,Scalar);
    Scalar eval_exact_rho(Scalar,Scalar);
    Scalar eval_exact_nu (Scalar,Scalar);
    Scalar eval_exact_p (Scalar,Scalar);
  };

  // ------------------------------------------------------
  // ---------- RANS: Spelart Alamaras (Free Shear) -------
  // ------------------------------------------------------
  template <typename Scalar>
  class   fans_sa_transient_free_shear : public manufactured_solution<Scalar>
  {

    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar u_0;
    Scalar u_x;
    Scalar u_y;
    Scalar v_0;
    Scalar v_x;
    Scalar v_y;
    Scalar rho_0;
    Scalar rho_x;
    Scalar rho_y;
    Scalar p_0;
    Scalar p_x;
    Scalar p_y;
    Scalar a_px;
    Scalar a_py;
    Scalar a_rhox;
    Scalar a_rhoy;
    Scalar a_ux;
    Scalar a_uy;
    Scalar a_vx;
    Scalar a_vy;
    Scalar mu;
    Scalar L;

    Scalar u_t;
    Scalar v_t;
    Scalar p_t;
    Scalar rho_t;
    Scalar a_ut;
    Scalar a_vt;
    Scalar a_pt;
    Scalar a_rhot;

    Scalar nu_sa_0;
    Scalar nu_sa_x;
    Scalar nu_sa_y;
    Scalar nu_sa_t;
    Scalar a_nusax;
    Scalar a_nusay;
    Scalar a_nusat;

    Scalar c_v1;
    Scalar c_b1;
    Scalar c_b2;
    Scalar c_w1;
    Scalar c_w2;
    Scalar c_w3;
    Scalar kappa;
    Scalar sigma;

    //Scalar cp;
    //Scalar cv;

    Scalar Gamma;
    Scalar R;

    Scalar Pr;
    Scalar Pr_t;

  public:
    fans_sa_transient_free_shear(); // constructor
    int init_var();

    // provide steady versions
    Scalar eval_q_rho_u (Scalar,Scalar);
    Scalar eval_q_rho_v (Scalar,Scalar);
    Scalar eval_q_rho_e (Scalar,Scalar);
    Scalar eval_q_rho   (Scalar,Scalar);
    Scalar eval_q_nu    (Scalar,Scalar);

    Scalar eval_exact_nu (Scalar,Scalar);

    // provide unsteady versions
    Scalar eval_q_rho_u (Scalar,Scalar,Scalar);
    Scalar eval_q_rho_v (Scalar,Scalar,Scalar);
    Scalar eval_q_rho_e (Scalar,Scalar,Scalar);
    Scalar eval_q_rho   (Scalar,Scalar,Scalar);
    Scalar eval_q_nu    (Scalar,Scalar,Scalar);

    Scalar eval_exact_u  (Scalar,Scalar);
    Scalar eval_exact_v  (Scalar,Scalar);
    Scalar eval_exact_p  (Scalar,Scalar);
    Scalar eval_exact_rho(Scalar,Scalar);
    Scalar eval_exact_nu (Scalar,Scalar,Scalar);
  };

  template <typename Scalar>
  class rans_sa : public manufactured_solution<Scalar>
  {

    Scalar cb1;
    Scalar cb2;
    Scalar cv1;
    Scalar cw2;
    Scalar cw3;
    Scalar sigma;
    Scalar kappa;
    Scalar re_tau;

    // adding for modified SA
    Scalar cv2;
    Scalar cv3;

    // parameters not for users
    Scalar etam;
    Scalar a1;
    Scalar b1;

  public:
    rans_sa(); // constructor
    int init_var();          // default problem values

    Scalar eval_q_u (Scalar); // velocity term
    Scalar eval_q_v (Scalar); // eddy viscosity term
    Scalar eval_exact_u(Scalar); // analytical
    Scalar eval_exact_v(Scalar);

    // member functions not exposed by API
    Scalar   u(Scalar);
    Scalar  du(Scalar);
    Scalar d2u();

    Scalar   nu(Scalar);
    Scalar  dnu(Scalar);
    Scalar d2nu(Scalar);

    Scalar production(Scalar);
    Scalar destruction(Scalar);
    Scalar transport(Scalar);

    Scalar cw1();
    Scalar chi(Scalar);
    Scalar fv1(Scalar);
    Scalar fv2(Scalar);
    Scalar   s(Scalar);
    Scalar  sb(Scalar);
    Scalar   r(Scalar);
    Scalar   g(Scalar);
    Scalar  fw(Scalar);
    Scalar dvt(Scalar);
    Scalar  vt(Scalar);

  };

  // ------------------------------------------------------
  // ---------- Integrated Radiative Intensity   ----------
  // ------------------------------------------------------
  template <typename Scalar>
  class radiation_integrated_intensity : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar no_gauss;
    std::vector<Scalar> vec_mean;
    std::vector<Scalar> vec_amp;
    std::vector<Scalar> vec_stdev;

  public:
    radiation_integrated_intensity(); // constructor
    int init_var();
    int check_vec();  // sanity check

    Scalar eval_q_u (Scalar);
    Scalar eval_exact_u(Scalar);
    Scalar phi(Scalar);

  };

  // ------------------------------------------------------
  // ---------- Laplace's Equation               ----------
  // ------------------------------------------------------
  template <typename Scalar>
  class laplace_2d : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar Ly;
    Scalar Lx;

  public:
    laplace_2d(); // constructor
    int init_var();

    Scalar eval_q_f(Scalar,Scalar);
    Scalar eval_exact_phi(Scalar,Scalar);

  };

  // ------------------------------------------------------
  // ------------- compressible navier stokes  ------------
  // ------------------------------------------------------

  template <typename Scalar>
  class navierstokes_2d_compressible : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar R;     // (or is this the ideal gas constant?) ratio of specific heat capacities, defined in constructor
    Scalar k;     // Boltzmanns constant

    Scalar u_0;
    Scalar u_x;
    Scalar u_y;
    Scalar v_0;
    Scalar v_x;
    Scalar v_y;
    Scalar rho_0;
    Scalar rho_x;
    Scalar rho_y;
    Scalar p_0;
    Scalar p_x;
    Scalar p_y;
    Scalar a_px;
    Scalar a_py;
    Scalar a_rhox;
    Scalar a_rhoy;
    Scalar a_ux;
    Scalar a_uy;
    Scalar a_vx;
    Scalar a_vy;
    Scalar Gamma;
    Scalar mu;
    Scalar L;

  public:
    navierstokes_2d_compressible(); // constructor
    int init_var();        // default problem values

    Scalar eval_q_rho_u (Scalar,Scalar);
    Scalar eval_q_rho_v (Scalar,Scalar);
    Scalar eval_q_rho_e (Scalar,Scalar);
    Scalar eval_q_rho   (Scalar,Scalar);

    Scalar eval_exact_u  (Scalar,Scalar); // analytical
    Scalar eval_exact_v  (Scalar,Scalar);
    Scalar eval_exact_p  (Scalar,Scalar);
    Scalar eval_exact_rho(Scalar,Scalar);

    Scalar eval_g_u(Scalar,Scalar,int);   // gradient of source term
    Scalar eval_g_v(Scalar,Scalar,int);
    Scalar eval_g_p(Scalar,Scalar,int);

    Scalar eval_g_rho(Scalar,Scalar,int);


  };

  template <typename Scalar>
  class navierstokes_3d_compressible : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar R;     // (or is this the ideal gas constant?) ratio of specific heat capacities, defined in constructor
    Scalar k;     // Boltzmanns constant

    Scalar u_0;
    Scalar u_x;
    Scalar u_y;
    Scalar u_z;
    Scalar v_0;
    Scalar v_x;
    Scalar v_y;
    Scalar v_z;
    Scalar w_0;
    Scalar w_x;
    Scalar w_y;
    Scalar w_z;
    Scalar rho_0;
    Scalar rho_x;
    Scalar rho_y;
    Scalar rho_z;
    Scalar p_0;
    Scalar p_x;
    Scalar p_y;
    Scalar p_z;
    Scalar a_px;
    Scalar a_py;
    Scalar a_pz;
    Scalar a_rhox;
    Scalar a_rhoy;
    Scalar a_rhoz;
    Scalar a_ux;
    Scalar a_uy;
    Scalar a_uz;
    Scalar a_vx;
    Scalar a_vy;
    Scalar a_vz;
    Scalar a_wx;
    Scalar a_wy;
    Scalar a_wz;
    Scalar mu;
    Scalar Gamma;
    Scalar L;

  public:
    navierstokes_3d_compressible(); //constructor
    int init_var();        // default problem values

    Scalar eval_q_rho_u  (Scalar,Scalar,Scalar); // source terms
    Scalar eval_q_rho_v  (Scalar,Scalar,Scalar);
    Scalar eval_q_rho_w  (Scalar,Scalar,Scalar);
    Scalar eval_q_rho_e  (Scalar,Scalar,Scalar);
    Scalar eval_q_rho    (Scalar,Scalar,Scalar);

    Scalar eval_exact_u  (Scalar,Scalar,Scalar); // analytical
    Scalar eval_exact_v  (Scalar,Scalar,Scalar);
    Scalar eval_exact_w  (Scalar,Scalar,Scalar);
    Scalar eval_exact_q  (Scalar,Scalar,Scalar);
    Scalar eval_exact_p  (Scalar,Scalar,Scalar);
    Scalar eval_exact_rho(Scalar,Scalar,Scalar);

    Scalar eval_g_u(Scalar,Scalar,Scalar,int);   // gradient of source term
    Scalar eval_g_v(Scalar,Scalar,Scalar,int);
    Scalar eval_g_w(Scalar,Scalar,Scalar,int);

    Scalar eval_g_p  (Scalar,Scalar,Scalar,int);
    Scalar eval_g_rho(Scalar,Scalar,Scalar,int);

  }; // done with navier stokes 3d class

// External nsctpl::manufactured_solution implementation pulled into MASA namespace
#include "nsctpl_fwd.hpp"

  template <typename Scalar>
  class navierstokes_4d_compressible_powerlaw
    : public  MASA::manufactured_solution<Scalar>,
      private MASA::nsctpl::manufactured_solution<Scalar, /*One-indexed gradients*/1>
  {

  public:

    // constructor
    navierstokes_4d_compressible_powerlaw();

    // default parameters set by nsctpl::zero and nsctpl::isothermal_channel
    friend void MASA::nsctpl::zero<>(navierstokes_4d_compressible_powerlaw<Scalar>&);
    friend void MASA::nsctpl::isothermal_channel<>(navierstokes_4d_compressible_powerlaw<Scalar>&);
    int init_var();

    // solution
    Scalar eval_exact_rho(Scalar x, Scalar y, Scalar z, Scalar t) { return this->rho(x,y,z,t); }
    Scalar eval_exact_u  (Scalar x, Scalar y, Scalar z, Scalar t) { return this->u  (x,y,z,t); }
    Scalar eval_exact_v  (Scalar x, Scalar y, Scalar z, Scalar t) { return this->v  (x,y,z,t); }
    Scalar eval_exact_w  (Scalar x, Scalar y, Scalar z, Scalar t) { return this->w  (x,y,z,t); }
    Scalar eval_exact_t  (Scalar x, Scalar y, Scalar z, Scalar t) { return this->T  (x,y,z,t); }
    Scalar eval_exact_p  (Scalar x, Scalar y, Scalar z, Scalar t) { return this->p  (x,y,z,t); }

    // source terms
    Scalar eval_q_rho  (Scalar x, Scalar y, Scalar z, Scalar t) { return this->Q_rho (x,y,z,t); }
    Scalar eval_q_rho_u(Scalar x, Scalar y, Scalar z, Scalar t) { return this->Q_rhou(x,y,z,t); }
    Scalar eval_q_rho_v(Scalar x, Scalar y, Scalar z, Scalar t) { return this->Q_rhov(x,y,z,t); }
    Scalar eval_q_rho_w(Scalar x, Scalar y, Scalar z, Scalar t) { return this->Q_rhow(x,y,z,t); }
    Scalar eval_q_rho_e(Scalar x, Scalar y, Scalar z, Scalar t) { return this->Q_rhoe(x,y,z,t); }

    // gradient of solution
    Scalar eval_g_rho(Scalar x, Scalar y, Scalar z, Scalar t, int i) { return this->grad_rho(x,y,z,t,i); }
    Scalar eval_g_u  (Scalar x, Scalar y, Scalar z, Scalar t, int i) { return this->grad_u  (x,y,z,t,i); }
    Scalar eval_g_v  (Scalar x, Scalar y, Scalar z, Scalar t, int i) { return this->grad_v  (x,y,z,t,i); }
    Scalar eval_g_w  (Scalar x, Scalar y, Scalar z, Scalar t, int i) { return this->grad_w  (x,y,z,t,i); }
    Scalar eval_g_t  (Scalar x, Scalar y, Scalar z, Scalar t, int i) { return this->grad_T  (x,y,z,t,i); }
    Scalar eval_g_p  (Scalar x, Scalar y, Scalar z, Scalar t, int i) { return this->grad_p  (x,y,z,t,i); }

  }; // done with navierstokes_4d_compressible_powerlaw class

  template <typename Scalar>
  class navierstokes_ablation_1d_steady : public manufactured_solution<Scalar>
  {
    using manufactured_solution<Scalar>::pi;
    using manufactured_solution<Scalar>::PI;

    Scalar R;
    Scalar k;

    Scalar u_0;
    Scalar u_x;
    Scalar a_ux;
    Scalar Gamma;
    Scalar mu;
    Scalar L;

    Scalar T_0;
    Scalar T_x;
    Scalar a_Tx;

    Scalar W_C;
    Scalar W_C3;

    Scalar rho_N_0;
    Scalar rho_N_x;
    Scalar rho_N2_0;
    Scalar rho_N2_x;
    Scalar a_rho_N2_x;

    Scalar rho_C3_0;
    Scalar rho_C3_x;
    Scalar a_rho_C3_x;

    Scalar rho_C_0;
    Scalar rho_C_x;
    Scalar a_rho_C_x;

    Scalar rho_an_C3;

    Scalar k_B;
    Scalar beta_C3;
    Scalar A_C3Enc;
    Scalar D_C;
    Scalar D_C3;
    Scalar m_C3;
    Scalar E_aC3nc;

    Scalar qr;
    Scalar alpha;
    Scalar a_rho_N_x;
    Scalar sigma;
    Scalar epsilon;

  public:
    navierstokes_ablation_1d_steady(); // constructor
    int init_var();        // default problem values

    // source terms
    // 2d -> 1d
    // note: using 'real' here instead of scalar to avoid compilation difficulties.
    template <typename Real> Real eval_q_rho_u     (Real x,Real y){return eval_q_rho_u(x);}
    template <typename Real> Real eval_q_e         (Real x,Real y,Real (*in_func)(Real)) {return eval_q_e(x,in_func);}
    template <typename Real> Real eval_q_C         (Real x,Real y){return eval_q_C(x);}
    template <typename Real> Real eval_q_C3        (Real x,Real y){return eval_q_C3(x);}
    template <typename Real> Real eval_q_rho_C     (Real x,Real y){return eval_q_rho_C(x);}
    template <typename Real> Real eval_q_rho_C3    (Real x,Real y){return eval_q_rho_C3(x);}
    template <typename Real> Real eval_q_u_boundary(Real x,Real y){return eval_q_u_boundary(x);}

    // 3d-> 1d
    template <typename Real> Real eval_q_rho_u     (Real x,Real y,Real z){return eval_q_rho_u(x);}
    template <typename Real> Real eval_q_e         (Real x,Real y,Real z,Real (*in_func)(Real)) {return eval_q_e(x,in_func);}
    template <typename Real> Real eval_q_C         (Real x,Real y,Real z){return eval_q_C(x);}
    template <typename Real> Real eval_q_C3        (Real x,Real y,Real z){return eval_q_C3(x);}
    template <typename Real> Real eval_q_rho_C     (Real x,Real y,Real z){return eval_q_rho_C(x);}
    template <typename Real> Real eval_q_rho_C3    (Real x,Real y,Real z){return eval_q_rho_C3(x);}
    template <typename Real> Real eval_q_u_boundary(Real x,Real y,Real z){return eval_q_u_boundary(x);}

    // manufactured analytical solutions
    // 3d-> 1d
    template <typename Real> Real eval_exact_u     (Real x,Real y,Real z) {return eval_exact_u(x);}
    template <typename Real> Real eval_exact_t     (Real x,Real y,Real z) {return eval_exact_t(x);}
    template <typename Real> Real eval_exact_rho   (Real x,Real y,Real z) {return eval_exact_rho(x);}
    template <typename Real> Real eval_exact_rho_C (Real x,Real y,Real z) {return eval_exact_rho_C(x);}
    template <typename Real> Real eval_exact_rho_C3(Real x,Real y,Real z) {return eval_exact_rho_C3(x);}

    // 2d-> 1d
    template <typename Real> Real eval_exact_u     (Real x,Real y) {return eval_exact_u(x);}
    template <typename Real> Real eval_exact_t     (Real x,Real y) {return eval_exact_t(x);}
    template <typename Real> Real eval_exact_rho   (Real x,Real y) {return eval_exact_rho(x);}
    template <typename Real> Real eval_exact_rho_C (Real x,Real y) {return eval_exact_rho_C(x);}
    template <typename Real> Real eval_exact_rho_C3(Real x,Real y) {return eval_exact_rho_C3(x);}

    // 1d (Actual) Source terms
    Scalar eval_q_rho_u (Scalar);
    Scalar eval_q_rho_e (Scalar);
    Scalar eval_q_e     (Scalar, Scalar (*)(Scalar));
    Scalar eval_q_C     (Scalar);
    Scalar eval_q_C3    (Scalar);
    Scalar eval_q_rho_C (Scalar);
    Scalar eval_q_rho_C3(Scalar);
    Scalar eval_q_u_boundary (Scalar);

    Scalar eval_exact_u     (Scalar);
    Scalar eval_exact_t     (Scalar);
    Scalar eval_exact_rho   (Scalar);
    Scalar eval_exact_rho_C (Scalar);
    Scalar eval_exact_rho_C3(Scalar);
  };


} // end MASA namespace


// DO NOT EDIT THIS LINE OR ANY BELOW IT
// automatically generated MMS classes:
  
// ------------------------------------------------------
// --------------- burgers_equation 
// ------------------------------------------------------
namespace MASA{
template <typename Scalar>
class burgers_equation : public manufactured_solution<Scalar>
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

private:
  Scalar nu;
  Scalar u_0;
  Scalar u_x;
  Scalar u_y;
  Scalar u_t;
  Scalar v_0;
  Scalar v_x;
  Scalar v_y;
  Scalar v_t;
  Scalar a_ux;
  Scalar a_uy;
  Scalar a_ut;
  Scalar a_vx;
  Scalar a_vy;
  Scalar a_vt;
  Scalar L;
  Scalar Lt;

public:
  burgers_equation();
  int init_var();

  Scalar eval_q_v_transient_viscous (Scalar,Scalar,Scalar);
  Scalar eval_q_u_transient_viscous (Scalar,Scalar,Scalar);

  //Scalar eval_q_v_transient_inviscid (Scalar,Scalar,Scalar);
  //Scalar eval_q_u_transient_inviscid (Scalar,Scalar,Scalar);
  Scalar eval_q_v (Scalar,Scalar,Scalar);
  Scalar eval_q_u (Scalar,Scalar,Scalar);

  Scalar eval_q_v_steady_viscous (Scalar,Scalar);
  Scalar eval_q_v_steady_inviscid (Scalar,Scalar);
  Scalar eval_q_u_steady_viscous (Scalar,Scalar);
  Scalar eval_q_u_steady_inviscid (Scalar,Scalar);
  Scalar eval_exact_u(Scalar,Scalar);
  Scalar eval_exact_v(Scalar,Scalar);
  Scalar eval_exact_u(Scalar,Scalar,Scalar);
  Scalar eval_exact_v(Scalar,Scalar,Scalar);
};}


// ------------------------------------------------------
// --------------- euler_transient_2d 
// ------------------------------------------------------
namespace MASA{
template <typename Scalar>
class euler_transient_2d : public manufactured_solution<Scalar>
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

private:
  Scalar rho_0;
  Scalar rho_x;
  Scalar rho_y;
  Scalar rho_t;
  Scalar u_0;
  Scalar u_x;
  Scalar u_y;
  Scalar u_t;
  Scalar v_0;
  Scalar v_x;
  Scalar v_y;
  Scalar v_t;
  Scalar L;
  Scalar p_0;
  Scalar p_x;
  Scalar p_y;
  Scalar p_t;
  Scalar a_vx;
  Scalar a_vy;
  Scalar a_vt;
  Scalar a_ux;
  Scalar a_uy;
  Scalar a_ut;
  Scalar a_rhox;
  Scalar a_rhoy;
  Scalar a_rhot;
  Scalar a_py;
  Scalar a_px;
  Scalar a_pt;
  Scalar Gamma;

public:
  euler_transient_2d();
  int init_var();
  Scalar eval_q_e (Scalar,Scalar,Scalar);
  Scalar eval_q_u (Scalar,Scalar,Scalar);
  Scalar eval_q_v (Scalar,Scalar,Scalar);
  Scalar eval_q_rho (Scalar,Scalar,Scalar);
  Scalar eval_exact_rho(Scalar,Scalar,Scalar);
  Scalar eval_exact_p(Scalar,Scalar,Scalar);
  Scalar eval_exact_u(Scalar,Scalar,Scalar);
  Scalar eval_exact_v(Scalar,Scalar,Scalar);
};}


// ------------------------------------------------------
// --------------- euler_transient_3d 
// ------------------------------------------------------
namespace MASA{
template <typename Scalar>
class euler_transient_3d : public manufactured_solution<Scalar>
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

private:
  Scalar rho_0;
  Scalar rho_x;
  Scalar rho_y;
  Scalar rho_z;
  Scalar rho_t;
  Scalar u_0;
  Scalar u_x;
  Scalar u_y;
  Scalar u_z;
  Scalar u_t;
  Scalar v_0;
  Scalar v_x;
  Scalar v_y;
  Scalar v_z;
  Scalar v_t;
  Scalar w_0;
  Scalar w_x;
  Scalar w_y;
  Scalar w_z;
  Scalar w_t;
  Scalar L;
  Scalar p_0;
  Scalar p_x;
  Scalar p_y;
  Scalar p_z;
  Scalar p_t;
  Scalar a_vx;
  Scalar a_vy;
  Scalar a_vz;
  Scalar a_vt;
  Scalar a_wx;
  Scalar a_wy;
  Scalar a_wz;
  Scalar a_wt;
  Scalar a_ux;
  Scalar a_uy;
  Scalar a_uz;
  Scalar a_ut;
  Scalar a_rhox;
  Scalar a_rhoy;
  Scalar a_rhoz;
  Scalar a_rhot;
  Scalar a_py;
  Scalar a_px;
  Scalar a_pz;
  Scalar a_pt;
  Scalar Gamma;

public:
  euler_transient_3d();
  int init_var();
  Scalar eval_q_e (Scalar,Scalar,Scalar,Scalar);
  Scalar eval_q_u (Scalar,Scalar,Scalar,Scalar);
  Scalar eval_q_v (Scalar,Scalar,Scalar,Scalar);
  Scalar eval_q_w (Scalar,Scalar,Scalar,Scalar);
  Scalar eval_q_rho (Scalar,Scalar,Scalar,Scalar);
  Scalar eval_exact_rho(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_exact_p(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_exact_u(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_exact_v(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_exact_w(Scalar,Scalar,Scalar,Scalar);
};}


// ------------------------------------------------------
// --------------- axi_euler_transient 
// ------------------------------------------------------
namespace MASA{
template <typename Scalar>
class axi_euler_transient : public manufactured_solution<Scalar>
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

private:
  Scalar rho_0;
  Scalar rho_r;
  Scalar rho_z;
  Scalar rho_t;
  Scalar a_rhoz;
  Scalar a_rhor;
  Scalar a_rhot;
  Scalar u_0;
  Scalar u_r;
  Scalar u_z;
  Scalar u_t;
  Scalar a_uz;
  Scalar a_ur;
  Scalar a_ut;
  Scalar w_0;
  Scalar w_r;
  Scalar w_z;
  Scalar w_t;
  Scalar a_wz;
  Scalar a_wr;
  Scalar a_wt;
  Scalar p_0;
  Scalar p_r;
  Scalar p_z;
  Scalar p_t;
  Scalar a_pz;
  Scalar a_pr;
  Scalar a_pt;
  Scalar L;
  Scalar Gamma;

public:
  axi_euler_transient();
  int init_var();
  Scalar eval_q_e (Scalar,Scalar,Scalar);
  Scalar eval_q_u (Scalar,Scalar,Scalar);
  Scalar eval_q_w (Scalar,Scalar,Scalar);
  Scalar eval_q_rho (Scalar,Scalar,Scalar);
  Scalar eval_exact_rho(Scalar,Scalar,Scalar);
  Scalar eval_exact_p(Scalar,Scalar,Scalar);
  Scalar eval_exact_u(Scalar,Scalar,Scalar);
  Scalar eval_exact_w(Scalar,Scalar,Scalar);
};}


// ------------------------------------------------------
// --------------- axi_cns_transient 
// ------------------------------------------------------
namespace MASA{
template <typename Scalar>
class axi_cns_transient : public manufactured_solution<Scalar>
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

private:
  Scalar rho_0;
  Scalar rho_r;
  Scalar rho_z;
  Scalar rho_t;
  Scalar a_rhoz;
  Scalar a_rhor;
  Scalar a_rhot;
  Scalar u_0;
  Scalar u_r;
  Scalar u_z;
  Scalar u_t;
  Scalar a_uz;
  Scalar a_ur;
  Scalar a_ut;
  Scalar w_0;
  Scalar w_r;
  Scalar w_z;
  Scalar w_t;
  Scalar a_wz;
  Scalar a_wr;
  Scalar a_wt;
  Scalar p_0;
  Scalar p_r;
  Scalar p_z;
  Scalar p_t;
  Scalar a_pz;
  Scalar a_pr;
  Scalar a_pt;
  Scalar L;
  Scalar Gamma;
  Scalar mu;
  Scalar k;
  Scalar R;

public:
  axi_cns_transient();
  int init_var();
  Scalar eval_q_e (Scalar,Scalar,Scalar);
  Scalar eval_q_u (Scalar,Scalar,Scalar);
  Scalar eval_q_w (Scalar,Scalar,Scalar);
  Scalar eval_q_rho (Scalar,Scalar,Scalar);
  Scalar eval_exact_rho(Scalar,Scalar,Scalar);
  Scalar eval_exact_p(Scalar,Scalar,Scalar);
  Scalar eval_exact_u(Scalar,Scalar,Scalar);
  Scalar eval_exact_w(Scalar,Scalar,Scalar);
};}


// MetaPhysicL-powered manufactured solutions:
#ifdef HAVE_METAPHYSICL

// ------------------------------------------------------
// --------------- ad_cns_2d_crossterms 
// ------------------------------------------------------
namespace MASA{
template <typename Scalar>
class ad_cns_2d_crossterms : public manufactured_solution<Scalar>
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

private:
  Scalar R;
  Scalar k;
  Scalar u_0;
  Scalar u_x;
  Scalar u_y;
  Scalar v_0;
  Scalar v_x;
  Scalar v_y;
  Scalar rho_0;
  Scalar rho_x;
  Scalar rho_y;
  Scalar p_0;
  Scalar p_x;
  Scalar p_y;
  Scalar a_px;
  Scalar a_py;
  Scalar a_rhox;
  Scalar a_rhoy;
  Scalar a_ux;
  Scalar a_uy;
  Scalar a_vx;
  Scalar a_vy;
  Scalar Gamma;
  Scalar mu;
  Scalar L;

public:
  ad_cns_2d_crossterms();
  int init_var();
  Scalar eval_q_u(Scalar,Scalar) const;
  Scalar eval_q_v(Scalar,Scalar) const;
  Scalar eval_q_e(Scalar,Scalar) const;
  Scalar eval_q_rho(Scalar,Scalar) const;
  Scalar eval_exact_u(Scalar,Scalar);
  Scalar eval_exact_v(Scalar,Scalar);
  Scalar eval_exact_p(Scalar,Scalar);
  Scalar eval_exact_rho(Scalar,Scalar);
};}


// ------------------------------------------------------
// --------------- ad_cns_3d_crossterms 
// ------------------------------------------------------
namespace MASA{
template <typename Scalar>
class ad_cns_3d_crossterms : public manufactured_solution<Scalar>
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

private:
  Scalar R;
  Scalar k;
  Scalar u_0;
  Scalar u_x;
  Scalar u_y;
  Scalar u_z;
  Scalar v_0;
  Scalar v_x;
  Scalar v_y;
  Scalar v_z;
  Scalar w_0;
  Scalar w_x;
  Scalar w_y;
  Scalar w_z;
  Scalar rho_0;
  Scalar rho_x;
  Scalar rho_y;
  Scalar rho_z;
  Scalar p_0;
  Scalar p_x;
  Scalar p_y;
  Scalar p_z;
  Scalar a_px;
  Scalar a_py;
  Scalar a_pz;
  Scalar a_rhox;
  Scalar a_rhoy;
  Scalar a_rhoz;
  Scalar a_ux;
  Scalar a_uy;
  Scalar a_uz;
  Scalar a_vx;
  Scalar a_vy;
  Scalar a_vz;
  Scalar a_wx;
  Scalar a_wy;
  Scalar a_wz;
  Scalar Gamma;
  Scalar mu;
  Scalar L;

public:
  ad_cns_3d_crossterms();
  int init_var();
  Scalar eval_q_u(Scalar,Scalar,Scalar) const;
  Scalar eval_q_v(Scalar,Scalar,Scalar) const;
  Scalar eval_q_w(Scalar,Scalar,Scalar) const;
  Scalar eval_q_e(Scalar,Scalar,Scalar) const;
  Scalar eval_q_rho(Scalar,Scalar,Scalar) const;
  Scalar eval_exact_u(Scalar,Scalar,Scalar);
  Scalar eval_exact_v(Scalar,Scalar,Scalar);
  Scalar eval_exact_w(Scalar,Scalar,Scalar);
  Scalar eval_exact_p(Scalar,Scalar,Scalar);
  Scalar eval_exact_rho(Scalar,Scalar,Scalar);
};}


// ------------------------------------------------------
// --------------- convdiff_steady_nosource_1d 
// ------------------------------------------------------
namespace MASA{
template <typename Scalar>
class convdiff_steady_nosource_1d : public manufactured_solution<Scalar>
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

private:
  Scalar a_ux;
  Scalar a_cx;
  Scalar nu;
  Scalar L;

public:
  convdiff_steady_nosource_1d();
  int init_var();
  Scalar eval_q_c(Scalar,Scalar) const;
  Scalar eval_exact_c(Scalar);
  Scalar eval_exact_u(Scalar);
};}


// ------------------------------------------------------
// --------------- navierstokes_3d_incompressible 
// ------------------------------------------------------
namespace MASA{
template <typename Scalar>
class navierstokes_3d_incompressible : public manufactured_solution<Scalar>
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

private:
  Scalar a;
  Scalar b;
  Scalar c;
  Scalar d;
  Scalar beta;
  Scalar gamma;
  Scalar nu;
  Scalar kx;
  Scalar kz;

public:
  navierstokes_3d_incompressible();
  int init_var();
  Scalar eval_q_u(Scalar,Scalar,Scalar);
  Scalar eval_q_v(Scalar,Scalar,Scalar);
  Scalar eval_q_w(Scalar,Scalar,Scalar);
  // Scalar helper_f(Scalar);
  // Scalar helper_g(Scalar);
  // Scalar helper_gt(Scalar);
  // Scalar helper_h(Scalar);
  Scalar eval_exact_u(Scalar,Scalar,Scalar);
  Scalar eval_exact_v(Scalar,Scalar,Scalar);
  Scalar eval_exact_w(Scalar,Scalar,Scalar);
  Scalar eval_exact_p(Scalar,Scalar,Scalar);
};}


// ------------------------------------------------------
// --------------- navierstokes_3d_incompressible_homogeneous 
// ------------------------------------------------------
namespace MASA{
template <typename Scalar>
class navierstokes_3d_incompressible_homogeneous : public manufactured_solution<Scalar>
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

private:
  Scalar a;
  Scalar b;
  Scalar c;
  Scalar d;
  Scalar beta;
  Scalar gamma;
  Scalar delta;
  Scalar nu;
  Scalar kx;
  Scalar ky;
  Scalar kz;

public:
  navierstokes_3d_incompressible_homogeneous();
  int init_var();
  Scalar eval_q_u(Scalar,Scalar,Scalar);
  Scalar eval_q_v(Scalar,Scalar,Scalar);
  Scalar eval_q_w(Scalar,Scalar,Scalar);
  Scalar eval_exact_u(Scalar,Scalar,Scalar);
  Scalar eval_exact_v(Scalar,Scalar,Scalar);
  Scalar eval_exact_w(Scalar,Scalar,Scalar);
  Scalar eval_exact_p(Scalar,Scalar,Scalar);
};}


// ------------------------------------------------------
// --------------- navierstokes_3d_transient_sutherland 
// ------------------------------------------------------
namespace MASA{
template <typename Scalar>
class navierstokes_3d_transient_sutherland : public manufactured_solution<Scalar>
{
  using manufactured_solution<Scalar>::pi;
  using manufactured_solution<Scalar>::PI;

private:
  Scalar L;
  Scalar Lt;
  Scalar rho_0;
  Scalar rho_x;
  Scalar a_rhox;
  Scalar rho_y;
  Scalar a_rhoy;
  Scalar rho_z;
  Scalar a_rhoz;
  Scalar rho_t;
  Scalar a_rhot;
  Scalar p_0;
  Scalar p_x;
  Scalar a_px;
  Scalar p_y;
  Scalar a_py;
  Scalar p_z;
  Scalar a_pz;
  Scalar p_t;
  Scalar a_pt;
  Scalar u_0;
  Scalar u_x;
  Scalar a_ux;
  Scalar u_y;
  Scalar a_uy;
  Scalar u_z;
  Scalar a_uz;
  Scalar u_t;
  Scalar a_ut;
  Scalar v_0;
  Scalar v_x;
  Scalar a_vx;
  Scalar v_y;
  Scalar a_vy;
  Scalar v_z;
  Scalar a_vz;
  Scalar v_t;
  Scalar a_vt;
  Scalar w_0;
  Scalar w_x;
  Scalar a_wx;
  Scalar w_y;
  Scalar a_wy;
  Scalar w_z;
  Scalar a_wz;
  Scalar w_t;
  Scalar a_wt;
  Scalar B_mu;
  Scalar A_mu;
  Scalar Gamma;
  Scalar R;
  Scalar Pr;

public:
  navierstokes_3d_transient_sutherland();
  int init_var();
  Scalar eval_q_e(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_q_rho(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_q_u(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_q_v(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_q_w(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_exact_rho(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_exact_p(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_exact_u(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_exact_v(Scalar,Scalar,Scalar,Scalar);
  Scalar eval_exact_w(Scalar,Scalar,Scalar,Scalar);
};}


#endif // HAVE_METAPHYSICL

// --l33t--
