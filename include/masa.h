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

// includes

#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>
#include <iostream>
#include <map>

using namespace std;

namespace MASA 
{

  class manufactured_solution
  {
  private: 
    string mms_name; // the name of the manufactured solution  
    double Tan;      // analytical solution 
    double Q;        // source term
    gradT;           // gradient 
    
  public: 
    manufactured_solution ();  // constructor
    ~manufactured_solution();  // destructor
    void return_name      ();  // method: returns name
    void set_name   (string);  // method: sets name
    
  }; // done with MMS class

  // masa map functions here
  // probably want to hide this from the user eventually

  int masa_map_solution (string, string);
  int masa_map_temporal (string, string);
  int masa_map_coeff    (string, string); 
  int masa_map          (string, string);
  int masa_map_dimension(string, string);

  // masa_shell
  //void masa_shell_choose_solution();

} //end MASA namespace
