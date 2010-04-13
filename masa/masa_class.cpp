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

using <masa.h>

namespace MASA {
  
  manufactured_solution::manufactured_solution()
  {
    
    mms_name="uninit"; // create bad name
    
  }
  
  manufactured_solution::~manufactured_solution()
  {
    
    //need to free everything here
    
  }
  
  manufactured_solution::return_name()
  {
    
    cout << endl << mms_name << endl;
    
  }
  
  manufactured_solution::set_name(string input_string)
  {
    
    mms_name=input_string;

  }

};
