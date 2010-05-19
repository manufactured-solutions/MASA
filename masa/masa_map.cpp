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

// these routines take a string, parse it
// and returns what we 'really' wanted the user to input

// ideas:
// need many more sanity/error checks
// how am I going to automate this all?

#include <masa_internal.h>

using namespace MASA;

// -------------------------------------------------------
// these are just a few utility functions
// -------------------------------------------------------

// all this does is convert a string to entirely lower case
// should simplify our user stuff
string uptolow(string str) 
{
  for (int i=0;i<strlen(str.c_str());i++) 
    if (str[i] >= 0x41 && str[i] <= 0x5A) 
      str[i] = str[i] + 0x20;
  return str;
}

// remove dashes from strings
string remove_line(string string1)
{
  int position = string1.find( "-" ); // find first space
  while ( position != string::npos )
    {
      string1.replace( position, 1, "" );
      position = string1.find( "-", position + 1 );
    }
  return string1;
}

string remove_whitespace(string string1)
{
  int position = string1.find( " " ); // find first space
  while ( position != string::npos )
    {
      string1.replace( position, 1, "" );
      position = string1.find( " ", position + 1 );
    }
  return string1;
}

// -------------------------------------------------------
// various parameter mappings
// -------------------------------------------------------

// user chooses 1d,2d, or 3d 
int MASA::masa_map_dimension(string input_string, string return_string)
{
  multimap<string, string> names;
  string temp_string;
  int ret_val=0;  // exits with zero upon failure, 1 for success

  // convert to lower case, remove dashes: -
  temp_string=uptolow(input_string);
  temp_string=remove_line(temp_string);
  temp_string=remove_whitespace(temp_string);

  // map 1-d
  names.insert(pair<string, string>("1",               "1d"));
  names.insert(pair<string, string>("1d",              "1d"));
  names.insert(pair<string, string>("onedimensional",  "1d"));
  names.insert(pair<string, string>("one",             "1d"));
  names.insert(pair<string, string>("oned",            "1d"));

  // map 2-d
  names.insert(pair<string, string>("2",               "2d"));
  names.insert(pair<string, string>("2d",              "2d"));
  names.insert(pair<string, string>("twodimensional",  "2d"));
  names.insert(pair<string, string>("two",             "2d"));
  names.insert(pair<string, string>("twod",            "2d"));

  // map 3-d
  names.insert(pair<string, string>("3",               "3d"));
  names.insert(pair<string, string>("3d",              "3d"));
  names.insert(pair<string, string>("threedimensional","3d"));
  names.insert(pair<string, string>("three",           "3d"));
  names.insert(pair<string, string>("threed",          "3d"));

  multimap<string, string>::iterator p;

  p = names.find(temp_string);
  if(p != names.end()) 
    { // found a name
      do {
	cout << temp_string << ", " << p->second;
	cout << endl;
	p++;
      } while (p != names.upper_bound(temp_string));
      ret_val=1;
    }
  else
    {
      cout << "Name not found.\n";
      ret_val = 0;
    }
      
  return ret_val;

}// end masa_map_dimension


// user chooses type of solution
int MASA::masa_map_solution(string input_string, string return_string)
{
  multimap<string, string> names;
  string temp_string;
  int ret_val=0;  // exits with zero upon failure, 1 for success

  // convert to lower case
  temp_string=uptolow(input_string);
  temp_string=remove_line(temp_string);
  temp_string=remove_whitespace(temp_string);

  // map to heat equation
  names.insert(pair<string, string>("he",               "heatequation"));
  names.insert(pair<string, string>("heeq",             "heatequation"));
  names.insert(pair<string, string>("heatequation",     "heatequation"));
  names.insert(pair<string, string>("heateq",           "heatequation"));
  names.insert(pair<string, string>("heat",             "heatequation"));

  // map to navier-stokes
  names.insert(pair<string, string>("ns",               "navier-stokes"));
  names.insert(pair<string, string>("navierstokes",     "navier-stokes"));
  names.insert(pair<string, string>("navierstokeseq",   "navier-stokes"));
  names.insert(pair<string, string>("navierstokesequation",     "navier-stokes"));

  // map euler equations
  names.insert(pair<string, string>("euler",            "euler"));
  names.insert(pair<string, string>("eulereq",          "euler"));
  names.insert(pair<string, string>("eulerequation",    "euler"));

  multimap<string, string>::iterator p;

  p = names.find(temp_string);
  if(p != names.end()) { // found a name
    do {
      cout << temp_string << ", " << p->second;
      cout << endl;
      p++;
    } while (p != names.upper_bound(temp_string));
    ret_val=1;
  }
  else{
    cout << "Name not found.\n";
    ret_val = 0;
  }

  return ret_val;

}// end masa_map_solution


// user chooses if solution has constant or variable coefficients
int MASA::masa_map_coeff(string input_string, string return_string)
{
  multimap<string, string> names;
  string temp_string;
  int ret_val=0;  // exits with zero upon failure, 1 for success

  // convert to lower case
  temp_string=uptolow(input_string);
  temp_string=remove_line(temp_string);
  temp_string=remove_whitespace(temp_string);

  // map to constant coefficients
  names.insert(pair<string, string>("const",                  "constant-coefficients"));
  names.insert(pair<string, string>("constant",               "constant-coefficients"));
  names.insert(pair<string, string>("constantcoeff",          "constant-coefficients"));
  names.insert(pair<string, string>("constantcoefficients",   "constant-coefficients"));
  names.insert(pair<string, string>("constantcoefficient",    "constant-coefficients"));
  names.insert(pair<string, string>("constcoefficient",       "constant-coefficients"));
  names.insert(pair<string, string>("constcoeff",             "constant-coefficients"));

  // map to variable coefficients
  names.insert(pair<string, string>("var",                    "variable-coefficients"));
  names.insert(pair<string, string>("variable",               "variable-coefficients"));
  names.insert(pair<string, string>("varcoeff",               "variable-coefficients"));
  names.insert(pair<string, string>("varcoefficient",         "variable-coefficients"));
  names.insert(pair<string, string>("varcoefficients",        "variable-coefficients"));
  names.insert(pair<string, string>("variablecoeff",          "variable-coefficients"));
  names.insert(pair<string, string>("variablecoefficient",    "variable-coefficients"));
  names.insert(pair<string, string>("variablecoefficients",   "variable-coefficients"));

  multimap<string, string>::iterator p;

  p = names.find(temp_string);
  if(p != names.end()) { // found a name
    do {
      cout << temp_string << ", " << p->second;
      cout << endl;
      p++;
    } while (p != names.upper_bound(temp_string));
    ret_val=1;
  }
  else{
    cout << "Name not found.\n";
    ret_val = 0;
  }

  return ret_val;

}// end masa_map_coeff

// user chooses steady or transient solution
int MASA::masa_map_temporal(string input_string, string return_string)
{
  multimap<string, string> names;
  string temp_string;
  int ret_val=0;  // exits with zero upon failure, 1 for success

  // convert to lower case
  temp_string=uptolow(input_string);
  temp_string=remove_line(temp_string);
  temp_string=remove_whitespace(temp_string);

  // map to steady state solutions
  names.insert(pair<string, string>("steady",                  "steady-state"));
  names.insert(pair<string, string>("steadystate",             "steady-state"));
  names.insert(pair<string, string>("stdystate",               "steady-state"));

  // map to time dependent function
  names.insert(pair<string, string>("transient",               "time-dependent"));
  names.insert(pair<string, string>("timedependent",           "time-dependent"));
  names.insert(pair<string, string>("temporal",                "time-dependent"));
  names.insert(pair<string, string>("timevarying",             "time-dependent"));

  multimap<string, string>::iterator p;

  p = names.find(temp_string);
  if(p != names.end()) { // found a name
    do {
      cout << temp_string << ", " << p->second;
      cout << endl;
      p++;
    } while (p != names.upper_bound(temp_string));
    ret_val=1;
  }
  else{
    cout << "Name not found.\n";
    ret_val = 0;
  }
  
  return ret_val;

}// end masa_map_temporal

// get string, parse it, and return the string name of the solution we want
int MASA::masa_map(string input_string, string return_string)
{
  int error; // simple error handler
  int n=0;
  int pos_old=0;
  string line;
  vector<string> str2;
  char* special_char="_";  
  
  // return array of strings
  // we break up the string between each underscore, and save that into a vector
  for(int pos=input_string.find(special_char);pos<string::npos;pos=input_string.find(special_char,pos+1))
    {
      line=input_string.substr(pos_old,pos-pos_old);
      str2.push_back(line);
      pos_old=pos+1;
      n++;
    }
  // we are still missing the last piece -- you break code by ending with an underscore?
  line=input_string.substr(pos_old);
  str2.push_back(line);
  n++;
  
  // loop over each element in the array, constructing pieces of the solution name
  for(int i=0;i<n;i++)
    {
      error=0;
      error+=masa_map_solution (str2[i],return_string);
      error+=masa_map_temporal (str2[i],return_string);
      error+=masa_map_dimension(str2[i],return_string);
      error+=masa_map_coeff    (str2[i],return_string);
      
      if(error>1) // error handler: this implies that redundant mappings exist
	{
	  cout << "\n MASA: Massive failure -- redundant mappings exist!\n";
	  cout << "\n Please email nick@ices.utexas.edu about this.\n";
	  exit(1);
	}
      
    }// end for loop
  
  return 0;
}
