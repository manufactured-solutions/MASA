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
// $Author$
// $Id$
//
// masa_map.cpp: various helper functions and input parsing
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h>
#include <cstring>
#include <locale>

using namespace MASA;

// -------------------------------------------------------
// these are just a few utility functions
// -------------------------------------------------------

// all this does is convert a string to entirely lower case
// should simplify our user stuff
void uptolow(std::string& str) 
{
  for (int i=0;i<str.length();i++) 
    str[i] = std::tolower(str[i]);
}

// remove dashes from strings
void remove_line(std::string& str)
{
  int position = str.find( "-" ); // find first space
  while ( position != std::string::npos )
    {
      str.replace( position, 1, "" );
      position = str.find( "-", position + 1 );
    }
}

void remove_whitespace(std::string& str)
{
  int position = str.find( " " ); // find first space
  while ( position != std::string::npos )
    {
      str.replace( position, 1, "" );
      position = str.find( " ", position + 1 );
    }
}

int MASA::masa_map(std::string* input_string)
{
  // fix up the input given
  std::string temp;
  temp = *input_string;
  uptolow(temp);
  remove_line(temp);
  remove_whitespace(temp);
  *input_string=temp;
  return 0;

}//done with masa_map
