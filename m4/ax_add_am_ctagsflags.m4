# SYNOPSIS
#
#   AX_ADD_AM_CTAGSFLAGS
#
# DESCRIPTION
#
#  Using the AX_AM_MACROS framework, specify AM_CTAGSFLAGS with nice
#  options for every Automake file containing @INC_AMINCLUDE@.
#
# COPYLEFT
#
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_ADD_AM_CTAGSFLAGS],[
AC_REQUIRE([AX_AM_MACROS])
AX_ADD_AM_MACRO([
# Added by AX_ADD_AM_CTAGSFLAGS macro
AM_CTAGSFLAGS=--c++-kinds=+p --fields=+iaS --extra=+qf
])
])
