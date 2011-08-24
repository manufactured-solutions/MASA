# SYNOPSIS
#
#   Builds python interfaces
#   AX_SWIG()
#
# DESCRIPTION
#
#   Provides --enable-swig for python interface generation
#
# LAST MODIFICATION
#
#   $Id: 
#
# COPYLEFT
#
#   Copyright (c) 2011 Nicholas Malaya <nick@ices.utexas.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.
#

AC_DEFUN([AX_SWIG],
[

AC_ARG_ENABLE([swig], AC_HELP_STRING([--enable-swig],[enable SWIG python interface generation]),
 	       MASA_SWIG=1 
	       AC_DEFINE(MASA_SWIG,1,[Define if SWIG python inteface generation enabled]),[])

if test "$MASA_SWIG" = "1"; then

   ax_pkg_swig(1,[echo 'found it'],[echo 'didnt find it'])

fi


])