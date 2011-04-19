# SYNOPSIS
#
#   Enables -Wall (Warn all)
#   AX_MASA_WARN_ALL()
#
# DESCRIPTION
#
#   Provides a MASA --enable-warn-all 
#
# LAST MODIFICATION
#
#   $Id$: 
#
# COPYLEFT
#
#   Copyright (c) 2011 Nicholas Malaya <nick@ices.utexas.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.
#

AC_DEFUN([AX_MASA_WARN_ALL],
[

AC_ARG_ENABLE([warn-all], AC_HELP_STRING([--enable-warn-all],[enable all warnings]),
 	       MASA_WARN_ALL=1 
	       AC_DEFINE(MASA_WARN_ALL,1,[Define if strict warnings enabled]),[])

if test "$MASA_WARN_ALL" = "1"; then

   AX_WARNINGS_SANITIZE
   AX_CFLAGS_WARN_ALL

   # vendor specific warning settings

   if test "$ax_cv_c_compiler_vendor" == "intel"; then
     CFLAGS="$CFLAGS -wd981 -wd383 -wd2259 -wd869 -wd1418"
   fi

   if test "$ax_cv_c_compiler_vendor" == "portland"; then
     CFLAGS="$CFLAGS -Minform=inform"
   fi

   AX_CXXFLAGS_WARN_ALL

   if test "$ax_cv_cxx_compiler_vendor" == "gnu"; then
     CXXFLAGS="$CXXFLAGS -D_GLIBCXX_DEBUG"
   fi

   if test "$ax_cv_cxx_compiler_vendor" == "intel"; then
     CXXFLAGS="$CXXFLAGS -wd981 -wd383 -wd2259 -wd869 -wd1418"
   fi

   if test "$ax_cv_cxx_compiler_vendor" == "portland"; then
     CXXFLAGS="$CXXFLAGS -Minform=inform"
   fi

fi

])
