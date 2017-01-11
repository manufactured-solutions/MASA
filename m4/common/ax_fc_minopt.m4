#
# SYNOPSIS
#
#   AX_FC_MINOPT
#
# DESCRIPTION
#
#   Try to disable F90 compiler optimization(s). Macro based in part on
#   ax_cc_maxopt.m4. In order to allow for the user to override these
#   options and to avoid defaults imposed via AC_PROG_CC, it is
#   important to include this macro *before* AC_PROC_CC.
#
#   Requires macro(s): AX_COMPILER_VENDOR
#
# LAST MODIFICATION
#
#   $Id: ax_fc_minopt.m4 28022 2012-02-16 20:23:18Z karl $
#
# LICENSE
#
#   Copyright (c) 2012 Karl W. Schulz    <karl@ices.utexas.edu>
#   Copyright (c) 2011 Nicholas Malaya   <nick@ices.utexas.edu>
#   Copyright (c) 2010 Karl W. Schulz    <karl@ices.utexas.edu>
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_FC_MINOPT],
[

# Disable compiler optimizaitons if none specivied via CFLAGS

if test "$ac_test_FCFLAGS" != "set"; then

  AX_COMPILER_VENDOR

  # Baseline disable

  FCFLAGS="-O0"

  # Vendor specific settings for reduced optimization or
  # floating-point accuracy.

  case $ax_cv_fc_compiler_vendor in
    gnu) 
    	 FCFLAGS="$FCFLAGS -fno-unsafe-math-optimizations -Dgnu_compiler"
         ;;

    intel) 
    	 FCFLAGS="$FCFLAGS -fp-model precise -Dintel_compiler"
	 ;;

    portland)
 	 FCFLAGS="$FCFLAGS -Kieee -Mnofpapprox -Dportland_compiler"
	 ;;	

  esac	 
  
  AC_MSG_NOTICE([disabling Fortran compiler optimizations ($FCFLAGS)])

else 

    delim="[Info]:"

    echo "------------------------------------------------------------"
    echo "$delim Overriding minimum optimization request with"
    echo "$delim user-provided FCFLAGS ($FCFLAGS)."    
    echo "$delim Floating-point accuracy may be affected."
    echo "------------------------------------------------------------"

fi
])
