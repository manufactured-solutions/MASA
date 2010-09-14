# SYNOPSIS
#
#   Summarizes configuration settings.
#
#   AX_SUMMARIZE_CONFIG([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Outputs a summary of relevant configuration settings.
#
# LAST MODIFICATION
#
#   2010-7-20
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version............... : $PACKAGE-$VERSION
echo
echo C++ compiler.................. : $CXX
echo C++ compiler flags............ : $CXXFLAGS
echo C compiler.................... : $CC
echo C compiler flags.............. : $CFLAGS
echo Fortran compiler.............. : $FC
echo Fortran compiler flags........ : $FCFLAGS
# echo Debug mode.................... : $enable_debug
echo Install dir................... : $prefix 
echo Build user.................... : $USER
echo Build host.................... : $BUILD_HOST
echo Configure date................ : $BUILD_DATE
echo Build architecture............ : $BUILD_ARCH
echo SVN revision number........... : $BUILD_VERSION
echo
echo Optional Features:
   if test "$HAVE_GCOV_TOOLS" = "0"; then
     echo '   'Enable gcov code coverage.... : no
   else
     echo '   'Enable gcov code coverage.... : yes
   fi

   if test "$MASA_STRICT_REGRESSION" = "1"; then
     echo '   'Enable absolute error tests.. : yes
   else
     echo '   'Enable absolute error tests.. : no
   fi

echo
echo '-------------------------------------------------------------------------------'

echo
echo Configure complete, now type \'make\' and then \'make install\'.
echo
echo To verify your verification library, type \'make check\' 
echo to run a suite of regression tests.

])
