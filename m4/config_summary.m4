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
#   $Id: config_summary.m4 40370 2013-07-11 13:16:42Z nick $
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version.................   : $PACKAGE-$VERSION
echo				      
echo C++ compiler....................   : $CXX
echo C++ compiler flags..............   : $CXXFLAGS
echo C compiler......................   : $CC
echo C compiler flags................   : $CFLAGS
echo Install dir.....................   : $prefix 
echo Build user......................   : $USER
echo Build host......................   : $BUILD_HOST
echo Configure date..................   : $BUILD_DATE
echo Build architecture..............   : $BUILD_ARCH
echo Source control revision.........   : $BUILD_VERSION
echo
echo Optional Features:
   if test "$HAVE_GCOV_TOOLS" = "0"; then
     echo '   'Enable gcov code coverage.... : no
   else
     echo '   'Enable gcov code coverage.... : yes
   fi
   if test "x$HAVE_METAPHYSICL" = "x1"; then
     echo '   'MetaPhysicL.................. : yes
     echo '      'METAPHYSICL_CPPFLAGS...... : $METAPHYSICL_CPPFLAGS
     echo '      'METAPHYSICL_LDFLAGS....... : $METAPHYSICL_LDFLAGS
     echo '      'METAPHYSICL_LIBS.......... : $METAPHYSICL_LIBS
   else
     echo '   'MetaPhysicL.................. : no
   fi
   if test "$MASA_STRICT_REGRESSION" = "1"; then
     echo '   'Enable absolute error tests.. : yes
   else
     echo '   'Enable absolute error tests.. : no
   fi

   if test "$MASA_WARN_ALL" = "1"; then
     echo '   'Enable all warnings.......... : yes
   else
     echo '   'Enable all warnings.......... : no
   fi

   if test "$SWIG_INTERFACES" = "1"; then
     echo '   'Enable python interfaces..... : yes
   else
     echo '   'Enable python interfaces..... : no
   fi

   if test "$FORTRAN_INTERFACES" = "1"; then
     echo '   'Enable fortran interfaces.... : yes
     echo Fortran compiler................   : $FC
     echo Fortran compiler flags..........   : $FCFLAGS
   else
     echo '   'Enable fortran interfaces.... : no
   fi

echo
echo '-------------------------------------------------------------------------------'

echo
echo Configure complete, now type \'make\' and then \'make install\'.
echo
echo To verify your verification library, type \'make check\' 
echo to run a suite of regression tests.

])
