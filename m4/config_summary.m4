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
#   2009-07-16
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version............... : $PACKAGE-$VERSION
echo
echo Fortran compiler.............. : $FC
echo Debug mode.................... : $enable_debug
echo Fortran compiler flags........ : $FCFLAGS
#echo libGRVY DIR................... : $GRVY_PREFIX
echo Install dir................... : $prefix 
echo Build user.................... : $USER
echo Build host.................... : $BUILD_HOST
echo Configure date................ : $BUILD_DATE
echo Build architecture............ : $BUILD_ARCH
echo SVN revision number........... : $BUILD_VERSION

#echo
#echo Optional Features:
# if test $LINK_PETSC -eq 0; then
#   echo '   'Link with PETSc............ : no
# else
#   echo '   'Link with PETSc............ : yes
#   echo '   'PETSC DIR.................. : $PETSC_DIR
# fi

echo
echo '-------------------------------------------------------------------------------'

echo
echo Configure complete, now type \'make\' and then \'make install\'.
echo

])
