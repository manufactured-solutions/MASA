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
#   2010-04-06
#  By Nick, who does not necessarily know what he is doing here.

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version............... : $PACKAGE-$VERSION
echo
echo C++ compiler.................. : $CXX
echo Fortran compiler support...... : $FC
echo Install dir................... : $prefix 
echo Boost dir..................... : $BOOST_ROOT
echo Build user.................... : $USER
echo Build host.................... : $BUILD_HOST
echo Configure date................ : $BUILD_DATE
echo Build architecture............ : $BUILD_ARCH
echo SVN revision number........... : $BUILD_REV
echo SVN status.................... : $BUILD_STATUS
echo
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
