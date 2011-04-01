# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)

# Release versioning
AC_INIT([masa], [0.30], [nick@ices.utexas.edu])
PACKAGE_DESCRIPTION="Manufactured Analytical Solutions Abstraction Library"
AC_SUBST([PACKAGE_DESCRIPTION])
PACKAGE_URL="https://red.ices.utexas.edu/projects/software/"
AC_SUBST([PACKAGE_URL])

GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=30
GENERIC_MICRO_VERSION=0

GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION 
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE

AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)
AC_SUBST(GENERIC_RELEASE)

# ------------------------------------------------------------------
# Optimization Settings:
# 
# As MASA is targeted to support verification efforts, we desire to
# maintain as much floating-point accuracy as possible (and are
# willing to sacrifice speed). We thus request minimal compiler
# optimization unless the user dictates otherwise via specific 
# CFLAGS/CXXFLAGS settings.
# ------------------------------------------------------------------

AX_CC_MINOPT	# request minimum optimization level for C
AX_CFLAGS_WARN_ALL
if test "$ax_cv_c_compiler_vendor" == "intel"; then
  CFLAGS="$CFLAGS -wd981 -wd383 -wd2259 -wd869 -wd1418"
fi

AX_CXX_MINOPT	# request minimum optimization level for CXX
AX_CXXFLAGS_WARN_ALL
if test "$ax_cv_cxx_compiler_vendor" == "intel"; then
  CXXFLAGS="$CXXFLAGS -wd981 -wd383 -wd2259 -wd869 -wd1418"
fi

# set suffix to .F90 to get CPP support for testing vendor Fortran compiler

AC_FC_SRCEXT([F90]) 

AX_FC_MINOPT	# request minimum optimization level for FC

# ------------------------------
# Checks for required programs
# ------------------------------

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_F77
AC_PROG_LIBTOOL
AX_WARNINGS_SANITIZE
AM_SANITY_CHECK

#---------------------------------------------
# Check for extra-strict regression tolerance
#---------------------------------------------

FC_STRICT_TOL_CPP=""

AC_ARG_ENABLE([strict-tests], AC_HELP_STRING([--enable-strict-tests],[enable strict absolute tolerance in regression tests]),
 	       MASA_STRICT_REGRESSION=1 
	       FC_STRICT_TOL_CPP="-DMASA_STRICT_REGRESSION"
	       AC_DEFINE(MASA_STRICT_REGRESSION,1,[Define if strict regression enabled]),[])

AC_SUBST(FC_STRICT_TOL_CPP)

#---------------------------------------------
# Check for exception handling
#---------------------------------------------

AC_ARG_ENABLE([masa-exceptions], AC_HELP_STRING([--enable-masa-exceptions],[enable exception handling in masa]),
 	       MASA_EXCEPTIONS=1 
	       AC_DEFINE(MASA_EXCEPTIONS,1,[Define if exception handling enabled]),[])

#--------------------------
# Checks for code coverage
#-------------------------

AX_CODE_COVERAGE

# Query configuration environment

AX_SUMMARIZE_ENV

# ----------------
# Doxygen support 
# ---------------- 

# DX_DOXYGEN_FEATURE(ON)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)

DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN(MASA, doxygen/masa.dox, docs) 

# Generate Output Files
AC_OUTPUT( Makefile src/Makefile examples/Makefile tests/Makefile doxygen/Makefile src/masa.h masa.spec masa.pc )

# FINAL SUMMARY

AX_SUMMARIZE_CONFIG 
