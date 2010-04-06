# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(MASA, 0.10.0, nick@ices.utexas.edu)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE(masa,0.10.0)

# Release versioning

GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=27
GENERIC_MICRO_VERSION=0

#GENERIC_API_VERSION=1.0
#AC_SUBST(GENERIC_API_VERSION)

#shared library versioning
#AC_SUBST(GENERIC_LIBRARY_VERSION)

#PACKAGE=$GENERIC_LIBRARY_NAME
#AC_SUBST(GENERIC_LIBRARY_NAME)
#GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION

GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_RELEASE)

# Required for some macros

#AX_AM_MACROS

#------------------------------
# Checks for required programs
#------------------------------

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_F77
AC_PROG_LIBTOOL
AM_SANITY_CHECK

#----------------------
# Checks for libraries.
#----------------------

#---------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
#---------------------------------------------------------------

AC_C_CONST
AC_HEADER_STDBOOL
AX_TLS

#----------------
# Check for boost
#----------------

AC_LANG([C++])
BOOST_REQUIRE([1.36])
BOOST_MATH
#BOOST_PROGRAM_OPTIONS
AC_CACHE_SAVE

#------------------------------
# Checks for library functions.
#------------------------------

AC_FUNC_MALLOC

# Query configuration environment

AX_SUMMARIZE_ENV

#----------------
# Doxygen support
#----------------

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN(MASA, doxygen/masa.dox, docs)

AC_OUTPUT( Makefile lib/Makefile masa/Makefile masa.pc masa.spec doc/Makefile)

# Final summary

AX_SUMMARIZE_CONFIG