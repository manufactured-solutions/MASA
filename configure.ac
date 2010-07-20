# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([MASA], [0.10], [nick@ices.utexas.edu])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE

# Release versioning
GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=10
GENERIC_MICRO_VERSION=0
GENERIC_LIBRARY_NAME=MASA

PACKAGE=$GENERIC_LIBRARY_NAME
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION 
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

AC_SUBST(GENERIC_LIBRARY_NAME)   
AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_RELEASE)

# ------------------------------
# Checks for required programs
# ------------------------------

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_F77
AC_PROG_LIBTOOL

AM_CONDITIONAL([USE_INTEL],[test "$CXX" = "icpc"])
AM_CONDITIONAL([USE_GNU],[test "$CXX" = "g++"])
AM_SANITY_CHECK

# ------------------------------
# Flags
# ------------------------------

AM_CPPFLAGS="-O3 -g"
AM_FCFLAGS="$FCFLAGS -g -DDEBUG"

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
AC_OUTPUT( Makefile src/Makefile examples/Makefile tests/Makefile doxygen/Makefile)

# FINAL SUMMARY
AX_SUMMARIZE_ENV
AX_SUMMARIZE_CONFIG 
