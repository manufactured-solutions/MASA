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

GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_RELEASE)

#------------------------------
# Checks for required programs
#------------------------------

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_F77
AC_PROG_LIBTOOL
AM_SANITY_CHECK

# Generate Output Files

AC_OUTPUT( Makefile src/Makefile examples/Makefile tests/Makefile)