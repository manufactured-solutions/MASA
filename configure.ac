# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)

#init
AC_INIT([masa],[0.50.0], [masa-dev@ices.utexas.edu])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

# snarf and provide versioning numbers
AX_SPLIT_VERSION
GENERIC_MAJOR_VERSION=$AX_MAJOR_VERSION
GENERIC_MINOR_VERSION=$AX_MINOR_VERSION
GENERIC_MICRO_VERSION=$AX_POINT_VERSION
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION 
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)
AC_SUBST(GENERIC_RELEASE)

# define
PACKAGE_DESCRIPTION="Manufactured Analytical Solutions Abstraction Library"
AC_SUBST([PACKAGE_DESCRIPTION])
PACKAGE_URL="https://red.ices.utexas.edu/projects/software/"
AC_SUBST([PACKAGE_URL])

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
AX_CXX_MINOPT	# request minimum optimization level for CXX

# ---------------------------------------------
# Enable/Disable -Wall (Warn all)
# ---------------------------------------------

AX_MASA_WARN_ALL

# ------------------------------
# Checks for required programs
# ------------------------------

AC_DISABLE_STATIC

AC_PROG_CC
AC_PROG_CXX

# ----------------------------
# We may need the math library
# ----------------------------
AC_LANG_PUSH([C])
AC_SEARCH_LIBS([cos], [m], [], [
    AC_MSG_ERROR([unable to find the cos() function])
])
AC_LANG_POP([C])

# ---------------------------------------------
# enable fortran interfaces
# ---------------------------------------------
FORTRAN_INTERFACES=0
AC_ARG_ENABLE([fortran-interfaces], AC_HELP_STRING([--enable-fortran-interfaces],[enable fortran support]),
[
dnl set suffix to .F90 to get CPP support for testing vendor Fortran compiler
AC_FC_SRCEXT([F90]) 
dnl request minimum optimization level for FC
AX_FC_MINOPT
AC_PROG_F77()
AC_PROG_FC()
FORTRAN_INTERFACES=1
AC_DEFINE(FORTRAN_INTERFACES,1,[Define if Fortran interfaces enabled])
],[])
AM_CONDITIONAL(FORT_ENABLED,[test "x$FORTRAN_INTERFACES" = x1]) 

dnl--------------------------
dnl Checks for third-party libraries
dnl--------------------------
masa_optional_INCLUDES=""

dnl Enable MetaPhysicL if a new enough version is found
AX_PATH_METAPHYSICL(0.0.0, no)
if (test x$HAVE_METAPHYSICL = x1); then
  masa_optional_INCLUDES="$METAPHYSICL_CPPFLAGS $masa_optional_INCLUDES"
fi

#
# now call libtool
#
AC_PROG_LIBTOOL
AM_SANITY_CHECK


# ---------------------------------------------
# enable python interfaces with SWIG
# ---------------------------------------------
SWIG_INTERFACES=0
AC_ARG_ENABLE([python-interfaces], AC_HELP_STRING([--enable-python-interfaces],[enable python support]),
[
AM_PATH_PYTHON(2.3)
dnl AC_PROG_SWIG(1.3.21)
AC_PROG_SWIG()
dnl SWIG_ENABLE_C
SWIG_PYTHON
SWIG_INTERFACES=1
AC_DEFINE(SWIG_INTERFACES,1,[Define if python interfaces enabled])
],[])

AM_CONDITIONAL(SWIG_ENABLED,[test "x$SWIG_INTERFACES" = x1]) 

# ---------------------------------------------
# Check for extra-strict regression tolerance
# ---------------------------------------------

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

# --------------------------
# Checks for code coverage
# -------------------------

AX_CODE_COVERAGE

# Query configuration environment

AX_SUMMARIZE_ENV

# ----------------
# Doxygen support 
# ---------------- 

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)

DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN(MASA, doxygen/masa.dox, docs) 

# Run scripts for regression tests
AC_CONFIG_FILES(tests/py_init.sh, [chmod +x tests/py_init.sh])
AC_CONFIG_FILES(tests/touch_swig_interface.sh, [chmod +x tests/touch_swig_interface.sh])

# Generate Output Files
AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile
                 tests/Makefile tutorial/Makefile
                 doxygen/Makefile doxygen/txt_common/about_vpath.page
                 src/masa.h masa.spec masa.pc])
AC_OUTPUT()

# FINAL SUMMARY

AX_SUMMARIZE_CONFIG 
