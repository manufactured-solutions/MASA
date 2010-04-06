## Process this file with automake to produce Makefile.in

AM_CFLAGS     = -I$(top_srcdir) -I$(top_builddir)
AM_CFLAGS    += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_FCFLAGS   += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_LDFLAGS         =
bin_PROGRAMS       = # Append below
check_PROGRAMS     = # Append below
lib_LTLIBRARIES    = # Append below
noinst_HEADERS     = # Append below
noinst_LTLIBRARIES = # Append below
pkginclude_HEADERS = # Append below
CLEANFILES         = # Append below

# Use MPI compilers
CC=@MPICC@
FC=@MPIFC@

# Headers installed with the library
pkginclude_HEADERS += masa.h

# Public MASA library containing C and Fortran interfaces
# Need distinct C and Fortran sublibraries for libtool --tag purposes
lib_LTLIBRARIES     += libmasa.la
libmasa_la_LDFLAGS   = -release $(PACKAGE_VERSION)
libmasa_la_SOURCES   = # Required even if empty
libmasa_la_LIBADD    = libmasa_c.la libmasa_f.la

# C sources to place into the public library
noinst_LTLIBRARIES   += libmasa_c.la
libmasa_c_la_SOURCES  = masa.c

# Fortran sources to place into the public library
noinst_LTLIBRARIES        += libmasa_f.la
libmasa_f_la_SOURCES       = masaf.f90
libmasa_f_la_LIBTOOLFLAGS  = --tag=FC

# Fortran module installation for masaf
moduledir   = @libdir@
module_DATA = masaf.mod
CLEANFILES += masaf.mod

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
