## Process this file with automake to produce Makefile.in

AM_CFLAGS       = -I$(top_srcdir) -I$(top_builddir)
AM_CFLAGS      += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS      = -I$(top_srcdir) -I$(top_builddir)
AM_FCFLAGS     += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS     += -I$(top_srcdir)/esio -I$(top_builddir)/esio
AM_LDFLAGS      = # Append below
bin_PROGRAMS    = # Append below
check_PROGRAMS  = # Append below
noinst_HEADERS  = # Append below
noinst_PROGRAMS = # Append below

# Use MPI compilers
CC=@MPICC@
FC=@MPIFC@

## Pick up HPCT library installation
AM_CFLAGS  += $(HPCT_CFLAGS)
AM_FCFLAGS += $(HPCT_CFLAGS)/../lib
AM_LDFLAGS += $(HPCT_LIBS)

## Sample test program
check_PROGRAMS += sample
sample_SOURCES  = sample.c
sample_LDADD    = ../esio/libesio.la

## Fortran test program to read/write a field
noinst_PROGRAMS        += simplef90
simplef90_SOURCES       = header.f90 simplef90.f90
simplef90_LDADD         = ../esio/libesio.la
simplef90_LIBTOOLFLAGS  = --tag=FC

TESTS_ENVIRONMENT =
TESTS = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
