INCLUDES      = -I$(top_srcdir)/src -I$(top_builddir)/src
AM_FCFLAGS    = $(INCLUDES)
AM_CFLAGS     = $(INCLUDES)
AM_CPPFLAGS   = $(INCLUDES)
LIBS          = -L$(top_builddir)/src -lmasa

exampledir    = @prefix@/examples

# License Information
BUILT_SOURCES = .license.stamp

# binaries to insall

bin_PROGRAMS  = MASAshell

noinst_PROGRAMS  = MASAshell display_solutions euler_example navierstokes 
noinst_PROGRAMS += rans_sa fans_sa fans_sa_finite_d radiation euler_chem 
noinst_PROGRAMS += euler_transient heat_example switch ablation_example
noinst_PROGRAMS += c_heat_example c_euler_example 
noinst_PROGRAMS += f_cns f_euler 
noinst_PROGRAMS += cp_gaussian

# ----------------- 
#    CPP
# ----------------- 

MASAshell_SOURCES         = masa_shell.cpp
display_solutions_SOURCES = display_solutions.cpp
heat_example_SOURCES      = heat-eq.cpp
euler_example_SOURCES     = euler_example.cpp
ablation_example_SOURCES  = ablation_example.cpp
euler_transient_SOURCES   = euler_transient.cpp
navierstokes_SOURCES      = compressible_navier_stokes.cpp
rans_sa_SOURCES           = rans_sa.cpp
fans_sa_SOURCES           = fans_sa.cpp
fans_sa_finite_d_SOURCES  = fans_sa_finite_d.cpp
radiation_SOURCES         = radiation.cpp
euler_chem_SOURCES        = euler_chem.cpp
switch_SOURCES            = switch.cpp
example_DATA              = euler_example.cpp heat-eq.cpp compressible_navier_stokes.cpp


# ----------------- 
#    smasa
# ----------------- 

cp_gaussian_SOURCES       = cp_gaussian.cpp


# ----------------- 
#    C
# ----------------- 

c_heat_example_SOURCES    = c_heat_example.c
c_euler_example_SOURCES   = c_euler_example.c

example_DATA             += c_heat_example.c c_euler_example.c

# ----------------- 
#    Fortran
# ----------------- 

f_euler_SOURCES           = f_euler.f90
f_euler_LDADD             = -lfmasa
f_euler_LIBTOOLFLAGS      = --tag=FC

f_cns_SOURCES             = f_cns.f90
f_cns_LDADD               = -lfmasa
f_cns_LIBTOOLFLAGS        = --tag=FC

example_DATA             += f_euler.f90

# ----------------- 
# License
# ----------------- 
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/lic_utils/update_license.pl $(top_srcdir)/LICENSE *.c *.cpp
	echo 'updated source license headers' >$@

EXTRA_DIST  = examples.h
CLEANFILES = .license.stamp *.mod *.gcda *.gcno
