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

noinst_PROGRAMS  = verify MASAshell display_solutions euler_example navierstokes_example heat_example switch rans_sa_example euler_chem euler_transient
noinst_PROGRAMS += c_heat_example c_euler_example 
noinst_PROGRAMS += f_euler f_cns

# ----------------- 
#    CPP
# ----------------- 

verify_SOURCES               = verify.cpp cg.cpp
MASAshell_SOURCES            = masa_shell.cpp
display_solutions_SOURCES    = display_solutions.cpp
heat_example_SOURCES         = heat-eq.cpp
euler_example_SOURCES        = euler_example.cpp
euler_transient_SOURCES      = euler_transient.cpp
navierstokes_example_SOURCES = compressible_navier_stokes_example.cpp
rans_sa_example_SOURCES      = rans_sa.cpp
euler_chem_SOURCES           = euler_chem.cpp
switch_SOURCES               = switch.cpp

example_DATA                 = euler_example.cpp heat-eq.cpp compressible_navier_stokes_example.cpp

# ----------------- 
#    C
# ----------------- 

c_heat_example_SOURCES       = c_heat_example.c
c_euler_example_SOURCES      = c_euler_example.c

example_DATA                += c_heat_example.c c_euler_example.c


# ----------------- 
#    Fortran
# ----------------- 

f_euler_SOURCES              = f_euler.f90
f_euler_LDADD                = -lfmasa
f_euler_LIBTOOLFLAGS         = --tag=FC

f_cns_SOURCES                = f_cns.f90
f_cns_LDADD                  = -lfmasa
f_cns_LIBTOOLFLAGS           = --tag=FC

example_DATA                += f_euler.f90

# ----------------- 
# License
# ----------------- 
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/lic_utils/update_license.pl $(top_srcdir)/LICENSE *.c *.cpp
	echo 'updated source license headers' >$@

CLEANFILES = .license.stamp *.mod *.gcda *.gcno
