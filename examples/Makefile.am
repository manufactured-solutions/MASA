INCLUDES      = -I$(top_srcdir)/src -I$(top_builddir)/src
AM_FCFLAGS    = $(INCLUDES)
AM_CFLAGS     = $(INCLUDES)
AM_CPPFLAGS   = $(INCLUDES)
LIBS          = -L$(top_builddir)/src -lmasa

# License Information
BUILT_SOURCES = .license.stamp

# cpp 
bin_PROGRAMS  = verify MASAshell euler_example navierstokes_example heat_example switch rans_sa_example #tester 

# c-code 
bin_PROGRAMS += c_heat_example c_euler_example 

#fortran 
bin_PROGRAMS += f_euler f_cns

# ----------------- 
#    CPP
# ----------------- 

verify_SOURCES               = verify.cpp cg.cpp
MASAshell_SOURCES            = masa_shell.cpp

heat_example_SOURCES         = heat-eq.cpp
euler_example_SOURCES        = euler_example.cpp
navierstokes_example_SOURCES = compressible_navier_stokes_example.cpp
rans_sa_example_SOURCES      = rans_sa.cpp
switch_SOURCES               = switch.cpp
#tester_SOURCES               = tester.cpp

# ----------------- 
#    C
# ----------------- 

c_heat_example_SOURCES       = c_heat_example.c
c_euler_example_SOURCES      = c_euler_example.c


# ----------------- 
#    Fortran
# ----------------- 

f_euler_SOURCES         = f_euler.f90 masa.f90
test_LIBTOOLFLAGS       = --tag=FC

f_cns_SOURCES         = f_cns.f90 masa.f90
test_LIBTOOLFLAGS       = --tag=FC

# ----------------- 
# License
# ----------------- 
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/lic_utils/update_license.pl $(top_srcdir)/LICENSE *.c *.cpp
	echo 'updated source license headers' >$@
