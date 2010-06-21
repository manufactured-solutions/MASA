INCLUDES      = -I$(top_srcdir)/src -I$(top_builddir)/src
AM_FCFLAGS    = $(INCLUDES)
AM_CFLAGS     = $(INCLUDES)
AM_CPPFLAGS   = $(INCLUDES)
LIBS          = -L$(top_builddir)/src -lmasa

bin_PROGRAMS  = tester MASAshell euler_example navierstokes_example heat_example switch
#bin_PROGRAMS += C_euler_example

# ----------------- 
#    CPP
# ----------------- 

MASAshell_SOURCES            = masa_shell.cpp
heat_example_SOURCES         = heat-eq.cpp
euler_example_SOURCES        = euler_example.cpp
navierstokes_example_SOURCES = compressible_navier_stokes_example.cpp
switch_SOURCES               = switch.cpp
tester_SOURCES               = tester.cpp

# ----------------- 
#    C
# ----------------- 

#C_euler_example_SOURCES       = ceuler_example.c