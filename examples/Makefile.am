INCLUDES                = -I$(top_srcdir)/src -I$(top_builddir)/src
bin_PROGRAMS            = tester MASAshell euler_example navierstokes_example heat_example

MASAshell_SOURCES       = masa_shell.cpp
MASAshell_LDADD         = -L$(top_srcdir)/src -lmasa

heat_example_SOURCES   = heat-eq.cpp
heat_example_LDADD     = -L$(top_srcdir)/src -lmasa

euler_example_SOURCES   = euler_example.cpp
euler_example_LDADD     = -L$(top_srcdir)/src -lmasa

navierstokes_example_SOURCES   = compressible_navier_stokes_example.cpp
navierstokes_example_LDADD     = -L$(top_srcdir)/src -lmasa

tester_SOURCES       = tester.cpp
tester_LDADD         = -L$(top_srcdir)/src -lmasa

#MASAshell_CFLAGS = -I $(top_srcdir) -I $(top_builddir)/src
