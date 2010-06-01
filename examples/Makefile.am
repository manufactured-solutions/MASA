INCLUDES                = -I$(top_srcdir)/src -I$(top_builddir)/src
bin_PROGRAMS            = tester MASAshell euler_example

MASAshell_SOURCES       = masa_shell.cpp
MASAshell_LDADD         = -L$(top_srcdir)/src -lmasa

euler_example_SOURCES   = euler_example.cpp
euler_example_LDADD     = -L$(top_srcdir)/src -lmasa

tester_SOURCES       = tester.cpp
tester_LDADD         = -L$(top_srcdir)/src -lmasa

#MASAshell_CFLAGS = -I $(top_srcdir) -I $(top_builddir)/src
