INCLUDES               = -I$(top_srcdir)/src -I$(top_builddir)/src

bin_PROGRAMS           = euler2d euler3d navierstokes3d

euler2d_SOURCES        = euler2d.cpp Euler_2d_source.C
euler2d_LDADD          = -L$(top_srcdir)/src -lmasa

euler3d_SOURCES        = euler3d.cpp Euler_3d_source.C
euler3d_LDADD          = -L$(top_srcdir)/src -lmasa

navierstokes3d_SOURCES = ns3d.cpp Navier-Stokes_3d_compressible_source.C
navierstokes3d_LDADD   = -L$(top_srcdir)/src -lmasa

TESTS_ENVIRONMENT      = 
TESTS                  = euler2d euler3d navierstokes3d
