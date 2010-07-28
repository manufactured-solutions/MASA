INCLUDES                = -I$(top_srcdir)/src -I$(top_builddir)/src
AM_CXXFLAGS             = $(INCLUDES) -O0

if USE_INTEL
  AM_CXXFLAGS          += -fp-model precise
endif

if USE_GNU
  AM_CXXFLAGS          += -fno-unsafe-math-optimizations -Wunsafe-loop-optimizations
endif

bin_PROGRAMS              = heat_const_steady \
			    heat_var_steady \
	                    heat_const_unsteady \
                            heat_var_unsteady \
	                    euler1d \
	                    euler2d \
                            euler3d \
	                    navierstokes2d \
	                    navierstokes3d \
	                    c_euler1d \
	                    c_euler2d \
                            c_euler3d \
	                    c_navierstokes2d \
	                    c_navierstokes3d 

#----------------------
# Test Programs --- CPP: 
#----------------------
heat_const_steady_SOURCES   = heatsc.cpp
heat_const_steady_LDADD     = -L$(top_srcdir)/src -lmasa

heat_var_steady_SOURCES     = heatsv.cpp
heat_var_steady_LDADD       = -L$(top_srcdir)/src -lmasa

heat_const_unsteady_SOURCES = heatuc.cpp
heat_const_unsteady_LDADD   = -L$(top_srcdir)/src -lmasa

heat_var_unsteady_SOURCES   = heatuv.cpp
heat_var_unsteady_LDADD     = -L$(top_srcdir)/src -lmasa

euler1d_SOURCES             = euler1d.cpp Euler_1d_source.C
euler1d_LDADD               = -L$(top_srcdir)/src -lmasa

euler2d_SOURCES             = euler2d.cpp Euler_2d_source.C
euler2d_LDADD               = -L$(top_srcdir)/src -lmasa

euler3d_SOURCES             = euler3d.cpp Euler_3d_source.C
euler3d_LDADD               = -L$(top_srcdir)/src -lmasa

navierstokes2d_SOURCES      = ns2d.cpp Navier-Stokes_2d_compressible_source.C
navierstokes2d_LDADD        = -L$(top_srcdir)/src -lmasa

navierstokes3d_SOURCES      = ns3d.cpp Navier-Stokes_3d_compressible_source.C
navierstokes3d_LDADD        = -L$(top_srcdir)/src -lmasa

#----------------------
# Test Programs --- C: 
#----------------------
c_euler1d_SOURCES             = c_euler1d.c
c_euler1d_LDADD               = -L$(top_srcdir)/src -lmasa

c_euler2d_SOURCES             = c_euler2d.c
c_euler2d_LDADD               = -L$(top_srcdir)/src -lmasa

c_euler3d_SOURCES             = c_euler3d.c
c_euler3d_LDADD               = -L$(top_srcdir)/src -lmasa

c_navierstokes2d_SOURCES      = c_ns2d.c
c_navierstokes2d_LDADD        = -L$(top_srcdir)/src -lmasa

c_navierstokes3d_SOURCES      = c_ns3d.c
c_navierstokes3d_LDADD        = -L$(top_srcdir)/src -lmasa


#----------------------
# Test Queue 
#----------------------

TESTS_ENVIRONMENT         = 
TESTS                     = heat_const_steady \
			    heat_var_steady \
	                    heat_const_unsteady \
                            heat_var_unsteady \
	                    euler1d \
	                    euler2d \
	                    euler3d \
                            navierstokes2d\
                            navierstokes3d\
	                    c_euler1d \
	                    c_euler2d \
	                    c_euler3d \
                            c_navierstokes2d\
                            c_navierstokes3d 
