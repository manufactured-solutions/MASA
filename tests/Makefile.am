INCLUDES                  = -I$(top_srcdir)/src -I$(top_builddir)/src

bin_PROGRAMS              = heat_const_steady \
			    heat_var_steady \
	                    heat_const_unsteady \
                            heat_var_unsteady \
	                    euler2d \
                            euler3d \
	                    navierstokes3d \
	                    navierstokes2d

# Test Programs:
heat_const_steady_SOURCES   = heatsc.cpp
heat_const_steady_LDADD     = -L$(top_srcdir)/src -lmasa

heat_var_steady_SOURCES     = heatsv.cpp
heat_var_steady_LDADD       = -L$(top_srcdir)/src -lmasa

heat_const_unsteady_SOURCES = heatuc.cpp
heat_const_unsteady_LDADD   = -L$(top_srcdir)/src -lmasa

heat_var_unsteady_SOURCES   = heatuv.cpp
heat_var_unsteady_LDADD     = -L$(top_srcdir)/src -lmasa

euler2d_SOURCES             = euler2d.cpp Euler_2d_source.C
euler2d_LDADD               = -L$(top_srcdir)/src -lmasa

euler3d_SOURCES             = euler3d.cpp Euler_3d_source.C
euler3d_LDADD               = -L$(top_srcdir)/src -lmasa

navierstokes2d_SOURCES      = ns2d.cpp Navier-Stokes_2d_compressible_source.C
navierstokes2d_LDADD        = -L$(top_srcdir)/src -lmasa

navierstokes3d_SOURCES      = ns3d.cpp Navier-Stokes_3d_compressible_source.C
navierstokes3d_LDADD        = -L$(top_srcdir)/src -lmasa

TESTS_ENVIRONMENT         = 
TESTS                     = heat_const_steady \
			    heat_var_steady \
	                    heat_const_unsteady \
                            heat_var_unsteady \
	                    euler2d \
	                    euler3d \
                            navierstokes3d \
                            navierstokes2d