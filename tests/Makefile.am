INCLUDES         = -I$(top_srcdir)/src -I$(top_builddir)/src
LIBS             = -L$(top_builddir)/src -lmasa
AM_CFLAGS        = $(INCLUDES)
AM_CXXFLAGS      = $(INCLUDES)
AM_F77FLAGS      = $(INCLUDES)
AM_FCFLAGS       = $(INCLUDES) $(FC_STRICT_TOL_CPP)
AM_LDFLAGS       = $(LIBS)

# License Information
BUILT_SOURCES = .license.stamp

check_PROGRAMS   = misc                \
	           poly                \
	           fail_cond           \
	           catch_exception     \
	           register            \
	           uninit              \
	           heat_const_steady   \
		   heat_var_steady     \
	           heat_const_unsteady \
                   heat_var_unsteady   \
	           euler1d             \
	           euler2d             \
                   euler3d             \
	           ns2d                \
	           ns3d                \
	           ns3d_phys           \
	           n2d3d               \
                   axi_euler           \
                   axi_cns             \
                   rans_sa             \
                   sod                 \
	           c_misc              \
	           c_heat1dsc          \
	           c_heat2dsc          \
	           c_heat3dsc          \
	           c_euler1d           \
	           c_euler2d           \
                   c_euler3d           \
	           c_navierstokes2d    \
	           c_navierstokes3d    \
                   f_misc              \
                   f_heat              \
                   f_euler1d           \
                   f_euler2d           \
                   f_euler3d           \
                   f_cns2d             \
                   f_cns3d

#----------------------
# C++ Test Sources
#----------------------

misc_SOURCES                = misc.cpp
poly_SOURCES                = poly.cpp
fail_cond_SOURCES           = fail_cond.cpp
catch_exception_SOURCES     = catch_exception.cpp
register_SOURCES            = register.cpp
uninit_SOURCES              = uninit.cpp
heat_const_steady_SOURCES   = heatsc.cpp
heat_var_steady_SOURCES     = heatsv.cpp
heat_const_unsteady_SOURCES = heatuc.cpp
heat_var_unsteady_SOURCES   = heatuv.cpp
axi_euler_SOURCES           = axi_euler.cpp
axi_cns_SOURCES             = axi_cns.cpp
euler1d_SOURCES             = euler1d.cpp
euler2d_SOURCES             = euler2d.cpp
euler3d_SOURCES             = euler3d.cpp
ns2d_SOURCES                = ns2d.cpp
ns3d_SOURCES                = ns3d.cpp
ns3d_phys_SOURCES           = ns3d_phys.cpp
n2d3d_SOURCES               = ns-2d-3d.cpp
rans_sa_SOURCES             = rans_sa.cpp
sod_SOURCES                 = sod.cpp

#----------------------
# C Test Sources
#----------------------

c_misc_SOURCES              = c_misc.c
c_heat1dsc_SOURCES          = c_heat1dsc.c
c_heat2dsc_SOURCES          = c_heat2dsc.c
c_heat3dsc_SOURCES          = c_heat3dsc.c
c_euler1d_SOURCES           = c_euler1d.c
c_euler2d_SOURCES           = c_euler2d.c
c_euler3d_SOURCES           = c_euler3d.c
c_navierstokes2d_SOURCES    = c_ns2d.c
c_navierstokes3d_SOURCES    = c_ns3d.c

#----------------------
# Fortran Test Sources
#----------------------

f_misc_SOURCES              = f_misc.F90
f_misc_LDADD                = -lfmasa

f_heat_SOURCES              = f_heat.F90 heat_source.c 
f_heat_LDADD                = -lfmasa

f_euler1d_SOURCES           = f_euler1d.F90 euler_source.c
f_euler1d_LDADD             = -lfmasa

f_euler2d_SOURCES           = f_euler_source.F90 f_euler2d.F90  euler_source.c 
f_euler2d_LDADD             = -lfmasa

f_euler3d_SOURCES           = f_euler3d.F90 euler_source.c
f_euler3d_LDADD             = -lfmasa

f_cns2d_SOURCES             = f_cns2d.F90 cns_source.c
f_cns2d_LDADD               = -lfmasa

f_cns3d_SOURCES             = f_cns3d.F90 cns_source.c
f_cns3d_LDADD               = -lfmasa

#----------------------
# Regression Tests
#----------------------

TESTS_ENVIRONMENT          = 
TESTS                      = init.sh              \
			     misc                 \
			     fail_cond            \
                             catch_exception      \
			     register             \
			     poly                 \
			     uninit               \
			     heat_const_steady    \
			     heat_var_steady      \
	                     heat_const_unsteady  \
                             heat_var_unsteady    \
	                     euler1d              \
	                     euler2d              \
	                     euler3d              \
                             ns2d                 \
                             ns3d                 \
                             ns3d_phys            \
	       	             n2d3d                \
                             axi_euler            \
                             axi_cns              \
                             rans_sa              \
                             sod                  \
			     c_init.sh            \
	                     c_misc               \
	                     c_heat1dsc           \
	                     c_heat2dsc           \
	                     c_heat3dsc           \
	                     c_euler1d            \
	                     c_euler2d            \
	                     c_euler3d            \
                             c_navierstokes2d     \
                             c_navierstokes3d     \
			     f_init.sh            \
		             f_misc               \
		             f_heat               \
		             f_euler1d            \
		             f_euler2d            \
		             f_euler3d            \
		             f_cns2d              \
		             f_cns3d              \
                             finalize.sh

EXTRA_DIST                 = init.sh \
			     c_init.sh \
	                     f_init.sh \
                             finalize.sh

# License Information
.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/lic_utils/update_license.pl $(top_srcdir)/LICENSE *.c *.cpp
	$(top_srcdir)/src/lic_utils/update_license.pl --c2f_comment $(top_srcdir)/LICENSE *.F90
	echo 'updated source license headers' >$@

CLEANFILES = .license.stamp *.mod *.gcda *.gcno
