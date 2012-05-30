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
	           vec                 \
	           pass_func           \
		   purge               \
		   laplace             \
	           heat_const_steady   \
		   heat_var_steady     \
	           heat_const_unsteady \
                   heat_var_unsteady   \
	           euler1d             \
	           euler2d             \
                   euler3d             \
                   euler_transient_1d  \
                   euler_chem_1d       \
	           cns2d               \
	           cns3d               \
	           cns3d_phys          \
	           cn2d3d              \
	           ablation            \
                   axi_euler           \
                   axi_cns             \
                   rans_sa             \
                   fans_sa_wall        \
                   radiation           \
                   sod                 \
                   cp_normal           \
	           c_misc              \
	           c_source            \
	           c_purge             \
	           c_array             \
	           c_laplace           \
	           c_heat1dsc          \
	           c_heat2dsc          \
	           c_heat3dsc          \
	           c_euler1d           \
	           c_euler2d           \
                   c_euler3d           \
                   c_euler_chem_1d     \
	           c_cns2d             \
	           c_cns3d             \
	           c_radiation         \
                   f_misc              \
                   f_array             \
                   f_purge             \
                   f_laplace           \
                   f_heat              \
                   f_euler1d           \
                   f_euler2d           \
                   f_euler3d           \
                   f_euler_chem_1d     \
                   f_cns2d             \
                   f_cns3d             \
                   f_radiation         \
		   ad_euler            \
		   ad_cns              \
		   ad_cns_crossterm

#----------------------
# C++ Test Sources
#----------------------

misc_SOURCES                = misc.cpp
poly_SOURCES                = poly.cpp
fail_cond_SOURCES           = fail_cond.cpp
catch_exception_SOURCES     = catch_exception.cpp
register_SOURCES            = register.cpp
uninit_SOURCES              = uninit.cpp
vec_SOURCES                 = vec.cpp
pass_func_SOURCES           = pass_func.cpp
purge_SOURCES               = purge.cpp
laplace_SOURCES             = laplace.cpp
heat_const_steady_SOURCES   = heatsc.cpp
heat_var_steady_SOURCES     = heatsv.cpp
heat_const_unsteady_SOURCES = heatuc.cpp
heat_var_unsteady_SOURCES   = heatuv.cpp
axi_euler_SOURCES           = axi_euler.cpp
ablation_SOURCES            = ablation.cpp
axi_cns_SOURCES             = axi_cns.cpp
euler1d_SOURCES             = euler1d.cpp
euler2d_SOURCES             = euler2d.cpp
euler3d_SOURCES             = euler3d.cpp
euler_transient_1d_SOURCES  = euler_transient_1d.cpp
euler_chem_1d_SOURCES       = euler_chem_1d.cpp
cns2d_SOURCES               = cns2d.cpp
cns3d_SOURCES               = cns3d.cpp
cns3d_phys_SOURCES          = cns3d_phys.cpp
cn2d3d_SOURCES              = cns-2d-3d.cpp
rans_sa_SOURCES             = rans_sa.cpp
fans_sa_wall_SOURCES        = fans_sa_wall.cpp
radiation_SOURCES           = radiation.cpp
sod_SOURCES                 = sod.cpp
cp_normal_SOURCES           = cp_normal.cpp

#----------------------
# C++ AD Test Sources
#----------------------
ad_euler_SOURCES            = ad_euler.cpp
ad_cns_SOURCES              = ad_cns.cpp
ad_cns_crossterm_SOURCES    = ad_cns_crossterm.cpp

#----------------------
# C Test Sources
#----------------------

c_misc_SOURCES              = c_misc.c
c_source_SOURCES            = c_source.c
c_purge_SOURCES             = c_purge.c
c_array_SOURCES             = c_array.c
c_laplace_SOURCES           = c_laplace.c
c_heat1dsc_SOURCES          = c_heat1dsc.c
c_heat2dsc_SOURCES          = c_heat2dsc.c
c_heat3dsc_SOURCES          = c_heat3dsc.c
c_euler1d_SOURCES           = c_euler1d.c
c_euler2d_SOURCES           = c_euler2d.c
c_euler3d_SOURCES           = c_euler3d.c
c_euler_chem_1d_SOURCES     = c_euler_chem_1d.c
c_cns2d_SOURCES             = c_cns2d.c
c_cns3d_SOURCES             = c_cns3d.c
c_radiation_SOURCES         = c_radiation.c

#----------------------
# Fortran Test Sources
#----------------------

f_misc_SOURCES              = f_misc.F90
f_misc_LDADD                = -lfmasa

f_array_SOURCES             = f_array.F90
f_array_LDADD               = -lfmasa

f_purge_SOURCES             = f_purge.F90
f_purge_LDADD               = -lfmasa

f_laplace_SOURCES           = f_laplace.F90
f_laplace_LDADD             = -lfmasa

f_heat_SOURCES              = f_heat.F90 heat_source.c 
f_heat_LDADD                = -lfmasa

f_euler1d_SOURCES           = f_euler_source.F90 f_euler1d.F90 euler_source.c
f_euler1d_LDADD             = -lfmasa

f_euler2d_SOURCES           = f_euler_source.F90 f_euler2d.F90  euler_source.c 
f_euler2d_LDADD             = -lfmasa

f_euler3d_SOURCES           = f_euler_source.F90 f_euler3d.F90 euler_source.c
f_euler3d_LDADD             = -lfmasa

f_euler_chem_1d_SOURCES     = f_euler_source.F90 f_euler_chem_1d.F90 euler_source.c
f_euler_chem_1d_LDADD       = -lfmasa

f_cns2d_SOURCES             = f_cns_source.F90 f_cns2d.F90 cns_source.c
f_cns2d_LDADD               = -lfmasa

f_cns3d_SOURCES             = f_cns_source.F90 f_cns3d.F90 cns_source.c
f_cns3d_LDADD               = -lfmasa

f_radiation_SOURCES         = f_radiation.F90
f_radiation_LDADD           = -lfmasa

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
			     vec                  \
			     pass_func            \
			     purge                \
			     laplace              \
			     heat_const_steady    \
			     heat_var_steady      \
	                     heat_const_unsteady  \
                             heat_var_unsteady    \
	                     euler1d              \
	                     euler2d              \
	                     euler3d              \
                             euler_transient_1d   \
                             euler_chem_1d        \
	                     cns2d                \
                             cns3d                \
                             cns3d_phys           \
	       	             cn2d3d               \
	       	             ablation             \
                             rans_sa              \
                             fans_sa_wall         \
                             radiation            \
                             sod                  \
                             cp_normal            \
			     c_init.sh            \
	                     c_misc               \
	                     c_array              \
	                     c_source             \
	                     c_purge              \
	                     c_laplace            \
	                     c_heat1dsc           \
	                     c_heat2dsc           \
	                     c_heat3dsc           \
	                     c_euler1d            \
	                     c_euler2d            \
	                     c_euler3d            \
                             c_euler_chem_1d      \
                             c_cns2d              \
                             c_cns3d              \
                             c_radiation          \
			     f_init.sh            \
		             f_misc               \
	                     f_array              \
		             f_purge              \
		             f_laplace            \
		             f_heat               \
		             f_euler1d            \
		             f_euler2d            \
		             f_euler3d            \
		             f_euler_chem_1d      \
		             f_cns2d              \
		             f_cns3d              \
		             f_radiation          \
			     ad_init.sh           \
			     ad_euler             \
			     ad_cns               \
			     ad_cns_crossterm     \
                             finalize.sh



EXTRA_DIST                 = init.sh              \
			     c_init.sh            \
	                     f_init.sh            \
	                     ad_init.sh           \
                             finalize.sh          \
			     tests.h              \
			     .license.stamp       

# License Information
.license.stamp: $(top_srcdir)/LICENSE
	for file in $(srcdir)/*.h $(srcdir)/*.c $(srcdir)/*.cpp ; do \
		$(top_srcdir)/src/lic_utils/update_license.pl \
			$(top_srcdir)/LICENSE $$file | tee -a $@ ;\
	done
	for file in $(srcdir)/*.F90 ; do \
		$(top_srcdir)/src/lic_utils/update_license.pl \
			--c2f_comment $(top_srcdir)/LICENSE $$file | tee -a $@ ;\
	done
	echo 'updated source license headers' >$@

CLEANFILES = .license.stamp *.mod *.gcda *.gcno
