AM_CPPFLAGS        = -I$(top_srcdir)/src -I$(top_builddir)/src $(masa_optional_INCLUDES)
AM_FCFLAGS         = $(FC_STRICT_TOL_CPP)
check_PROGRAMS     =
dist_check_SCRIPTS =
EXTRA_DIST         = tests.h init.sh ad_init.sh c_init.sh \
                     f_init.sh finalize.sh py_init.sh \
                     touch_swig_interface.sh test_swig.py

#--------------
# C++ Binaries
#--------------
TESTS_CXX       =
check_PROGRAMS += $(TESTS_CXX)

TESTS_CXX                   +=  misc
misc_SOURCES                 =  misc.cpp
misc_LDADD                   =  ../src/libmasa.la

TESTS_CXX                   +=  poly
poly_SOURCES                 =  poly.cpp
poly_LDADD                   =  ../src/libmasa.la

TESTS_CXX                   +=  fail_cond
fail_cond_SOURCES            =  fail_cond.cpp
fail_cond_LDADD              =  ../src/libmasa.la

TESTS_CXX                   +=  catch_exception
catch_exception_SOURCES      =  catch_exception.cpp
catch_exception_LDADD        =  ../src/libmasa.la

TESTS_CXX                   +=  register
register_SOURCES             =  register.cpp
register_LDADD               =  ../src/libmasa.la

TESTS_CXX                   +=  uninit
uninit_SOURCES               =  uninit.cpp
uninit_LDADD                 =  ../src/libmasa.la

TESTS_CXX                   +=  vec
vec_SOURCES                  =  vec.cpp
vec_LDADD                    =  ../src/libmasa.la

TESTS_CXX                   +=  pass_func
pass_func_SOURCES            =  pass_func.cpp
pass_func_LDADD              =  ../src/libmasa.la

TESTS_CXX                   +=  purge
purge_SOURCES                =  purge.cpp
purge_LDADD                  =  ../src/libmasa.la

TESTS_CXX                   +=  laplace
laplace_SOURCES              =  laplace.cpp
laplace_LDADD                =  ../src/libmasa.la

TESTS_CXX                   +=  heat_steady_const
heat_steady_const_SOURCES    =  heatsc.cpp
heat_steady_const_LDADD      =  ../src/libmasa.la

TESTS_CXX                   +=  heat_steady_var
heat_steady_var_SOURCES      =  heatsv.cpp
heat_steady_var_LDADD        =  ../src/libmasa.la

TESTS_CXX                   +=  heat_unsteady_const
heat_unsteady_const_SOURCES  =  heatuc.cpp
heat_unsteady_const_LDADD    =  ../src/libmasa.la

TESTS_CXX                   +=  heat_unsteady_var
heat_unsteady_var_SOURCES    =  heatuv.cpp
heat_unsteady_var_LDADD      =  ../src/libmasa.la

TESTS_CXX                   +=  axi_euler
axi_euler_SOURCES            =  axi_euler.cpp
axi_euler_LDADD              =  ../src/libmasa.la

TESTS_CXX                   +=  ablation
ablation_SOURCES             =  ablation.cpp
ablation_LDADD               =  ../src/libmasa.la

TESTS_CXX                   +=  axi_cns
axi_cns_SOURCES              =  axi_cns.cpp
axi_cns_LDADD                =  ../src/libmasa.la

TESTS_CXX                   +=  euler1d
euler1d_SOURCES              =  euler1d.cpp
euler1d_LDADD                =  ../src/libmasa.la

TESTS_CXX                   +=  euler2d
euler2d_SOURCES              =  euler2d.cpp
euler2d_LDADD                =  ../src/libmasa.la

TESTS_CXX                   +=  euler3d
euler3d_SOURCES              =  euler3d.cpp
euler3d_LDADD                =  ../src/libmasa.la

TESTS_CXX                   +=  euler_transient_1d
euler_transient_1d_SOURCES   =  euler_transient_1d.cpp
euler_transient_1d_LDADD     =  ../src/libmasa.la

TESTS_CXX                   +=  euler_chem_1d
euler_chem_1d_SOURCES        =  euler_chem_1d.cpp
euler_chem_1d_LDADD          =  ../src/libmasa.la

TESTS_CXX                   +=  cns2d
cns2d_SOURCES                =  cns2d.cpp
cns2d_LDADD                  =  ../src/libmasa.la

TESTS_CXX                   +=  cns3d
cns3d_SOURCES                =  cns3d.cpp
cns3d_LDADD                  =  ../src/libmasa.la

TESTS_CXX                   +=  cns3d_phys
cns3d_phys_SOURCES           =  cns3d_phys.cpp
cns3d_phys_LDADD             =  ../src/libmasa.la

TESTS_CXX                   +=  cns-2d-3d
cns_2d_3d_SOURCES            =  cns-2d-3d.cpp
cns_2d_3d_LDADD              =  ../src/libmasa.la

TESTS_CXX                   +=  rans_sa
rans_sa_SOURCES              =  rans_sa.cpp
rans_sa_LDADD                =  ../src/libmasa.la

TESTS_CXX                   +=  fans_sa_wall
fans_sa_wall_SOURCES         =  fans_sa_wall.cpp
fans_sa_wall_LDADD           =  ../src/libmasa.la

TESTS_CXX                   +=  radiation
radiation_SOURCES            =  radiation.cpp
radiation_LDADD              =  ../src/libmasa.la

TESTS_CXX                   +=  sod
sod_SOURCES                  =  sod.cpp
sod_LDADD                    =  ../src/libmasa.la

TESTS_CXX                   +=  cp_normal
cp_normal_SOURCES            =  cp_normal.cpp
cp_normal_LDADD              =  ../src/libmasa.la


#-----------------
# C++ AD Binaries
#-----------------
TESTS_CXX_AD    =
check_PROGRAMS += $(TESTS_CXX_AD)

TESTS_CXX_AD             +=  ad_euler
ad_euler_SOURCES          =  ad_euler.cpp
ad_euler_LDADD            =  ../src/libmasa.la

TESTS_CXX_AD             +=  ad_cns
ad_cns_SOURCES            =  ad_cns.cpp
ad_cns_LDADD              =  ../src/libmasa.la

TESTS_CXX_AD             +=  ad_ins
ad_ins_SOURCES            =  ad_ins.cpp
ad_ins_LDADD              =  ../src/libmasa.la

TESTS_CXX_AD             +=  ad_cns_crossterm
ad_cns_crossterm_SOURCES  =  ad_cns_crossterm.cpp
ad_cns_crossterm_LDADD    =  ../src/libmasa.la

TESTS_CXX_AD             +=  ad_divgradsimple
ad_divgradsimple_SOURCES  =  ad_divgradsimple.cpp
ad_divgradsimple_LDADD    =  ../src/libmasa.la

#------------
# C Binaries
#------------
TESTS_C         =
check_PROGRAMS += $(TESTS_C)

TESTS_C                 +=  c_misc
c_misc_SOURCES           =  c_misc.c
c_misc_LDADD             =  ../src/libmasa.la

TESTS_C                 +=  c_source
c_source_SOURCES         =  c_source.c
c_source_LDADD           =  ../src/libmasa.la

TESTS_C                 +=  c_purge
c_purge_SOURCES          =  c_purge.c
c_purge_LDADD            =  ../src/libmasa.la

TESTS_C                 +=  c_array
c_array_SOURCES          =  c_array.c
c_array_LDADD            =  ../src/libmasa.la

TESTS_C                 +=  c_laplace
c_laplace_SOURCES        =  c_laplace.c
c_laplace_LDADD          =  ../src/libmasa.la

TESTS_C                 +=  c_heat1dsc
c_heat1dsc_SOURCES       =  c_heat1dsc.c
c_heat1dsc_LDADD         =  ../src/libmasa.la

TESTS_C                 +=  c_heat2dsc
c_heat2dsc_SOURCES       =  c_heat2dsc.c
c_heat2dsc_LDADD         =  ../src/libmasa.la

TESTS_C                 +=  c_heat3dsc
c_heat3dsc_SOURCES       =  c_heat3dsc.c
c_heat3dsc_LDADD         =  ../src/libmasa.la

TESTS_C                 +=  c_euler1d
c_euler1d_SOURCES        =  c_euler1d.c
c_euler1d_LDADD          =  ../src/libmasa.la

TESTS_C                 +=  c_euler2d
c_euler2d_SOURCES        =  c_euler2d.c
c_euler2d_LDADD          =  ../src/libmasa.la

TESTS_C                 +=  c_euler3d
c_euler3d_SOURCES        =  c_euler3d.c
c_euler3d_LDADD          =  ../src/libmasa.la

TESTS_C                 +=  c_euler_chem_1d
c_euler_chem_1d_SOURCES  =  c_euler_chem_1d.c
c_euler_chem_1d_LDADD    =  ../src/libmasa.la

TESTS_C                 +=  c_cns2d
c_cns2d_SOURCES          =  c_cns2d.c
c_cns2d_LDADD            =  ../src/libmasa.la

TESTS_C                 +=  c_cns3d
c_cns3d_SOURCES          =  c_cns3d.c
c_cns3d_LDADD            =  ../src/libmasa.la

TESTS_C                 +=  c_radiation
c_radiation_SOURCES      =  c_radiation.c
c_radiation_LDADD        =  ../src/libmasa.la



#------------------
# Fortran Binaries
#------------------
if FORT_ENABLED
TESTS_F         =
check_PROGRAMS += $(TESTS_F)

TESTS_F                 +=  f_misc
f_misc_SOURCES           =  f_misc.F90
f_misc_LDADD             =  ../src/libfmasa.la

TESTS_F                 +=  f_array
f_array_SOURCES          =  f_array.F90
f_array_LDADD            =  ../src/libfmasa.la

TESTS_F                 +=  f_purge
f_purge_SOURCES          =  f_purge.F90
f_purge_LDADD            =  ../src/libfmasa.la

TESTS_F                 +=  f_laplace
f_laplace_SOURCES        =  f_laplace.F90
f_laplace_LDADD          =  ../src/libfmasa.la

TESTS_F                 +=  f_heat
f_heat_SOURCES           =  f_heat.F90 heat_source.c
f_heat_LDADD             =  ../src/libfmasa.la

TESTS_F                 +=  f_euler1d
f_euler1d_SOURCES        =  f_euler_source.F90 f_euler1d.F90 euler_source.c
f_euler1d_LDADD          =  ../src/libfmasa.la

TESTS_F                 +=  f_euler2d
f_euler2d_SOURCES        =  f_euler_source.F90 f_euler2d.F90 euler_source.c
f_euler2d_LDADD          =  ../src/libfmasa.la

TESTS_F                 +=  f_euler3d
f_euler3d_SOURCES        =  f_euler_source.F90 f_euler3d.F90 euler_source.c
f_euler3d_LDADD          =  ../src/libfmasa.la

TESTS_F                 += f_euler_chem_1d
f_euler_chem_1d_SOURCES  = f_euler_source.F90 f_euler_chem_1d.F90 euler_source.c
f_euler_chem_1d_LDADD    = ../src/libfmasa.la

TESTS_F                 +=  f_cns2d
f_cns2d_SOURCES          =  f_cns_source.F90 f_cns2d.F90 cns_source.c
f_cns2d_LDADD            =  ../src/libfmasa.la

TESTS_F                 +=  f_cns3d
f_cns3d_SOURCES          =  f_cns_source.F90 f_cns3d.F90 cns_source.c
f_cns3d_LDADD            =  ../src/libfmasa.la

TESTS_F                 +=  f_radiation
f_radiation_SOURCES      =  f_radiation.F90
f_radiation_LDADD        =  ../src/libfmasa.la

endif

#------------------
# Regression suite
#------------------

TESTS_ENVIRONMENT   =
TESTS               = init.sh
TESTS              += $(TESTS_CXX)

dist_check_SCRIPTS += ad_init.sh
TESTS              += ad_init.sh
TESTS              += $(TESTS_CXX_AD)

dist_check_SCRIPTS += c_init.sh
TESTS              += c_init.sh
TESTS              += $(TESTS_C)

# If optional fortran interfaces are enabled:
if FORT_ENABLED
dist_check_SCRIPTS += f_init.sh
TESTS              += f_init.sh
TESTS              += $(TESTS_F)
endif

dist_check_SCRIPTS += finalize.sh
TESTS              += finalize.sh

# If optional SWIG python interfaces are enabled:
if SWIG_ENABLED
TESTS              += py_init.sh touch_swig_interface.sh
dist_check_SCRIPTS += py_init.sh test_swig.py
endif


#---------------------
# License Information
#---------------------
BUILT_SOURCES  = .license.stamp
EXTRA_DIST    += .license.stamp
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
