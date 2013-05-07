h_sources  = masa.h masa_internal.h smasa.h ad_masa.h                        \
             nsctpl_fwd.hpp nsctpl.hpp                                       \
	     dualnumber.h numberarray.h dualnumberarray.h compare_types.h    \
	     raw_type.h shadownumber.h dualshadowarray.h dualshadow.h        \
	     testable.h

cc_sources = masa_core.cpp masa_class.cpp masa_map.cpp heat.cpp euler.cpp    \
             cns.cpp sod.cpp cmasa.cpp axi_euler.cpp axi_cns.cpp rans_sa.cpp \
             euler_chem.cpp euler_transient.cpp radiation.cpp fans_sa.cpp    \
             ablation.cpp cp_normal.cpp nsctpl.cpp laplace.cpp

cc_sources += burgers_equation.cpp
cc_sources += euler_transient_2d.cpp
cc_sources += euler_transient_3d.cpp
cc_sources += axi_euler_transient.cpp
cc_sources += axi_cns_transient.cpp
cc_sources += navierstokes_3d_transient_sutherland.cpp
cc_sources += ad_cns_2d_crossterms.cpp
cc_sources += ad_cns_3d_crossterms.cpp
cc_sources += convdiff_steady_nosource_1d.cpp
cc_sources += navierstokes_3d_incompressible.cpp
# do not edit this line! --l33t--

# License Information

BUILT_SOURCES = .license.stamp

lib_LTLIBRARIES         = libmasa.la libfmasa.la
library_includedir      = $(includedir)
library_include_HEADERS = masa.h

#-----------------------
# MASA C/C++ library
#-----------------------

libmasa_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libmasa_la_SOURCES      = $(cc_sources) $(h_sources)

#-----------------------
# MASA Fortran interface
#-----------------------

libfmasa_la_LDFLAGS     = $(all_libraries) -release $(GENERIC_RELEASE)
libfmasa_la_SOURCES     = masa.f90

#-----------------------
# Standalone binaries
#-----------------------

bin_PROGRAMS                           = @PACKAGE_NAME@_version
@PACKAGE_NAME@_version_SOURCES         = version.cpp masa.h
@PACKAGE_NAME@_version_LDADD           = -L$(top_builddir)/src -lmasa
@PACKAGE_NAME@_version_DEPENDENCIES    = libmasa.la

#-----------------------
# Licensing updates
#-----------------------

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/lic_utils/update_license.pl -S=$(top_srcdir)/src $(top_srcdir)/LICENSE $(cc_sources) $(h_sources)
	$(top_srcdir)/src/lic_utils/update_license.pl -S=$(top_srcdir)/src --c2f_comment $(top_srcdir)/LICENSE masa.f90
	echo 'updated source license headers' >$@

#------------------------
# install fo'tran modules
#------------------------

install-data-local: libfmasa.la
	mkdir -p $(DESTDIR)$(libdir)
	$(INSTALL_DATA) $(top_builddir)/src/masa.mod $(DESTDIR)$(libdir)/masa.mod

uninstall-local: libfmasa.la
	rm -rf $(DESTDIR)$(libdir)/masa.mod

EXTRA_DIST  = lic_utils/update_license.pl .license.stamp
CLEANFILES  = .license.stamp masa.mod *.gcda *.gcno
