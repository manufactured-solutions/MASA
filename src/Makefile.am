h_sources  = masa.h masa_internal.h smasa.h                                  \
             nsctpl_fwd.hpp nsctpl.hpp
cc_sources = masa_core.cpp masa_class.cpp masa_map.cpp heat.cpp euler.cpp    \
             cns.cpp sod.cpp cmasa.cpp axi_euler.cpp axi_cns.cpp rans_sa.cpp \
             euler_chem.cpp euler_transient.cpp radiation.cpp fans_sa.cpp    \
             ablation.cpp cp_normal.cpp nsctpl.cpp

# License Information

BUILT_SOURCES = .license.stamp

INCLUDES                =
AM_FCFLAGS              = $(INCLUDES)
AM_CFLAGS               = $(INCLUDES)
AM_CXXFLAGS             = $(INCLUDES)

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
	echo 'updated source license headers' >$@

#------------------------
# install fo'tran modules
#------------------------

install-data-local: libfmasa.la
	mkdir -p $(DESTDIR)$(libdir)
	$(INSTALL_DATA) $(top_builddir)/src/masa.mod $(DESTDIR)$(libdir)/masa.mod

uninstall-local: libfmasa.la
	rm -rf $(INSTALL_DATA) $(top_builddir)/src/masa.mod $(DESTDIR)$(libdir)/masa.mod

EXTRA_DIST  = lic_utils/update_license.pl
CLEANFILES  = .license.stamp masa.mod *.gcda *.gcno
