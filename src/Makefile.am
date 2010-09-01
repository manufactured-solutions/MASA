h_sources               = masa.h cmasa.h masa_internal.h
cc_sources              = masa_core.cpp masa_class.cpp masa_map.cpp heat.cpp euler.cpp \
			  cns.cpp sod.cpp cmasa.cpp axi_euler.cpp axi_cns.cpp

# License Information
BUILT_SOURCES = .license.stamp

INCLUDES                = 
AM_FCFLAGS              = $(INCLUDES)
AM_CFLAGS               = $(INCLUDES)
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = libmasa.la
library_includedir      = $(includedir)
library_include_HEADERS = masa.h cmasa.h
libmasa_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libmasa_la_SOURCES      = $(cc_sources) $(h_sources)

bin_PROGRAMS            = version
version_SOURCES         = version.cpp
version_LDADD           = -L$(top_builddir)/src -lmasa
#version_DEPENDENCIES    = libmasa.la

EXTRA_DIST              = lic_utils/update_license.pl

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(cc_sources) $(h_sources)
	echo 'updated source license headers' >$@

CLEANFILES = .license.stamp
