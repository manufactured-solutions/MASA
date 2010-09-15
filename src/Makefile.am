h_sources               = masa.h masa_internal.h
cc_sources              = masa_core.cpp masa_class.cpp masa_map.cpp heat.cpp euler.cpp \
			  cns.cpp sod.cpp cmasa.cpp axi_euler.cpp axi_cns.cpp rans_sa.cpp

# License Information

BUILT_SOURCES = .license.stamp

INCLUDES                = 
AM_FCFLAGS              = $(INCLUDES)
AM_CFLAGS               = $(INCLUDES)
AM_CXXFLAGS             = $(INCLUDES)

lib_LTLIBRARIES         = libmasa.la libfmasa.la
library_includedir      = $(includedir)
library_include_HEADERS = masa.h

# Main C/C++ library

libmasa_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libmasa_la_SOURCES      = $(cc_sources) $(h_sources)

# Fortran Interface library

libfmasa_la_LDFLAGS     = $(all_libraries) -release $(GENERIC_RELEASE)
libfmasa_la_SOURCES     = masa.f90

# Standalone binaries

bin_PROGRAMS            = version
version_SOURCES         = version.cpp masa.h libmasa.la
version_LDADD           = -L$(top_builddir)/src -lmasa
version_DEPENDENCIES    = libmasa.la

#--------------------------                                                                                                    
# Fortran Module Interface                                                                                                     
#--------------------------                                                                                                    

noinst_PROGRAMS         = test
standalonedir           = @prefix@/bin

test_SOURCES            = main.f90
test_LDADD              = -L$(top_builddir)/src -lfmasa -lmasa
test_DEPENDENCIES       = libmasa.la

main.f90:

EXTRA_DIST              = lic_utils/update_license.pl

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(cc_sources) $(h_sources)
	echo 'updated source license headers' >$@

moduledir   = @libdir@
module_DATA = masa.mod
CLEANFILES  = .license.stamp masa.mod