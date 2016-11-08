CLEANFILES  = .license.stamp *.gcda *.gcno

h_sources  = masa.h masa_internal.h smasa.h ad_masa.h                        \
             nsctpl_fwd.hpp nsctpl.hpp

cc_sources = masa_core.cpp masa_class.cpp masa_map.cpp heat.cpp euler.cpp    \
             cns.cpp sod.cpp cmasa.cpp axi_euler.cpp axi_cns.cpp rans_sa.cpp \
             euler_chem.cpp euler_transient.cpp radiation.cpp fans_sa.cpp    \
             ablation.cpp cp_normal.cpp nsctpl.cpp laplace.cpp

cc_sources += burgers_equation.cpp
cc_sources += euler_transient_2d.cpp
cc_sources += euler_transient_3d.cpp
cc_sources += axi_euler_transient.cpp
cc_sources += axi_cns_transient.cpp
cc_sources += ad_cns_2d_crossterms.cpp
cc_sources += ad_cns_3d_crossterms.cpp
cc_sources += convdiff_steady_nosource_1d.cpp
cc_sources += navierstokes_3d_incompressible.cpp
cc_sources += navierstokes_3d_incompressible_homogeneous.cpp
cc_sources += navierstokes_3d_transient_sutherland.cpp
# do not edit this line! --l33t--

# License Information

BUILT_SOURCES = .license.stamp

lib_LTLIBRARIES         = libmasa.la
library_includedir      = $(includedir)
library_include_HEADERS = masa.h

#-----------------------
# MASA C/C++ library
#-----------------------

libmasa_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libmasa_la_CPPFLAGS     = $(masa_optional_INCLUDES)
libmasa_la_SOURCES      = $(cc_sources) $(h_sources)

#----------------------
# MASA Fortran library
# (Installs standalone)
#----------------------
if FORT_ENABLED
  CLEANFILES             += masa.mod 
  lib_LTLIBRARIES        += libfmasa.la
  libfmasa_la_LDFLAGS     = $(all_libraries) -release $(GENERIC_RELEASE)
  libfmasa_la_CPPFLAGS    = $(masa_optional_INCLUDES)
  libfmasa_la_SOURCES     = masa.f90 $(libmasa_la_SOURCES)
endif

#-----------------------
# MASA (SWIG) Python interface (optional)
#-----------------------

if SWIG_ENABLED
  BUILT_SOURCES        += masa_wrap.c masa.py
  SWIG_SRC              = $(srcdir)/masa.i

  pkgpython_PYTHON      = masa.py __init__.py
  pkgpyexec_LTLIBRARIES = _masa.la
  _masa_la_SOURCES      = masa_wrap.c
  _masa_la_CFLAGS       = $(AX_SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
  _masa_la_LDFLAGS      = -module
  _masa_la_LIBADD       = libmasa.la

  # Produce two built sources from a single masa.i source
  masa_wrap.c : $(SWIG_SRC)
	$(SWIG) $(AX_SWIG_PYTHON_OPT) $(srcdir)/masa.i
  masa.py     : masa_wrap.c
	@if test ! -f $@; then rm -f $+; $(MAKE) $(AM_MAKEFLAGS) $+; fi

  CLEANFILES  += masa_wrap.c masa.py masa.pyc
else
  # Stubs for "make distcheck" to be happy
  masa_wrap.c :
	@touch $@
  masa.py     :
	@touch $@

endif

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

if FORT_ENABLED

install-data-local: libfmasa.la
	mkdir -p $(DESTDIR)$(libdir)
	$(INSTALL_DATA) $(top_builddir)/src/masa.mod $(DESTDIR)$(libdir)/masa.mod

uninstall-local: libfmasa.la
	rm -rf $(DESTDIR)$(libdir)/masa.mod
endif

EXTRA_DIST   = lic_utils/update_license.pl .license.stamp masa.i
