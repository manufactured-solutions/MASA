## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = foreign

# Build in these directories:
SUBDIRS = lib masa tests doc

# Distribute these items:
DIST_SUBDIRS = lib masa tests doc
EXTRA_DIST= bootstrap $(DX_CONFIG)

# Eliminate .svn directories in dist tarball
dist-hook:
	find $(distdir)/ -name '.svn' -type d -exec rm {} \;

# Tools in the auxillary directory
AUX_DIST = build-aux/install-sh    \
           build-aux/mkinstalldirs \
           build-aux/missing

# Support for pkgconfig
pkgconfigdir   = $(libdir)/pkgconfig
pkgconfig_DATA = masa.pc

# Files to be deleted by make clean variants
MOSTLYCLEANFILES = $(DX_CLEANFILES)
DISTCLEANFILES = aminclude.am
MAINTAINERCLEANFILES =     \
	aclocal.m4             \
	autom4te.cache/*       \
	$(AUX_DIST)            \
	config.log             \
	config.status          \
	config.sub             \
	configure              \
	Makefile.in            \
	src/config.h.in        \
	src/Makefile.in        \
	src/stamp.h.in

# Additional tasks run with the 'all' target
all-local:  libtool

# Documentation hooks in with the 'html' target
html-local: doxygen-run

# Make sure that libtool is up to date
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# Required for AX_PROG_DOXYGEN
include doxygen.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
