include $(top_srcdir)/doxygen/aminclude.am

# Build in these directories:
AUTOMAKE_OPTIONS = foreign
SUBDIRS          = src examples tests doxygen

# Distribute these directories:
DISTDIRS = src examples
EXTRA_DIST= AUTHORS README LICENSE CHANGES dist_version

# Revision control support for external distribution

if SVN_CHECKOUT
dist_version: FORCE
	@SVN_REVISION@ > $(top_srcdir)/dist_version

FORCE:
endif