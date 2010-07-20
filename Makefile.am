include $(top_srcdir)/doxygen/aminclude.am

# Build in these directories:
AUTOMAKE_OPTIONS = foreign
SUBDIRS          = src examples tests doxygen

BUILT_SOURCES = dist_version

# Distribute these directories:
DIST_SUBDIRS = src examples

# generated documentation will be included in distributed archive
EXTRA_DIST= bootstrap $(DX_CONFIG)