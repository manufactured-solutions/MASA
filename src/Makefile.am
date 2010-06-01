h_sources               = masa.h  masa_internal.h
cc_sources              = masa_core.cpp masa_class.cpp masa_map.cpp masa_heat.cpp masa_euler.cpp masa_cns.cpp

INCLUDES                = 
AM_FCFLAGS              = $(INCLUDES)
AM_CFLAGS               = $(INCLUDES)
AM_CPPFLAGS             = $(INCLUDES) 

if USE_INTEL
  AM_CPPFLAGS          += -fp-model precise
endif

lib_LTLIBRARIES         = libmasa.la
library_includedir      = $(includedir)
library_include_HEADERS = masa.h
libmasa_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
libmasa_la_SOURCES      = $(cc_sources) $(h_sources)

#----------------------------
# Standalone Binaries
#----------------------------

