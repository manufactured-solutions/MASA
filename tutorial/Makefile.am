AM_CPPFLAGS   = -I$(top_srcdir)/src -I$(top_builddir)/src
LIBS          = -L$(top_builddir)/src -lmasa -lm

# License Information
BUILT_SOURCES = .license.stamp

noinst_PROGRAMS = laplacian

# ----------------- 
#    CPP
# ----------------- 

# ----------------- 
#    C
# ----------------- 

laplacian_SOURCES         = laplacian.c laplacian_utils.c cg.c laplacian.h

# ----------------- 
#    Pf'ortran
# ----------------- 

#---------------------------------
# Embedded license header support
#---------------------------------

.license.stamp: $(top_srcdir)/LICENSE
	for file in $(srcdir)/*.c $(srcdir)/*.cpp $(srcdir)/*.h; do \
		$(top_srcdir)/src/lic_utils/update_license.pl \
			$(top_srcdir)/LICENSE $$file | tee -a $@ ;\
	done
	for file in $(srcdir)/*.f90 $(srcdir)/*.F90 ; do \
		$(top_srcdir)/src/lic_utils/update_license.pl \
			--c2f_comment $(top_srcdir)/LICENSE $$file | tee -a $@ ;\
	done
	touch .license.stamp

CLEANFILES  = .license.stamp *.mod *.gcda *.gcno
