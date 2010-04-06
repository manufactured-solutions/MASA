# SYNOPSIS
#
#   Test for HDF5
#
#   AM_PATH_HDF5([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Provides a --with-hdf5=DIR option. Searches --with-hdf5,
#   $HDF5_DIR, and the usual places for HDF5 headers and libraries.
#
#   On success, sets HDF5_CFLAGS, HDF5_LIBS, and
#   #defines HAVE_HDF5.  When ACTION-IF-NOT-FOUND is not specified,
#   the default behavior is for configure to fail.
#
# LAST MODIFICATION
#
#   2010-03-28
#
# COPYLEFT
#
#   Copyright (c) 2010 Karl W. Schulz <karl@ices.utexas.edu>
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_HDF5],
[

AC_ARG_VAR(HDF5_HOME,[root directory of HDF5 installation])

AC_ARG_WITH(hdf5, 
  [AS_HELP_STRING([--with-hdf5[=DIR]],[root directory of HDF5 installation (default = HDF5_DIR)])],
  [with_hdf5=$withval
if test "${with_hdf5}" != yes; then
    HDF5_PREFIX=$withval
fi
],[
with_hdf5=$withval
if test "x${HDF5_DIR}" != "x"; then
   HDF5_PREFIX=${HDF5_DIR}
fi
])

if test "${with_hdf5}" != no ; then

    if test -d "${HDF5_PREFIX}/lib" ; then
       HDF5_LIBS="-L${HDF5_PREFIX}/lib -lhdf5 -Wl,-rpath,${HDF5_PREFIX}/lib"
    fi

    if test -d "${HDF5_PREFIX}/include" ; then
        HDF5_CFLAGS="-I${HDF5_PREFIX}/include"
    fi

    ac_HDF5_save_CFLAGS="$CFLAGS"
    ac_HDF5_save_CPPFLAGS="$CPPFLAGS"
    ac_HDF5_save_LDFLAGS="$LDFLAGS"
    ac_HDF5_save_LIBS="$LIBS"

    CFLAGS="${HDF5_CFLAGS} ${CFLAGS}"
    CPPFLAGS="${HDF5_CFLAGS} ${CPPFLAGS}"
    LDFLAGS="${HDF5_LIBS} ${LDFLAGS}"
    AC_LANG_PUSH([C])
    AC_CHECK_HEADER([hdf5.h],[found_header=yes],[found_header=no])
    AC_CHECK_LIB([hdf5],[H5Fopen],[found_library=yes],[found_library=no])
    AC_LANG_POP([C])

    CFLAGS="$ac_HDF5_save_CFLAGS"
    CPPFLAGS="$ac_HDF5_save_CPPFLAGS"
    LDFLAGS="$ac_HDF5_save_LDFLAGS"
    LIBS="$ac_HDF5_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$found_library" = yes; then
            succeeded=yes
        fi
    fi

    if test "$succeeded" = no; then
        ifelse([$2],,AC_MSG_ERROR([HDF5 not found.  Try either --with-hdf5 or setting HDF5_DIR.]),
            [$2])
    else
        AC_DEFINE(HAVE_HDF5,1,[Define if HDF5 is available])
        AC_SUBST(HDF5_CFLAGS)
        AC_SUBST(HDF5_LIBS)
	AC_SUBST(HDF5_PREFIX)
        ifelse([$1],,,[$1])
    fi

fi

])
