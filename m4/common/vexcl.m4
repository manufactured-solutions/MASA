# SYNOPSIS
#
#   Test for VEXCL
#
#   AX_PATH_VEXCL( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-vexcl=DIR option. Searches --with-vexcl,
#   $VEXCL_DIR, and the usual places for VEXCL headers and libraries.
#
#   On success, sets VEXCL_CPPFLAGS, VEXCL_LIBS, and #defines HAVE_VEXCL.
#   Also defines automake conditional VEXCL_ENABLED.  Assumes package
#   is optional unless overridden with $2=yes.
#
# LAST MODIFICATION
#
#   $Id: vexcl.m4 -1   $
#
# COPYLEFT
#
#   Copyright (c) 2013 Roy H. Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2012 Paul T. Bauman <pbauman@ices.utexas.edu>
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

AC_DEFUN([AX_PATH_VEXCL],
[

AC_ARG_VAR(VEXCL_DIR,[root directory of VEXCL installation])

AC_ARG_WITH(vexcl,
  [AS_HELP_STRING([--with-vexcl[=DIR]],[root directory of VEXCL installation (default = VEXCL_DIR)])],
  [with_vexcl=$withval
if test "${with_vexcl}" != yes; then
    VEXCL_PREFIX=$withval
elif test "x${VEXCL_DIR}" != "x"; then
   VEXCL_PREFIX=${VEXCL_DIR}
else
    VEXCL_PREFIX=/usr
fi
],[
with_vexcl=$withval
if test "x${VEXCL_DIR}" != "x"; then
   VEXCL_PREFIX=${VEXCL_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_VEXCL=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_vexcl}" != no ; then

    # If we can see the vexcl headers, then we know where to get them
    # and we'll need C++11 to compile them
    if test -e "${VEXCL_PREFIX}/vexcl/vexcl.hpp" ; then
       VEXCL_CPPFLAGS="-I${VEXCL_PREFIX}"
       AX_CXX_COMPILE_STDCXX_11(noext)
    fi

    ac_VEXCL_save_CPPFLAGS="$CPPFLAGS"

    CPPFLAGS="${VEXCL_CPPFLAGS} ${CPPFLAGS}"

    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([vexcl/vexcl.hpp],[found_header=yes],[found_header=no])
    AC_LANG_POP([C++])

    # Make sure we have OpenCL support
    AX_CHECK_CL([C++])

    # Newer VexCL requires new Boost linkage.
    # There doesn't seem to be any way to simply *test* for Boost with
    # boost.m4, so at this point we'll assume that if we've seen VexCL
    # then you really wanted to use VexCL.
    if test x$found_header = xyes; then
      BOOST_REQUIRE([1.47]) # Chrono introduced in 1.47
      BOOST_CHRONO
      BOOST_DATE_TIME
      BOOST_FILESYSTEM
      BOOST_SYSTEM
      BOOST_THREADS

      VEXCL_CPPFLAGS="$VEXCL_CPPFLAGS $BOOST_CPPFLAGS"
      VEXCL_LDFLAGS="$VEXCL_LDFLAGS $BOOST_CHRONO_LDFLAGS $BOOST_DATE_TIME_LDFLAGS $BOOST_FILESYSTEM_LDFLAGS $BOOST_SYSTEM_LDFLAGS"
      VEXCL_LIBS="$VEXCL_LIBS $BOOST_CHRONO_LIBS $BOOST_DATE_TIME_LIBS $BOOST_FILESYSTEM_LIBS $BOOST_SYSTEM_LIBS"
    fi

    #-----------------------
    # Minimum version check skipped - there's no versioning
    # information in vexcl headers as of 0.7.0
    #----------------------

    CPPFLAGS="$ac_VEXCL_save_CPPFLAGS"

    succeeded=yes
    if test "$found_header" != yes; then
       succeeded=no
    fi
    if test "$no_cl" = yes; then
       succeeded=no
    else
       VEXCL_CPPFLAGS="$VEXCL_CPPFLAGS $CL_CFLAGS"
       VEXCL_LIBS="$VEXCL_LIBS $CL_LIBS"
    fi

    if test "$succeeded" = no; then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([VEXCL not found.  Try either --with-vexcl or setting VEXCL_DIR.])
       else
          AC_MSG_NOTICE([optional VEXCL library not found])
          VEXCL_CPPFLAGS="" # VEXCL_CPPFLAGS empty on failure
          VEXCL_LDFLAGS="" # VEXCL_LDFLAGS empty on failure
          VEXCL_LIBS="" # VEXCL_LIBS empty on failure
       fi
    else
        HAVE_VEXCL=1
        AC_DEFINE(HAVE_VEXCL,1,[Define if VEXCL is available])
        AC_SUBST(VEXCL_CPPFLAGS)
        AC_SUBST(VEXCL_LDFLAGS)
        AC_SUBST(VEXCL_LIBS)
        AC_SUBST(VEXCL_PREFIX)
    fi

    AC_SUBST(HAVE_VEXCL)

# fi

AM_CONDITIONAL(VEXCL_ENABLED,test x$HAVE_VEXCL = x1)

])
