# ===========================================================================
#         http://autoconf-archive.cryp.to/ax_check_enable_debug.html
# ===========================================================================
#
# SYNOPSIS
#
#   Check for the presence of an --enable-debug option to configure and
#   allow/avoid compiled debugging flags appropriately.
#
#   AX_CHECK_ENABLE_DEBUG([enable by default=yes/info/no])
#
# DESCRIPTION
#
#   Check for the presence of an --enable-debug option to configure, with the
#   specified default value used when the option is not present.  Return the
#   value in the variable $ax_enable_debug.
#
#   Specifying 'yes' adds '-g -O0' to the compilation flags for all languages.
#   Specifying 'info' adds '-g' to the compilation flags.  Otherwise, nothing
#   is added.  Define NDEBUG if debug disabled.  If debug not enabled, ensure
#   AC_PROG_* will not add debugging flags.  Should be invoked prior to any
#   AC_PROG_* compiler checks.
#
# LAST MODIFICATION
#
#   2009-11-02
#
# COPYLEFT
#
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_CHECK_ENABLE_DEBUG],[
    AC_BEFORE([$0],[AC_PROG_CC])dnl
    AC_BEFORE([$0],[AC_PROG_CXX])dnl
    AC_BEFORE([$0],[AC_PROG_F77])dnl
    AC_BEFORE([$0],[AC_PROG_FC])dnl

    AC_MSG_CHECKING(whether to enable debugging)
    ax_enable_debug_default="m4_tolower(m4_normalize(ifelse([$1],,[no],[$1])))"
    AC_ARG_ENABLE(debug,
        [AS_HELP_STRING([--enable-debug],[enable compiler debug flags])],
        [],enable_debug=$ax_enable_debug_default)
    if test "x$enable_debug" = "xyes"; then
        AC_MSG_RESULT(yes)
        CFLAGS="${CFLAGS} -g -O0"
        CXXFLAGS="${CXXFLAGS} -g -O0"
        FFLAGS="${FFLAGS} -g -O0"
        FCFLAGS="${FCFLAGS} -g -O0"
        OBJCFLAGS="${OBJCFLAGS} -g -O0"
    else
        if test "x$enable_debug" = "xinfo"; then
            AC_MSG_RESULT(info)
            CFLAGS="${CFLAGS} -g"
            CXXFLAGS="${CXXFLAGS} -g"
            FFLAGS="${FFLAGS} -g"
            FCFLAGS="${FCFLAGS} -g"
            OBJCFLAGS="${OBJCFLAGS} -g"
        else
            AC_MSG_RESULT(no)
            dnl Ensure AC_PROG_CC/CXX/F77/FC/OBJC will not enable debug flags
            dnl by setting any unset environment flag variables
            if test "x${CFLAGS+set}" != "xset"; then
                CFLAGS=""
            fi
            if test "x${CXXFLAGS+set}" != "xset"; then
                CXXFLAGS=""
            fi
            if test "x${FFLAGS+set}" != "xset"; then
                FFLAGS=""
            fi
            if test "x${FCFLAGS+set}" != "xset"; then
                FCFLAGS=""
            fi
            if test "x${OBJCFLAGS+set}" != "xset"; then
                OBJCFLAGS=""
            fi
        fi
        dnl assert.h is a NOP if NDEBUG is defined, so define it.
        AC_DEFINE(NDEBUG,,[Define if debugging is disabled])
    fi
    ax_enable_debug=$enable_debug
])
