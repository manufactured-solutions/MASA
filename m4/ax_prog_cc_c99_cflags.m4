# ===========================================================================
#         http://autoconf-archive.cryp.to/ax_prog_cc_c99_cflags.html
# ===========================================================================
#
# SYNOPSIS
#
#   Use AC_PROG_CC_C99 to determine how to enable C99 mode, appending
#   any required flags to CFLAGS.
#
#   AX_PROG_CC_C99_CFLAGS([ACTION-IF-AVAILABLE], [ACTION-IF-UNAVAILABLE])
#
# DESCRIPTION
#
#   Use AC_PROG_CC_C99 to determine how to enable C99 mode, but append
#   any required compiler flags to CFLAGS instead of modifying CC.
#   This allows other tools that use CFLAGS, e.g mpicc via ACX_MPI, to
#   use the flag.  CC is not modified.  The default actions are to
#   do nothing.
#
# LAST MODIFICATION
#
#   2009-06-08
#
# COPYLEFT
#
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PROG_CC_C99_CFLAGS],[
    AC_LANG_PUSH([C])

    ax_prog_cc_c99_cflags_CCSAVED=${CC}
    AC_PROG_CC_C99([$1],[$2])
    ax_prog_cc_c99_cflags_CFLAGSAPPEND=${CC#$ax_prog_cc_c99_cflags_CCSAVED}
    CFLAGS="${CFLAGS} ${ax_prog_cc_c99_cflags_CFLAGSAPPEND}"
    CC=${ax_prog_cc_c99_cflags_CCSAVED}

    AC_LANG_POP([C])
])
