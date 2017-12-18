# Process this file with autoconf to produce a configure script.
AC_INIT(craig, 2.1, axl.bernal@gmail.com)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/craigTrain.cpp)
AC_CONFIG_TESTDIR(tests)
AC_LANG(C++)

# New options
opt_default="yes"
AC_ARG_ENABLE(opt, [  --enable-opt=[no/yes] turn on optimization
                       [default=$opt_default]],, enable_opt=$opt_default)

if test "x$enable_opt" = "xyes"; then
CXXFLAGS="$CXXFLAGS -DNDEBUG -O3"
AC_MSG_RESULT(yes)
else

CXXFLAGS="$CXXFLAGS -g -DDEBUG"
AC_MSG_RESULT(no)
fi

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL

# Checks for libraries.
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_LIBTOOL

# Checks for header files.
AC_HEADER_STDC
AC_CONFIG_HEADERS(config)
AC_CHECK_HEADERS(ctype.h unistd.h)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

m4_include([acx_mpi.m4])
ACX_MPI

mpi_default="no"
AC_ARG_ENABLE(mpi, [  --enable-mpi=[no/yes]   turn on mpi support
                       [default=$mpi_default]],, enable_mpi=$mpi_default)

if test "x$enable_mpi" = "xyes"; then
CXXFLAGS="$CXXFLAGS -DWANT_MPI"
AC_MSG_RESULT(yes)
else
CXX="$acx_mpi_save_CXX"
AC_MSG_RESULT(no)
fi

# Building Makefile(s)
AM_INIT_AUTOMAKE([-Wall foreign])

# Set output files
AC_CONFIG_FILES([Makefile src/Makefile lless/Makefile regex/Makefile sparsevec/Makefile models/Makefile python/bin/Makefile python/lib/Makefile perl/bin/Makefile perl/lib/Makefile tests/Makefile doc/Makefile tests/atlocal datasets/Makefile])

AC_OUTPUT

