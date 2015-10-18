#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ(2.59)
AC_INIT([imcmc], [0.0.4], [yhxu@nao.cas.cn])
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE([disable])
#LT_PREREQ([2.2])
LT_INIT
AC_CONFIG_MACRO_DIR([m4])

AC_LANG_CPLUSPLUS
# Checks for programs.
AC_PROG_CC
AC_PROG_CXX([mpic++])
AC_PROG_CXXCPP


# Decide if this is a debug or optimized build.
AC_ARG_ENABLE(  debug,
                [AS_HELP_STRING([--enable-debug], [compile with debugging symbols (default no)])] )

# Set default compiler flags (if none were set by the user).
if test x"$enable_debug" = xyes; then
    # Debugging build
    AC_DEFINE([DEBUG], 1, [Enable debugging output])
    CFLAGS=${CFLAGS:--g -O2}
    CXXFLAGS=${CXXFLAGS:-$CFLAGS}
else
    enable_debug=no
    # Optimized build
    AX_CC_MAXOPT=no

    # This is broken for autoconf 2.64 and up, and I don't understand why.  It
    # gives errors about AC_PROG_CXX being expanded before required, even
    # though AX_CXX_MAXOPT is identical to AX_CC_MAXOPT...
    #AX_CXX_MAXOPT
    CXXFLAGS=${CXXFLAGS-$CFLAGS}
fi

# Decide whether or not to keep assert statements.
AC_ARG_ENABLE(
    assert,
    [AS_HELP_STRING([--disable-assert], [ignore assert statements])]
)
if test x"$enable_assert" = xno; then
    AC_DEFINE([NDEBUG], 1, [Ignore assert statements])
fi

# Decide whether or not to display verbose runtime messages.
AC_ARG_ENABLE(
    verbose,
    [AS_HELP_STRING([--enable-verbose], [display verbose runtime messages (default no)])]
)
if test x"$enable_verbose" = xyes; then
    AC_DEFINE([VERBOSE], [1], [Enable verbose output]),
else
    enable_verbose=no
fi


# Checks for common header files.
#AC_CHECK_HEADERS([limits.h pthread.h sys/stat.h unistd.h])
#AC_CHECK_HEADERS([limits.h sys/stat.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_HEADER_STDC
AC_C_CONST
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_C_VOLATILE
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STAT
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([memset strerror strtol pow sqrt])

# Check for GSL.
dnl GSL
AC_ARG_WITH(    gsllib,
                [  --with-gsllib=DIR       DIR=-LDIR, has libgsl and libgslcblas],
                gsllib="$withval/",gsllib= )

if test $gsllib; then
    AC_MSG_RESULT( Checking User specified GSL )
    if test $rpathset; then
        LDFLAGS="$LDFLAGS -L$gsllib -Wl,-rpath -Wl,$gsllib"
    else
        LDFLAGS="$LDFLAGS -L$gsllib"
    fi
    AC_CHECK_LIB(   gsl, main, ,
                    AC_MSG_ERROR(   "libgsl not in location specified by --with-gsllib: $gsllib"),
                    [ -lm -lgsl -L$gsllib ] )
    AC_CHECK_LIB(   gslcblas, main, ,
                    AC_MSG_ERROR("libgslcblas not in location specified by --with-gsllib: $gsllib"),
                    [ -lm -lgslcblas -L$gsllib ] )
else
    AC_MSG_RESULT(  Checking GSL )
    AC_CHECK_LIB(   gsl, gsl_rng_alloc,
                    ,
                    AC_MSG_ERROR("libgsl not found; use --with-gsllib"),
                    [ -lgsl -lgslcblas ] )
    AC_CHECK_LIB(   gslcblas, cblas_dsdot,
                    ,
                    AC_MSG_ERROR("libgslcblas not found; use --with-gsllib"),
                    [ -lgsl -lgslcblas ] )
fi

AC_ARG_WITH(    gslinc,
                [  --with-gslinc=VAL       VAL=-IVAL, has gsl include files ],
	            [gslinc="$withval/" CPPFLAGS="$CPPFLAGS -I$gslinc"],
                [gslinc=]   )

AC_MSG_RESULT(  Checking GSL Headers )
AC_CHECK_HEADERS(   [   gsl/gsl_errno.h \
                        gsl/gsl_math.h \
                        gsl/gsl_rng.h \
                        gsl/gsl_randist.h ],
                    ,
                    [   if test $gslinc; then
                            AC_MSG_ERROR("gsl headers not in location specified by --with-gslinc")
                        elif test $gsllib; then
                            AC_CHECK_HEADER([   $gsllib/../include/gsl/gsl_errno.h\
                                                $gsllib/../include/gsl_math.h \
                                                $gsllib/../include/gsl_rng.h \
                                                $gsllib/../include/gsl_randist.h ],
                                            [gslinc="$gsllib/../include/" CPPFLAGS="$CPPFLAGS -I$gslinc"],
                                            AC_MSG_ERROR("Can't find gsl includes in user specified location; use --with-gslinc"))
                        else
                            AC_MSG_ERROR("Can't find gsl includes; use --with-gslinc")
                        fi  ] )

#	Checking for mpi
dnl MPI
AC_ARG_WITH(    mpilib,
                [  --with-mpilib=DIR       DIR=-LDIR, has libmpi_cxx ],
                mpilib="$withval/",mpilib= )

if test $mpilib; then
    AC_MSG_RESULT( Checking User specified openmpi )
    if test $rpathset; then
        LDFLAGS="$LDFLAGS -L$mpilib -Wl,-rpath -Wl,$mpilib"
    else
        LDFLAGS="$LDFLAGS -L$mpilib"
    fi
    AC_CHECK_LIB(   mpi, main,
                    ,
                    AC_MSG_ERROR("libmpi not in location specified by --with-mpilib: $mpilib"),
                    [ -lmpi -L$mpilib ] )
    AC_CHECK_LIB(   mpi_cxx, main,
                    ,
                    AC_MSG_ERROR("libmpi_cxx not in location specified by --with-mpilib: $mpilib"),
                    [ -lmpi_cxx -L$mpilib ] )
else
    AC_MSG_RESULT(  Checking OpenMPI )
    AC_CHECK_LIB(   mpi, main, ,
                    AC_MSG_ERROR("libmpi not found; use --with-mpilib"), [-lmpi -lmpi_cxx] )
    AC_CHECK_LIB(   mpi_cxx, main, ,
                    AC_MSG_ERROR("libmpi_cxx not found; use --with-mpilib"), [-lmpi -lmpi_cxx] )
fi

AC_ARG_WITH(    mpiinc,
                [  --with-mpiinc=VAL       VAL=-IVAL, has mpi include files ],
        	    [ mpiinc="$withval/" CPPFLAGS="$CPPFLAGS -I$mpiinc" ],
                [ mpiinc= ] )

AC_MSG_RESULT(  Checking mpi Headers )
AC_CHECK_HEADERS(   [mpi.h],
                    ,
                    [   if test $mpiinc; then
                            AC_MSG_ERROR("mpi headers not in location specified by --with-mpiinc")
                        else
                            AC_MSG_ERROR("Can't find mpi includes; use --with-mpiinc")
                        fi  ] )

# Write Makefiles and configuration header.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES(    [ Makefile
                      src/Makefile
                      examples/Makefile  ]   )
AC_OUTPUT

# Summarize configuration.
cat << EOF
----------------------------------------------
    =======================================
    ==   imcmc configuration summary:    ==
    =======================================

Debug?              $enable_debug
Verbose?            $enable_verbose
#Python bindings?   $enable_python

Installation prefix = $prefix
CC                  = $CC
CXX                 = $CXX
CPPFLAGS            = $CPPFLAGS
CFLAGS              = $CFLAGS
CXXFLAGS            = $CXXFLAGS
LDFLAGS             = $LDFLAGS
LIBS                = $LIBS
OPENMP_CFLAGS       = $OPENMP_CFLAGS
SWIG                = $SWIG
PYTHON_CPPFLAGS     = $PYTHON_CPPFLAGS
PYTHON_LDFLAGS      = $PYTHON_LDFLAGS
----------------------------------------------
EOF
