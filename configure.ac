#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(cord, 0.3.4, sergey.astanin@polito.it)
AC_CONFIG_SRCDIR([amesh2d.h])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

## Checks for libraries.

AC_CHECK_LIB([m], [floor],,AC_MSG_ERROR([Math library not found]))
AC_CHECK_LIB([popt], [poptGetNextOpt],,AC_MSG_ERROR([popt library not found]))

# BLAS selection
AC_ARG_WITH([atlas],
	[AC_HELP_STRING([--with-atlas=[PATH]],
		[use BLAS implementation from ATLAS (default)])],
	[],
	[with_atlas=no])
AC_ARG_WITH([gotoblas],
	[AC_HELP_STRING([--with-gotoblas=[PATH]],
		[use GotoBLAS])],
	[],
	[with_gotoblas=no])
AC_ARG_WITH([fblas],
	[AC_HELP_STRING([--with-fblas=[PATH]],
		[use legacy BLAS implementation from netlib.org])],
	[],
	[with_fblas=no])
# use any BLAS if not sepcified
if test "x$with_atlas" == "xno" -a "x$with_fblas" == "xno" -a "x$with_gotoblas" == "xno" ; then
	with_atlas=try
# do not try GotoBLAS by default
#	with_gotoblas=try
	with_fblas=try
fi
# check for ATLAS libs
if test "x$with_atlas" != "xno"; then
	if test "x$with_atlas" != "xyes" -a "x$with_atlas" != "xtry" ; then
		LDFLAGS="${LDFLAGS} -L${with_atlas}"
	fi
	AC_CHECK_LIB([cerbla],[main],,
		AC_MSG_WARN([stand alone cerbla library not found]))
	AC_CHECK_LIB([atlas],[main],
		[
		 LIBS="-latlas ${LIBS}"
		 AC_DEFINE([HAVE_LIBATLAS],[1],[Define if you have libatlas])
		 AC_CHECK_LIB([f77blas],[main],,
		 	AC_MSG_ERROR([f77blas (ATLAS) library not found]))
		 AC_CHECK_LIB([cblas],[main],,
		 	AC_MSG_ERROR([cblas (ATLAS) library not found]))
		 with_atlas=yes
		],
		[
		 if test "x$with_atlas" != "xtry" ; then
			 AC_MSG_ERROR([ATLAS library not found])
		 fi
		 with_atlas=no
		])
fi
# check for GotoBLAS lib
if test "x$with_gotoblas" != "xno" -a "x$with_atlas" == "xno" ; then
	if test "x$with_gotoblas" != "xyes" -a "x$with_gotoblas" != "xtry" ; then
		LDFLAGS="${LDFLAGS} -L${with_gotoblas}"
	fi
	AC_CHECK_LIB([goto],[main],
		[
		 LIBS="-lgoto ${LIBS}"
		 AC_DEFINE([HAVE_LIBGOTO],[1],[Define if you have libgoto])
		 with_gotoblas=yes
		],
		[
		 if test "x$with_gotoblas" != "xtry" ; then
			AC_MSG_ERROR([GotoBLAS library not found])
		 fi
		 with_gotoblas=no
		])
	# in this case we need gslcblas too
	AC_CHECK_LIB([gslcblas],[main],,
		AC_MSG_ERROR([gslcblas (GSL) library not found]))
fi

# check for legacy BLAS libs
if test "x$with_fblas" != "xno"  -a "x$with_atlas" == "xno" -a "x$with_gotoblas" == "xno" ; then
	if test "x$with_fblas" != "xyes" -a "x$with_fblas" != "xtry" ; then
		LDFLAGS="${LDFLAGS} -L${with_fblas}"
	fi
	AC_F77_LIBRARY_LDFLAGS
	LIBS="${LIBS} ${FLIBS}"
	AC_SUBST(FLIBS)
	AC_CHECK_LIB([blas],[main],,
		AC_CHECK_LIB([fblas],[main],,
			AC_MSG_ERROR([legacy BLAS library not found])))
	# in this case we need gslcblas too
	AC_CHECK_LIB([gslcblas],[main],,
		AC_MSG_ERROR([gslcblas (GSL) library not found]))
fi

AC_CHECK_LIB([lapack],[main],
	     [
	      LIBS="-llapack ${LIBS}"
	      AC_DEFINE([HAVE_LIBLAPACK],[1],[Define if you have LAPACK])
	      with_lapack=yes],
	     [AC_MSG_WARN([LAPACK library not found])
	      with_lapack=no])

AC_CHECK_LIB([gsl],[gsl_block_alloc],,
	     AC_MSG_ERROR([Required GSL library not found]))

AC_ARG_WITH([umfpack],
	[AS_HELP_STRING([--with-umfpack[=PATH]],
		[use direct linear solvers from UMFPACK])],
	[],
	[with_umfpack=try])
if test "x$with_umfpack" != "xno" ; then
	if test "x$with_umfpack" != "xyes" ; then
		LDFLAGS="${LDFLAGS} -L${with_umfpack}"
	fi
	AC_CHECK_LIB([amd],[amd_l_valid],
	[
	 LIBS="-lamd ${LIBS}"
	 AC_DEFINE([HAVE_LIBAMD],[1],[Define if you have libamd])
	 AC_CHECK_LIB([umfpack],[umf_l_malloc],,
		[ # libumfpack not found
		if test "x$with_umfpack" != "xtry" ; then
		 AC_MSG_ERROR([--with-umfpack given, but libumfpack not found])
		fi
		with_umfpack=no
		])
	 with_umfpack=yes
	],
	[ # libamd not found
	if test "x$with_umfpack" != "xtry" ; then
	 AC_MSG_ERROR([--with-umfpack given, but libamd not found])
	fi
	with_umfpack=no
	])
fi

AC_CHECK_LIB([blitz],[main],,
	AC_MSG_ERROR([Required library Blitz++ not found]))

AC_ARG_WITH([hdf5],
	[AS_HELP_STRING([--with-hdf5],
		[use HDF5 library to save/load data])],
	[],
	[with_hdf5=try])
if test "x$with_hdf5" != "xno" ; then
	EXTRA_HDF5_LIBS=""
	AC_CHECK_LIB([sz], [main], EXTRA_HDF5_LIBS=" -lsz")
	AC_CHECK_LIB([z], [main], EXTRA_HDF5_LIBS="-lz${EXTRA_HDF5_LIBS}")
	AC_CHECK_LIB([hdf5],[main],
		     [with_hdf5=yes
		      LIBS="-lhdf5 ${EXTRA_HDF5_LIBS} ${LIBS}"
		      AC_DEFINE([HAVE_LIBHDF5],[1],[Define if you have libhdf5])
		      AC_CHECK_LIB([hdf5_cpp],[main],,
		       [AC_MSG_ERROR([hdf5_cpp not found, HDF5 is incomplete])])
		      AC_CHECK_LIB([hdf5_hl],[main],,
		       [AC_MSG_ERROR([hdf5_hl not found, HDF5 is incomplete])])
		     ],
		     [if test "x$with_hdf5" != "xtry" ; then
		      	AC_MSG_ERROR([--with-hdf5 was given, but test failed])
		      fi
		      with_hdf5=no
		     ])
fi

AC_SUBST(LIBS)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h],,AC_MSG_ERROR([Header files not found]))
if test "x$with_umfpack" != "xno" ; then
AC_CHECK_HEADERS([ufsparse/umfpack.h],,
	[AC_CHECK_HEADERS([umfpack/umfpack.h],,
		AC_MSG_ERROR([UMFPACK header not found]))])
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor sqrt])

# Build modes
AC_ARG_ENABLE([warnings],
	AS_HELP_STRING([--enable-warnings],
		[compile in pedantic mode with warnings]),
	[W_FLAGS="-Wall -Woverloaded-virtual -pedantic"],)
AC_ARG_ENABLE([more-warnings],
	AS_HELP_STRING([--enable-more-warnings],
		[compile in pedantic mode with lots of warnings]),
	[W_FLAGS="-Wall -Woverloaded-virtual -pedantic" EXTRA_W_FLAGS="-W -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fshort-enums -fno-common"],)
CXXFLAGS="${CXXFLAGS} ${W_FLAGS} ${EXTRA_W_FLAGS}"

INSTALL_PREFIX=${prefix}
CONFDIR=${prefix}/share/${PACKAGE}
AC_DEFINE_UNQUOTED([PREFIX],["$INSTALL_PREFIX"],[Installation prefix])
AC_DEFINE_UNQUOTED([CONFDIR],["$CONFDIR"],[Default config file location])

echo ; echo ; echo ;
echo "    Popt library: yes" ;
echo -n "    BLAS library: " ;
if test "x$with_atlas" == "xyes" ; then
	echo "ATLAS" ;
elif test "x$with_gotoblas" == "xyes" ; then
	echo "GotoBLAS" ;
else
	echo "Fortran BLAS" ;
fi
echo -n "    LAPACK library: ";
if test "x$with_lapack" == "xyes" ; then
	echo "yes" ;
else
	echo "no" ;
fi
echo "    GNU Scientific library: yes" ;
echo -n "    UMFPACK linear solvers: " ;
if test "x$with_umfpack" == "xyes" ; then
	echo "yes" ;
else
	echo "no" ;
fi
echo "    Blitz++ library: yes" ;
echo -n "    HDF5 data files: " ;
if test "x$with_hdf5" == "xyes" ; then
	echo "yes" ;
else
	echo "no" ;
fi
echo ; echo ; echo ;

AC_CONFIG_FILES([Makefile iniparser/Makefile])
AC_OUTPUT
