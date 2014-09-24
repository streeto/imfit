AC_PREREQ([2.69])
AC_INIT([imfit], [1.0.3], [erwin@mpe.mpg.de])
AC_CANONICAL_SYSTEM
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([1.10 -Wall])

AC_PROG_CXX
AC_PROG_CC
AM_PROG_AR
AC_LANG([C++])

LT_INIT

AC_OPENMP
AM_CONDITIONAL([USE_OPENMP], [test "$OPENMP_CXXFLAGS" != ""])

PKG_CHECK_MODULES([fftw], [fftw3])
AC_ARG_ENABLE([fftw_threads], [AS_HELP_STRING([--enable-fftw-threads],
		[enable threads in FFTW (default is no)])], 
	[ENABLE_FFTW_THREADS=1],
	[ENABLE_FFTW_THREADS=0])
AM_CONDITIONAL([USE_FFTW_THREADS], [test "$ENABLE_FFTW_THREADS" -eq 1])
if test "$ENABLE_FFTW_THREADS" -eq 1; then
	# This is probably not the most portable way.
	fftw_LIBS="$fftw_LIBS -lfftw3_threads"
	AC_SUBST([fftw_LIBS])
	echo Enabled FFTW threads.
fi

PKG_CHECK_MODULES([cfitsio], [cfitsio])

AC_ARG_WITH([gsl], [AS_HELP_STRING([--with-gsl],
		[enable some functions which depend on gsl (default is no)])],
	[WITH_GSL=yes],
	[WITH_GSL=no])
AS_IF([test x"$WITH_GSL" == xyes],
	[PKG_CHECK_MODULES([gsl], [gsl], [HAVE_GSL=yes], [HAVE_GSL=no])],
	[HAVE_GSL=no])
AM_CONDITIONAL([USE_GSL], [test "x$HAVE_GSL" == xyes])

AC_ARG_WITH([nlopt], [AS_HELP_STRING([--with-nlopt],
		[enable Nelder-Mead simplex algorithm using nlopt (default is no)])],
	[WITH_NLOPT=yes],
	[WITH_NLOPT=no])
AS_IF([test x"$WITH_NLOPT" == xyes],
	[PKG_CHECK_MODULES([nlopt], [nlopt], [HAVE_NLOPT=yes], [HAVE_NLOPT=no])],
	[HAVE_NLOPT=no])
AM_CONDITIONAL([USE_NLOPT], [test "x$HAVE_NLOPT" == xyes])

AC_CONFIG_FILES([config.h])
AC_CONFIG_FILES([imfit.pc Makefile])
AC_OUTPUT


