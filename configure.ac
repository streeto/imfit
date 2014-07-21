AC_PREREQ([2.69])
AC_INIT([imfit], [1.0.2], [erwin@mpe.mpg.de])
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
AC_ARG_ENABLE([fftw_threads], [--enable-fftw-threads], [ENABLE_FFTW_THREADS=1], [ENABLE_FFTW_THREADS=0])
AM_CONDITIONAL([USE_FFTW_THREADS], [test "$ENABLE_FFTW_THREADS" -eq 1])
if test "$ENABLE_FFTW_THREADS" -eq 1; then
	# This is probably not the most portable way.
	fftw_LIBS="$fftw_LIBS -lfftw3_threads"
	AC_SUBST([fftw_LIBS])
fi

PKG_CHECK_MODULES([cfitsio], [cfitsio])

AC_ARG_ENABLE([gsl], [--enable-gsl], [ENABLE_GSL=1], [ENABLE_GSL=0])
PKG_CHECK_MODULES([gsl], [gsl], [HAVE_GSL=1], [HAVE_GSL=0])
AM_CONDITIONAL([USE_GSL], [test "$HAVE_GSL" -eq 1 && test "$ENABLE_GSL" -eq 1])

AC_ARG_ENABLE([nlopt], [--enable-nlopt], [ENABLE_NLOPT=1], [ENABLE_NLOPT=0])
PKG_CHECK_MODULES([nlopt], [nlopt], [HAVE_NLOPT=1], [HAVE_NLOPT=0])
AM_CONDITIONAL([USE_NLOPT], [test "$HAVE_NLOPT" -eq 1 && test "$ENABLE_NLOPT" -eq 1])

AC_CONFIG_FILES([config.h])
AC_CONFIG_FILES([imfit.pc Makefile])
AC_OUTPUT


