# Process this file with autoconf to produce a configure script 

AC_INIT([h1jet],[1.0.0],[a.lind@sussex.ac.uk])
AC_MSG_NOTICE([ ---------- Configuring H1jet])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

AC_MSG_NOTICE([ ---------- Checking for Fortran compiler])
AC_PROG_FC([gfortran ifort lf95 g95])
AC_FC_LIBRARY_LDFLAGS

dnl AC_MSG_NOTICE([ ---------- Checking for CHAPLIN library])
dnl AC_CHECK_LIB([chaplin], [HPL1], [], 
dnl [
dnl 	AC_MSG_ERROR([Could not find CHAPLIN library])
dnl ])

AC_MSG_NOTICE([ ---------- Checking for HOPPET and LHAPDF])
AC_MSG_CHECKING([for FCFLAGS and LIBS from HOPPET])
AC_SUBST([HOPPET_FCFLAGS], [`hoppet-config --fflags`])
AC_SUBST([HOPPET_LIBS], [`hoppet-config --libs`])
AC_MSG_RESULT([yes])
AC_MSG_CHECKING([for LIBS from LHAPDF])
AC_SUBST([LHAPDF_LIBS], [`lhapdf-config --libs`])
AC_MSG_RESULT([yes])


AC_MSG_CHECKING([for USERFILE])
AC_SUBST([USERFILE])
AS_IF([test ! -z "$USERFILE"], AC_MSG_RESULT([yes]), AC_MSG_RESULT([no]))
AM_CONDITIONAL([USERFILE_SET], [test ! -z "$USERFILE"])

AC_SUBST([USERFILE_OBJ], [$(patsubst %.f90,%.o,$USERFILE)])

AC_CONFIG_COMMANDS([mkdirbin], [$MKDIR_P bin])

AC_MSG_NOTICE([ ---------- Creating Makefiles])
dnl AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
	Makefile
	src/Makefile
])

AC_OUTPUT
AC_MSG_NOTICE([ ---------- make.])
