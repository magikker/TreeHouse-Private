#                                               -*- Autoconf -*-    
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_CONFIG_MACRO_DIR([m4])
AC_INIT([habitat], [1.0], [sjm@cse.tamu.edu])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE
AC_PROG_LIBTOOL

# Setup compilers, use empty default flags
: ${CFLAGS=""}
AC_PROG_CC

: ${CXXFLAGS=""}
AC_PROG_CXX


# Checks for header files.
AC_LANG([C++])
AC_HEADER_STDC
AC_CHECK_HEADERS([stdint.h stdlib.h string.h sys/time.h unistd.h float.h sstream strstream])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor gettimeofday pow sqrt strchr strtoull])

# 
# Normally none of the AC_CONFIG_FILES lines below should be commented out. 
# Please make sure it's checked into git this way so you don't annout the next person!
#

#
# To build a distribution (make dist) for only a subset of the directories:
# a). Comment out the SUBDIRSAC_CONFIG_FILES lines for directories that will NOT be built
# c). Perform similar activities in Makefile.am (in the root directory only)
# d). Run autoreconf in the source directory
# e). make dist from the build directory (may be the same as the source directory)
#

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([TreeZip/Makefile])
#AC_CONFIG_FILES([Phlash/Makefile])
#AC_CONFIG_FILES([SetOp/Makefile])
#AC_CONFIG_FILES([BProf/Makefile])
#AC_CONFIG_FILES([Molly/Makefile])
#AC_CONFIG_FILES([HashCS/Makefile])
#AC_CONFIG_FILES([Example/Makefile])
#AC_CONFIG_FILES([ncl/Makefile])
AC_CONFIG_FILES([TreeHouse/Makefile])

AC_OUTPUT
