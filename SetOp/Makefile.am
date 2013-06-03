comdir = ${top_srcdir}/Common
randir = ${comdir}/RandomLib
tclap  = ${top_srcdir}/tclap

AM_CXXFLAGS = -Wall -O2 -funroll-loops -fexpensive-optimizations

AM_CXXFLAGS += -I${top_srcdir}
AM_CXXFLAGS += -I${comdir}

# real program stuff
bin_PROGRAMS = setop

#begin setop sources
setop_SOURCES =

#begin Common sources
setop_SOURCES += ${comdir}/hashfunc.hh ${comdir}/hashfunc.cc
setop_SOURCES += ${comdir}/label-map.hh ${comdir}/label-map.cc
setop_SOURCES += ${comdir}/hash.hh ${comdir}/hash.cc
setop_SOURCES += ${comdir}/newick.h ${comdir}/newick.c
setop_SOURCES += ${comdir}/parsing.h ${comdir}/parsing.cc
setop_SOURCES += ${comdir}/global.h ${comdir}/global.cc

# common sources: Mersenne Twister
setop_SOURCES += ${comdir}/Random.cpp
setop_SOURCES += ${randir}/Random.hpp
setop_SOURCES += ${randir}/RandomAlgorithm.hpp
setop_SOURCES += ${randir}/RandomCanonical.hpp
setop_SOURCES += ${randir}/RandomEngine.hpp
setop_SOURCES += ${randir}/RandomMixer.hpp
setop_SOURCES += ${randir}/RandomPower2.hpp
setop_SOURCES += ${randir}/RandomSeed.hpp
setop_SOURCES += ${randir}/RandomType.hpp

#begin setop sources
setop_SOURCES += setop.cc	
