comdir = ${top_srcdir}/Common
randir = ${comdir}/RandomLib
tclap  = ${top_srcdir}/tclap

AM_CXXFLAGS = -Wall -O2 -funroll-loops -fexpensive-optimizations

AM_CXXFLAGS += -I${top_srcdir}
AM_CXXFLAGS += -I${comdir}

# real program stuff
bin_PROGRAMS = bprof

# Common sources
bprof_SOURCES =
bprof_SOURCES += ${comdir}/hashfunc.hh ${comdir}/hashfunc.cc
bprof_SOURCES += ${comdir}/label-map.hh ${comdir}/label-map.cc
bprof_SOURCES += ${comdir}/hash.hh ${comdir}/hash.cc
bprof_SOURCES += ${comdir}/newick.h ${comdir}/newick.c
bprof_SOURCES += ${comdir}/SCNode.h ${comdir}/SCNode.cc
bprof_SOURCES += ${comdir}/SCTree.h ${comdir}/SCTree.cc
bprof_SOURCES += ${comdir}/global.h ${comdir}/global.cc
bprof_SOURCES += ${comdir}/parsing.h ${comdir}/parsing.cc
bprof_SOURCES += ${comdir}/buildtree.h ${comdir}/buildtree.cc

# common sources: Mersenne Twister
bprof_SOURCES += ${comdir}/Random.cpp
bprof_SOURCES += ${randir}/Random.hpp
bprof_SOURCES += ${randir}/RandomAlgorithm.hpp
bprof_SOURCES += ${randir}/RandomCanonical.hpp
bprof_SOURCES += ${randir}/RandomEngine.hpp
bprof_SOURCES += ${randir}/RandomMixer.hpp
bprof_SOURCES += ${randir}/RandomPower2.hpp
bprof_SOURCES += ${randir}/RandomSeed.hpp
bprof_SOURCES += ${randir}/RandomType.hpp

#begin bprof sources
bprof_SOURCES += bprof.cc	
