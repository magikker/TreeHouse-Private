comdir = ${top_srcdir}/Common
randir = ${comdir}/RandomLib
tclap  = ${top_srcdir}/tclap

AM_CXXFLAGS = -Wall -O2 -funroll-loops -fexpensive-optimizations -g

AM_CXXFLAGS += -I${top_srcdir}
AM_CXXFLAGS += -I${comdir}

# real program stuff
bin_PROGRAMS = phlash

# Begin Phlash files here
phlash_SOURCES =

# Begin Common files here
phlash_SOURCES += ${comdir}/hashfunc.hh ${comdir}/hashfunc.cc
phlash_SOURCES += ${comdir}/label-map.hh ${comdir}/label-map.cc
phlash_SOURCES += ${comdir}/hash.hh ${comdir}/hash.cc
phlash_SOURCES += ${comdir}/newick.h ${comdir}/newick.c
phlash_SOURCES += ${comdir}/global.h ${comdir}/global.cc
phlash_SOURCES += ${comdir}/parsing.h ${comdir}/parsing.cc

# common sources: Mersenne Twister
phlash_SOURCES += ${comdir}/Random.cpp
phlash_SOURCES += ${randir}/Random.hpp
phlash_SOURCES += ${randir}/RandomAlgorithm.hpp
phlash_SOURCES += ${randir}/RandomCanonical.hpp
phlash_SOURCES += ${randir}/RandomEngine.hpp
phlash_SOURCES += ${randir}/RandomMixer.hpp
phlash_SOURCES += ${randir}/RandomPower2.hpp
phlash_SOURCES += ${randir}/RandomSeed.hpp
phlash_SOURCES += ${randir}/RandomType.hpp

# Begin Phlash files here
phlash_SOURCES += traverse.cc traverse.hh
phlash_SOURCES += sample.cc sample.hh
phlash_SOURCES += compute.cc compute.hh
phlash_SOURCES += phlash.cc
