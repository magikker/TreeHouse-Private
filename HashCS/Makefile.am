comdir = ${top_srcdir}/Common
randir = ${comdir}/RandomLib
tclap  = ${top_srcdir}/tclap

AM_CXXFLAGS = -Wall -O2 -funroll-loops -fexpensive-optimizations

AM_CXXFLAGS += -I${top_srcdir}
AM_CXXFLAGS += -I${comdir}

# real program stuff
bin_PROGRAMS = hashcs

#begin hashcs sources
hashcs_SOURCES =

#begin common sources
hashcs_SOURCES += ${comdir}/label-map.hh ${comdir}/label-map.cc
hashcs_SOURCES += ${comdir}/hashfunc.hh ${comdir}/hashfunc.cc
hashcs_SOURCES += ${comdir}/hash.hh ${comdir}/hash.cc
hashcs_SOURCES += ${comdir}/SCNode.h ${comdir}/SCNode.cc
hashcs_SOURCES += ${comdir}/SCTree.h ${comdir}/SCTree.cc

#begin Common Newick parser
hashcs_SOURCES += ${comdir}/newick.h ${comdir}/newick.c

#begin Common Mersenne Twister
hashcs_SOURCES += ${comdir}/Random.cpp ${randir}/Random.hpp
hashcs_SOURCES += ${randir}/RandomAlgorithm.hpp
hashcs_SOURCES += ${randir}/RandomCanonical.hpp
hashcs_SOURCES += ${randir}/RandomEngine.hpp
hashcs_SOURCES += ${randir}/RandomMixer.hpp
hashcs_SOURCES += ${randir}/RandomPower2.hpp
hashcs_SOURCES += ${randir}/RandomSeed.hpp
hashcs_SOURCES += ${randir}/RandomType.hpp

#begin Common other
hashcs_SOURCES += ${comdir}/global.h ${comdir}/global.cc
hashcs_SOURCES += ${comdir}/buildtree.h ${comdir}/buildtree.cc
hashcs_SOURCES += ${comdir}/parsing.h ${comdir}/parsing.cc

hashcs_SOURCES += hashmap.hh hashmap.cc
hashcs_SOURCES += hashcs.cc	
