CC = g++

comdir = ${top_srcdir}/Common
randir = ${comdir}/RandomLib
tclap  = ${top_srcdir}/tclap
treezipdir = ${top_srcdir}/TreeZip
#ncl  = ${top_srcdir}/ncl

#AM_LDFLAGS = ${ncl}/libncl.la
AM_CXXFLAGS = -std=c++11 -Wall -O2 -funroll-loops -fexpensive-optimizations -g -I.
AM_CXXFLAGS += -I${top_srcdir}
AM_CXXFLAGS += -I${comdir}
AM_CXXFLAGS += -I${treezipdir}
AM_CXXFLAGS += `pkg-config --cflags --libs libxml++-2.6`
#AM_CXXFLAGS += -I${ncl}

AM_CFLAGS =  $(AM_CXXFLAGS)

# real program stuff
bin_PROGRAMS = treehouse

#common sources
treehouse_SOURCES =  
treehouse_SOURCES += ${comdir}/hashfunc.hh ${comdir}/randomc.h ${comdir}/hashfunc.cc
treehouse_SOURCES += ${comdir}/label-map.hh ${comdir}/label-map.cc
treehouse_SOURCES += ${comdir}/hash.hh ${comdir}/hash.cc

# common sources: Newick parser
treehouse_SOURCES += ${comdir}/newick.h ${comdir}/newick.c

# common sources: Mersenne Twister
treehouse_SOURCES += ${comdir}/Random.cpp
treehouse_SOURCES += ${randir}/Random.hpp
treehouse_SOURCES += ${randir}/RandomAlgorithm.hpp
treehouse_SOURCES += ${randir}/RandomCanonical.hpp
treehouse_SOURCES += ${randir}/RandomEngine.hpp
treehouse_SOURCES += ${randir}/RandomMixer.hpp
treehouse_SOURCES += ${randir}/RandomPower2.hpp
treehouse_SOURCES += ${randir}/RandomSeed.hpp
treehouse_SOURCES += ${randir}/RandomType.hpp

# common sources: Tree Building
treehouse_SOURCES += ${comdir}/SCNode.h ${comdir}/SCNode.cc
treehouse_SOURCES += ${comdir}/SCTree.h ${comdir}/SCTree.cc
treehouse_SOURCES += ${comdir}/buildtree.h ${comdir}/buildtree.cc
treehouse_SOURCES += ${comdir}/global.h ${comdir}/global.cc
treehouse_SOURCES += ${comdir}/parsing.h ${comdir}/parsing.cc

# TreeZip Files
#treehouse_SOURCES += ${treezipdir}/normalizer.h ${treezipdir}/normalizer.cpp
treehouse_SOURCES += ${treezipdir}/nexus_compress.hh ${treezipdir}/nexus_compress.cc
treehouse_SOURCES += ${treezipdir}/nexus_decompress.hh ${treezipdir}/nexus_decompress.cc
treehouse_SOURCES += ${treezipdir}/compressfunc.h ${treezipdir}/compressfunc.cc
treehouse_SOURCES += ${treezipdir}/bmerge.h ${treezipdir}/bmerge.cc
treehouse_SOURCES += ${treezipdir}/tmerge.h ${treezipdir}/tmerge.cc

#begin treehouse sources

treehouse_SOURCES += pql.h
treehouse_SOURCES += pqlsymbol.h
treehouse_SOURCES += BipartitionTable.h
treehouse_SOURCES += Taxon.h
treehouse_SOURCES += pqlLexer.h pqlLexer.c
treehouse_SOURCES += pqlParser.h pqlParser.c
treehouse_SOURCES += Bipartition.h Bipartition.cc
treehouse_SOURCES += THGlobals.h THGlobals.cc
#treehouse_SOURCES += HashTableSearch.h HashTableSearch.cc
treehouse_SOURCES += AnalysisFunctions.h AnalysisFunctions.cc
treehouse_SOURCES += UtilityFunctions.h UtilityFunctions.cc
treehouse_SOURCES += SearchFunctions.h SearchFunctions.cc
treehouse_SOURCES += InterfaceFunctions.h InterfaceFunctions.cc
treehouse_SOURCES += UserFunctions.h UserFunctions.cc
treehouse_SOURCES += TreeHouse.cc
treehouse_LDADD = -lantlr3c -lreadline -lncurses
