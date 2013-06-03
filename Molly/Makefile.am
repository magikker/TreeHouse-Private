comdir = ${top_srcdir}/Common
randir = ${comdir}/RandomLib
tclap  = ${top_srcdir}/tclap

AM_CXXFLAGS = -Wall -O2 -funroll-loops -fexpensive-optimizations

AM_CXXFLAGS += -I${top_srcdir}
AM_CXXFLAGS += -I${comdir}

# real program stuff
bin_PROGRAMS = convertdb molly

#Begin ConvertDB sources here
convertdb_SOURCES =

#Begin Common Files
convertdb_SOURCES += ${comdir}/hashfunc.hh ${comdir}/hashfunc.cc
convertdb_SOURCES += ${comdir}/hash.hh ${comdir}/hash.cc
convertdb_SOURCES += ${comdir}/global.h ${comdir}/global.cc

# common sources: Mersenne Twister
convertdb_SOURCES += ${comdir}/Random.cpp
convertdb_SOURCES += ${randir}/Random.hpp
convertdb_SOURCES += ${randir}/RandomAlgorithm.hpp
convertdb_SOURCES += ${randir}/RandomCanonical.hpp
convertdb_SOURCES += ${randir}/RandomEngine.hpp
convertdb_SOURCES += ${randir}/RandomMixer.hpp
convertdb_SOURCES += ${randir}/RandomPower2.hpp
convertdb_SOURCES += ${randir}/RandomSeed.hpp
convertdb_SOURCES += ${randir}/RandomType.hpp

#Begin covertdb source files
convertdb_SOURCES += convertdb.cc

# Begin Molly sources here
molly_SOURCES = 

#Begin Common Files
molly_SOURCES += ${comdir}/hashfunc.hh ${comdir}/hashfunc.cc
molly_SOURCES += ${comdir}/hash.hh ${comdir}/hash.cc
molly_SOURCES += ${comdir}/global.h ${comdir}/global.cc

# common sources: Mersenne Twister
molly_SOURCES += ${comdir}/Random.cpp
molly_SOURCES += ${randir}/Random.hpp
molly_SOURCES += ${randir}/RandomAlgorithm.hpp
molly_SOURCES += ${randir}/RandomCanonical.hpp
molly_SOURCES += ${randir}/RandomEngine.hpp
molly_SOURCES += ${randir}/RandomMixer.hpp
molly_SOURCES += ${randir}/RandomPower2.hpp
molly_SOURCES += ${randir}/RandomSeed.hpp
molly_SOURCES += ${randir}/RandomType.hpp

#Begin molly source files
molly_SOURCES += molly.cc

#eof
