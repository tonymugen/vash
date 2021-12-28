###############################################################
#
# Makefile for vash
#
###############################################################
CXX = g++
INSTALLDIR = /usr/local
SRCDIR = src
TSTDIR = tests
RANDIR = $(SRCDIR)/bayesicUtilities
# library
LIBOUT = vash.a

CXXFLAGS = -O3 -march=native -std=gnu++11 -pthread
VASHOBJ = gvarHash.o random.o

all : $(LIBOUT)
.PHONY : all

install : $(LIBOUT)
	-cp -v $(LIBOUT) $(INSTALLDIR)/lib
	-cp -v $(SRCDIR)/*.hpp $(INSTALLDIR)/include
	-cp -v $(RANDIR)/random.hpp $(INSTALLDIR)/include
.PHONY : install

$(LIBOUT) : $(VASHOBJ)
	ar crv $(LIBOUT) $(VASHOBJ)

random.o : $(RANDIR)/random.cpp $(RANDIR)/random.hpp
	$(CXX) -c $(RANDIR)/random.cpp $(CXXFLAGS)

gvarHash.o : $(SRCDIR)/gvarHash.hpp $(SRCDIR)/gvarHash.cpp random.o
	$(CXX) -c $(SRCDIR)/gvarHash.cpp $(CXXFLAGS)

.PHONY : clean
clean :
	-rm *.o $(LIBOUT)

