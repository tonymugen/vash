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
LIBOUT = libvash.a
# LD blocks binary
LDBLK = ldblocks

CXXFLAGS = -O3 -march=native -std=gnu++14 -pthread
VASHOBJ = gvarHash.o random.o

all :  $(LDBLK) $(LIBOUT)
.PHONY : all

install : $(LIBOUT)
	-cp $(LIBOUT) $(INSTALLDIR)/lib
	-cp $(SRCDIR)/*.hpp $(INSTALLDIR)/include
	-cp $(RANDIR)/random.hpp $(INSTALLDIR)/include
	-cp $(LDBLK) $(INSTALLDIR)/bin
.PHONY : install

$(LDBLK) : ldblocks.cpp $(VASHOBJ)
	$(CXX) ldblocks.cpp $(VASHOBJ) -o $(LDBLK) $(CXXFLAGS)

$(LIBOUT) : $(VASHOBJ)
	ar crv $(LIBOUT) $(VASHOBJ)

random.o : $(RANDIR)/random.cpp $(RANDIR)/random.hpp
	$(CXX) -c $(RANDIR)/random.cpp $(CXXFLAGS)

gvarHash.o : $(SRCDIR)/gvarHash.hpp $(SRCDIR)/gvarHash.cpp random.o
	$(CXX) -c $(SRCDIR)/gvarHash.cpp $(CXXFLAGS)

.PHONY : clean
clean :
	-rm *.o $(LIBOUT) $(LDBLK)

