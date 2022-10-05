#*******************************************#
#                                           #
#  Program: Makefile                        #
#  Version: 2.0                             #
#  By: Mette Olufsen                        #
#  Date: 14. Jan. 1997                      #
#                                           #
#  A makefile that ensures that all modules #
#  are linked together in the right order.  #
#*******************************************#

## JAM's Mac ##
#CXX=/usr/local/bin/g++-11
#FC=/usr/local/bin/gfortran-11

## Euclid-02 ##
CXX=/usr/bin/g++
FC=/usr/bin/gfortran

CXXFLAGS=-O2 -g -Wall -D_REENTRANT -fPIC
FLIBS=-lgfortran
FFLAGS=-O2 -g -Wall -fPIC

LIBS=$(FLIBS) -lm

LDFLAGS=-O2

OBJS1=tools.o main.o junction.o arteries.o
OBJS2=impedance_sub.o new_match.o f90_tools.o

MAIN=main

all: $(MAIN)

$(MAIN): $(OBJS1) $(OBJS2)
	@$(CXX) -o $(MAIN) $(LDFLAGS) $(OBJS1) $(OBJS2) $(LIBS)

main.o: main.C main.h
	@$(CXX) -c $(CXXFLAGS) main.C

junction.o: junction.C junction.h arteries.h tools.h
	@$(CXX) -c $(CXXFLAGS) junction.C

arteries.o: arteries.C arteries.h junction.h tools.h main.h nr3.h ludcmp.h qrdcmp.h roots_multidim.h
	@$(CXX) -c $(CXXFLAGS) arteries.C

tools.o: tools.C tools.h
	@$(CXX) -c $(CXXFLAGS) tools.C

new_match.o: new_match.f90 f90_tools.o
	@$(FC) -c $(FFLAGS) new_match.f90

f90_tools.o: f90_tools.f90
	@$(FC) -c $(FFLAGS) f90_tools.f90

impedance_sub.o: impedance_sub.f90 f90_tools.o new_match.o
	@$(FC) -c $(FFLAGS) impedance_sub.f90

clean:
	@-rm -f *.o *.mod *Admit* main *.tmp

veryclean: clean
	@-rm $(MAIN) *~ main*.2d

