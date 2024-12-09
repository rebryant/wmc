CC=gcc
CXX=g++
OPT=-O2
#OPT=-O0
CFLAGS=-g $(OPT) -Wno-nullability-completeness -I $(IDIR)
IDIR = ../../include
LDIR = ../../lib
CPPFLAGS=-g $(OPT) -Wno-nullability-completeness -std=c++11 -I $(IDIR)
GDIR = ~/repos/d4pp/d4v2/3rdParty/glucose-3.0
GINC = -I $(GDIR)

MYLIBS =  $(LDIR)/wmc_arithmetic.a $(LDIR)/wmc_util.a $(GDIR)/glucose.a
LIBS = $(MYLIBS) -lz

all:  nnfcount

$(GDIR)/Solver.h:
	cd $(GDIR); make lib -f Makefile-compile

$(GDIR)/glucose.a:
	cd $(GDIR); make lib -f Makefile-compile

cnf.o: cnf.hh cnf.cpp $(GDIR)/Solver.h
	$(CXX) $(CPPFLAGS) $(GINC) -c cnf.cpp

nnfcount: nnfcount.cpp cnf.o counters.o $(MYLIBS)
	$(CXX) $(CPPFLAGS) $(GINC) -o nnfcount nnfcount.cpp cnf.o counters.o $(LIBS)

.SUFFIXES: .c .cpp .o

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $<

.cpp.o:
	$(CXX) $(CPPFLAGS) $(INC) -c $<

clean:
	rm -f *.o *~
	rm -f tpp nnfcount
	rm -rf *.dSYM

