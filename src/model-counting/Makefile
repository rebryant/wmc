CC=gcc
CXX=g++
OPT=-O2
#OPT=-O0
CFLAGS=-g $(OPT) -Wno-nullability-completeness -I $(IDIR)
IDIR = ../../include
LDIR = ../../lib
CPPFLAGS=-g $(OPT) -Wno-nullability-completeness -std=c++11 -I $(IDIR)

MYLIBS =  $(LDIR)/wmc_arithmetic.a $(LDIR)/wmc_util.a 
LIBS = $(MYLIBS) -lz

all:  nnfcount

cnf_info.o: cnf_info.hh cnf_info.cpp 
	$(CXX) $(CPPFLAGS) -c cnf_info.cpp

nnfcount: nnfcount.cpp cnf_info.o counters.o $(MYLIBS)
	$(CXX) $(CPPFLAGS) $(GINC) -o nnfcount nnfcount.cpp cnf_info.o counters.o $(LIBS)

.SUFFIXES: .c .cpp .o

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $<

.cpp.o:
	$(CXX) $(CPPFLAGS) $(INC) -c $<

clean:
	rm -f *.o *~
	rm -f nnfcount
	rm -rf *.dSYM

