CC=gcc
CXX=g++
OPT=-O2
#OPT=-O0
CFLAGS=-g $(OPT) -Wno-nullability-completeness -I $(IDIR)
IDIR = ../../include
LDIR = ../../lib
CPPFLAGS=-g $(OPT) -Wno-nullability-completeness -std=c++11 -I $(IDIR)

MYLIBS =  $(LDIR)/wmc_arithmetic.a $(LDIR)/wmc_util.a 
LIBS = $(MYLIBS) -lz -lgmpxx -lgmp -lmpfr -lmpfi

# ARM specific things
LOCAL=/opt/homebrew
MYALIBS = $(LDIR)/wmc_arithmetic_arm.a $(LDIR)/wmc_util.a

ALIBS = $(MYALIBS) -lz $(LOCAL)/lib/libgmpxx.a $(LOCAL)/lib/libgmp.a $(LOCAL)/lib/libmpfr.a $(LOCAL)/lib/libmpfi.a
ACPPFLAGS=$(CPPFLAGS) -I $(LOCAL)/include
 
all:  nnfcount

arm: nnfcount-arm


cnf_info.o: cnf_info.hh cnf_info.cpp 
	$(CXX) $(CPPFLAGS) -c cnf_info.cpp

egraph.o: egraph.hh egraph.cpp $(IDIR)/erd.hh
	$(CXX) $(CPPFLAGS) -c egraph.cpp

nnfcount: nnfcount.cpp cnf_info.o counters.o egraph.o $(MYLIBS)
	$(CXX) $(CPPFLAGS) $(GINC) -o nnfcount nnfcount.cpp cnf_info.o counters.o egraph.o $(LIBS)

nnfcount-arm: nnfcount.cpp cnf_info.o counters.o egraph.cpp $(MYALIBS)
	$(CXX) $(ACPPFLAGS) -o nnfcount-arm nnfcount.cpp cnf_info.o counters.o egraph.cpp $(ALIBS)


.SUFFIXES: .c .cpp .o

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $<

.cpp.o:
	$(CXX) $(CPPFLAGS) $(INC) -c $<

clean:
	rm -f *.o *~
	rm -f nnfcount
	rm -rf *.dSYM

