2CC=gcc
CXX=g++
#OPT = -O2
OPT = -O0
CFLAGS=-g $(OPT) -fPIC -Wno-nullability-completeness $(DFLAGS)
CXXFLAGS=-std=c++11
INC = -I..
LDIR = ../../lib
IDIR = ../../include

LFILE = wmc_arithmetic.a

OFILES = q25.o analysis.o
IFILES = q25.h analysis.h

all: $(LFILE) $(IFILES)
	cp -p $(LFILE) $(LDIR)
	cp -p $(IFILES) $(IDIR)

$(LFILE): $(OFILES)
	ar cr $(LFILE) $(OFILES)

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $<

.cxx.o:
	$(CXX) $(CFLAGS) $(INC) -c $<

dtest: dtest.c q25.o
	$(CXX) $(CFLAGS) -o dtest dtest.c q25.o

atest: atest.c q25.o analysis.o
	$(CXX) $(CFLAGS) -o atest atest.c q25.o analysis.o

clean:
	 rm -f *.a *.o *~
