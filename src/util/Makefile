2CC=gcc
CXX=g++
#OPT = -O2
OPT = -O0
CFLAGS=-g $(OPT) -fPIC -Wno-nullability-completeness $(DFLAGS)
CXXFLAGS=-std=c++11
INC = -I..
LDIR = ../../lib
IDIR = ../../include

LFILE = wmc_util.a

OFILES = report.o
IFILES = report.h

all: $(LFILE) $(IFILES)
	cp -p $(LFILE) $(LDIR)
	cp -p $(IFILES) $(IDIR)

$(LFILE): $(OFILES)
	ar cr $(LFILE) $(OFILES)

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $<

.cxx.o:
	$(CXX) $(CFLAGS) $(INC) -c $<

clean:
	 rm -f *.a *.o *~
