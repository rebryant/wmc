CC=gcc
CXX=g++
OPT = -O2
#OPT = -O0
CFLAGS=-g $(OPT) -Wno-nullability-completeness $(DFLAGS)
CXXFLAGS=-std=c++11
INC = -I..
LDIR = ../../lib
IDIR = ../../include

LFILE = wmc_arithmetic.a

OFILES = q25.o analysis.o
IFILES = q25.h analysis.h

GLIB = -lz -lgmpxx -lgmp


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
	$(CC) $(CFLAGS) -o dtest dtest.c q25.o

atest: atest.c q25.o analysis.o
	$(CC) $(CFLAGS) -o atest atest.c q25.o analysis.o

find_product: find_product.c q25.o analysis.o
	$(CC) $(CFLAGS) -o find_product find_product.c q25.o analysis.o

tiny_product: tiny_product.c q25.o analysis.o
	$(CC) $(CFLAGS) -o tiny_product tiny_product.c q25.o analysis.o

random_product: random_product.c q25.o analysis.o
	$(CC) $(CFLAGS) -o random_product random_product.c q25.o analysis.o

mptest: mptest.c q25.o
	$(CC) $(CFLAGS) -o mptest mptest.c q25.o $(GLIB)

clean:
	 rm -f *.a *.o *~
