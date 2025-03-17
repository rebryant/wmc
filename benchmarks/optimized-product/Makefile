CC=gcc
OPT = -O0
CFLAGS = -g $(OPT) -Wno-nullability-completeness
LIB = -lgmp -lmpfr 
ALOCAL=/opt/homebrew
ALIB = $(ALOCAL)/lib/libgmp.a $(ALOCAL)/lib/libmpfr.a


all: product
arm: product-arm

product: product.c
	$(CC) $(CFLAGS) -o product product.c $(LIB)

product-arm: product.c
	$(CC) -I $(ALOCAL)/include $(CFLAGS) -o product-arm product.c $(ALIB)
