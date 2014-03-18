#
# Makefile for Project-- Reed-Solomon encoder/decoder
# 
# Revision History:
#  Jun 02, 2011    Nnoduka Eruchalu     Initial Revision
#  Mar 17, 2014    Nnoduka Eruchalu     Updated outputs to go to build dir.
#

CC     = g++
CFLAGS = -g -Wall -ansi -pedantic
ODIR = build

_OBJS = main.o reedSolomon.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

rs: $(OBJS)
	$(CC) $(OBJS) -o $(ODIR)/rs

$(ODIR)/reedSolomon.o: reedSolomon.cpp reedSolomon.h primitives.h
	$(CC) $(CFLAGS) -c -o $@ reedSolomon.cpp

$(ODIR)/main.o: main.cpp reedSolomon.h
	$(CC) $(CFLAGS) -c -o $@ main.cpp

clean:
	rm -f $(ODIR)/*.o $(ODIR)/rs



