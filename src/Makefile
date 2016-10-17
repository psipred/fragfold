######################################################################
#	Makefile for fragfold                           #
######################################################################

CC		= cc
LIBS		= -lm
CFLAGS	= -O3 -ffast-math -march=core2

all:	fragfold

install:
	/bin/cp fragfold ../bin

clean:
	/bin/rm -f fragfold
	/bin/rm -f ../bin/fragfold

fragfold:
	$(CC) $(CFLAGS) fragfold.c $(LIBS) -o fragfold
