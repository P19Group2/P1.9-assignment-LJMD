# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -g -O3 -fomit-frame-pointer -fno-fast-math -fopenmp -D_OPENMP
LDLIBS=-lm

# list of source files
SRC=ljmd.c force.c velverlet.c read_input.c write_output.c utilities.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)

default: ../ljmd-serial.x

clean:
	rm -f  *.o *.x

# linker rule
../ljmd-serial.x: $(OBJ)
	$(CC) -o $@  $(CFLAGS)  $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<
