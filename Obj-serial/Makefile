# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -fPIC  -g -O3 -fomit-frame-pointer -fno-fast-math
LDFLAGS=-shared
LDLIBS=-lm

# list of source files
SRC=ljmd.c force.c velverlet.c read_input.c write_output.c utilities.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)

default: ../ljmd-serial.so

clean:
	rm -f *.mod *.o ../ljmd-serial.so

# linker rule
../ljmd-serial.so: $(OBJ)
	$(CC) -o $@  $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
ljmd.o: ljmd.c force.c velverlet.c read_input.c write_output.c utilities.c
