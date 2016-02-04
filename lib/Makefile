# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CDFLAGS=-Wall -fPIC  -g -O3 -fomit-frame-pointer -fno-fast-math
LDFLAGS=-shared
LDLIBS=-lm

# list of source files
SRC=ljmd.c force.c velverlet.c read_input.c write_output.c utilities.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)

clean:
	rm -f  *.o *.so *.x

lib: ljmd-serial.so

#linker rule
ljmd-serial.so: $(OBJ)
	$(CC) -o $@  $(LDFLAGS) $^ $(LDLIBS)

#compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CDFLAGS) $<

