SHELL := /bin/bash

CC=gcc
TESTSRC=force.c velverlet.c utilities.c ljmd_test.c
VPATH=../src
TESTOBJ=$(TESTSRC:%.c=%.o)

default: test
test: ljmd_test.x

ljmd_test.x: $(TESTOBJ)
	$(CC) -o $@ $^ -lm

%.o: ../src/%.c
	$(CC) -c $<
clean:
	rm -f *.o *.x
test: test
	./ljmd_test.x

