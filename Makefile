# -*- Makefile -*-
SHELL=/bin/sh
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=Obj-serial/%.o)
############################################
TESTSRC=./src/force.c ./src/velverlet.c ./src/utilities.c ./test/ljmd_test.c
TESTOBJ=$(TESTSRC:%.c=%.o)

default: serial
test: ./test/ljmd_test.x

./test/ljmd_test.x: $(TESTOBJ)
	gcc -o $@ $^ -lm

force.o: ./src/force.c
	gcc -c $<

utilities.o: ./src/utilities.c
	gcc -c $<

velverlet.o: ./src/velverlet.c
	gcc -c $<

ljmd_test.o: ./test/ljmd_test.c
	gcc -c $<

test:
	./test/ljmd_test.x

serial:
	$(MAKE) $(MFLAGS) -C Obj-$@

clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean

