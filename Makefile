# -*- Makefile -*-
CC=gcc
default: serial 
.PHONY: test
serial:
	$(MAKE) $(MFLAGS) -C Obj-$@
python:
	$(MAKE) lib -C lib
test:
	$(MAKE) test -C test
clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean
	$(MAKE) -C test clean
	$(MAKE) -C lib clean

