CC=g++
CPPFLAGS= -Wall -O3
VFS_OBJS=vfs.o
LIBCRE=crelib/libcre.a

all: vfs

vfs: $(VFS_OBJS) $(LIBCRE)
	$(LINK.cpp) $^ -o $@

$(LIBCRE):crelib/*
	cd crelib && make
