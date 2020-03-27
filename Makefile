EXE=agl
OFILES=agl.o
LIBS=$(HOME)/lib
INCLUDE=$(HOME)/include
CFLAGS=-g
#CFLAGS=-O3
CC=cc $(CFLAGS) -L$(LIBS) -I$(INCLUDE)

$(EXE) : $(OFILES)
	$(CC) -o $@ $< -lbiop -lgen -lm -lxml2

.c.o :
	$(CC) -c -o $@ $<
