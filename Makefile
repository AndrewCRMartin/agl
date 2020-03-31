EXE=agl
OFILES=agl.o findfields.o
LIBS=$(HOME)/lib
INCLUDE=$(HOME)/include
CFLAGS=-g
#CFLAGS=-O3
CC=cc $(CFLAGS) -L$(LIBS) -I$(INCLUDE)

$(EXE) : $(OFILES)
	$(CC) -o $@ $(OFILES) -lbiop -lgen -lm -lxml2

.c.o :
	$(CC) -c -o $@ $<

clean :
	\rm -f $(OFILES)

distclean : clean
	\rm -f $(EXE)
