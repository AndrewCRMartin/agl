EXE=../agl
OFILES=agl.o findfields.o whereami/whereami.o
HFILES=agl.h findfields.h whereami/whereami.h

LIBDIR=$(HOME)/lib
INCLUDE=$(HOME)/include
LIBS=-lbiop -lgen -lm -lxml2 
#LIBS=-lbiop -lgen -lm -lxml2 -lasan

#CFLAGS=-g
CC=cc $(CFLAGS) -L$(LIBDIR) -I$(INCLUDE)
#CC=cc $(CFLAGS) -L$(LIBDIR) -I$(INCLUDE) -fsanitize=address

#CFLAGS=-O3
#CC=cc $(CFLAGS) -L$(LIBDIR) -I$(INCLUDE)

$(EXE) : $(OFILES)
	$(CC) -o $@ $(OFILES) $(LIBS)

agl.o : agl.c $(HFILES)
	$(CC) -c -o $@ $<

findfields.o : findfields.c $(HFILES)
	$(CC) -c -o $@ $<

whereami/whereami.o : whereami/whereami.c whereami/whereami.h
	$(CC) -c -o $@ $<

clean :
	\rm -f $(OFILES)

distclean : clean
	\rm -f $(EXE)
