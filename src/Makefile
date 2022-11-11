EXE=../agl
OFILES=agl.o findfields.o whereami/whereami.o
HFILES=agl.h findfields.h whereami/whereami.h

LIBS=$(HOME)/lib
INCLUDE=$(HOME)/include
CFLAGS=-g
#CC=cc $(CFLAGS) -L$(LIBS) -I$(INCLUDE) -fsanitize=address
#CFLAGS=-O3
CC=cc $(CFLAGS) -L$(LIBS) -I$(INCLUDE)

$(EXE) : $(OFILES)
	$(CC) -o $@ $(OFILES) -lbiop -lgen -lm -lxml2
#	$(CC) -o $@ $(OFILES) -lbiop -lgen -lm -lxml2 -lasan

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
