HEADERS = cuba.h
LIB = libcuba.a

lib: $(LIB) $(HEADERS)

srcdir = .
common = $(srcdir)/src/common
vegas = $(srcdir)/src/vegas
suave = $(srcdir)/src/suave
divonne = $(srcdir)/src/divonne
cuhre = $(srcdir)/src/cuhre

COMMON_DEPS = $(common)/stddecl.h $(common)/ChiSquare.c $(common)/Erf.c
COMMON_CDEPS = $(COMMON_DEPS) $(common)/CSample.c $(common)/Parallel.c $(common)/sock.h

VEGAS_C = $(vegas)/Vegas.c
VEGAS_DEPS = $(common)/Random.c \
  $(vegas)/decl.h $(vegas)/common.c \
  $(vegas)/Grid.c $(vegas)/Integrate.c

Vegas.o: config.h $(VEGAS_C) $(VEGAS_DEPS) $(COMMON_CDEPS) 
	$(CC) $(CFLAGS) -I$(vegas) -DNOUNDERSCORE -c -o Vegas.o $(VEGAS_C)
	$(AR) $(ARFLAGS) $(LIB) Vegas.o

SUAVE_C = $(suave)/Suave.c
SUAVE_DEPS = $(common)/Random.c \
  $(suave)/decl.h $(suave)/common.c \
  $(suave)/Grid.c $(suave)/Fluct.c $(suave)/Sample.c \
  $(suave)/Integrate.c

Suave.o: config.h $(SUAVE_C) $(SUAVE_DEPS) $(COMMON_CDEPS)
	$(CC) $(CFLAGS) -I$(suave) -DNOUNDERSCORE -c -o Suave.o $(SUAVE_C)
	$(AR) $(ARFLAGS) $(LIB) Suave.o

DIVONNE_C = $(divonne)/Divonne.c
DIVONNE_DEPS = $(common)/Random.c \
  $(divonne)/decl.h $(divonne)/common.c \
  $(divonne)/KorobovCoeff.c $(divonne)/Rule.c $(divonne)/Sample.c \
  $(divonne)/FindMinimum.c $(divonne)/Iterate.c $(divonne)/Explore.c \
  $(divonne)/Split.c $(divonne)/Integrate.c
DIVONNE = $(DIVONNE_C) $(DIVONNE_DEPS)

Divonne.o: config.h $(DIVONNE_C) $(DIVONNE_DEPS) $(COMMON_CDEPS)
	$(CC) $(CFLAGS) -I$(divonne) -DNOUNDERSCORE -c -o Divonne.o $(DIVONNE_C)
	$(AR) $(ARFLAGS) $(LIB) Divonne.o

CUHRE_C = $(cuhre)/Cuhre.c
CUHRE_DEPS = \
  $(cuhre)/decl.h $(cuhre)/common.c \
  $(cuhre)/Rule.c $(cuhre)/Integrate.c

Cuhre.o: config.h $(CUHRE_C) $(CUHRE_DEPS) $(COMMON_CDEPS)
	$(CC) $(CFLAGS) -I$(cuhre) -DNOUNDERSCORE -c -o Cuhre.o $(CUHRE_C)
	$(AR) $(ARFLAGS) $(LIB) Cuhre.o


Fork.o: $(common)/Fork.c $(common)/stddecl.h $(common)/sock.h
	$(CC) $(CFLAGS) -DNOUNDERSCORE -c -o Fork.o $(common)/Fork.c
	$(AR) $(ARFLAGS) $(LIB) Fork.o

Global.o: $(common)/Global.c $(common)/stddecl.h
	$(CC) $(CFLAGS) -DNOUNDERSCORE -c -o Global.o $(common)/Global.c
	$(AR) $(ARFLAGS) $(LIB) Global.o


Data.o: $(common)/Data.c $(common)/stddecl.h
	$(CC) $(CFLAGS) -c -o Data.o $(common)/Data.c
	$(AR) $(ARFLAGS) $(LIB) Data.o

$(LIB): $(LIB)(Vegas.o)     \
	$(LIB)(Suave.o)     \
	$(LIB)(Divonne.o)   \
	$(LIB)(Cuhre.o)     \
	$(LIB)(Fork.o)      \
	$(LIB)(Global.o)    \
	$(LIB)(Data.o)
	-$(RANLIB) $(LIB)

clean:
	-$(RM) demo-c.out \
	  partview.moc partview.o makefile.partview \
	  cubal.h cubaq.h *.o */*.o */*.tm.c *.tm

distclean: clean
	-$(RM) $(LIB) demo-c demo-fortran \
	  config.h config.log config.status makefile

