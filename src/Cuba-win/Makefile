HEADERS = cuba.h
LIB = libcuba.a

lib: $(LIB) $(HEADERS)

srcdir = .
common = $(srcdir)/src/common
vegas = $(srcdir)/src/vegas
suave = $(srcdir)/src/suave
divonne = $(srcdir)/src/divonne
cuhre = $(srcdir)/src/cuhre
demo = $(srcdir)/demo
tools = $(srcdir)/tools


COMMON_DEPS = $(common)/stddecl.h $(common)/ChiSquare.c $(common)/Erf.c
COMMON_CDEPS = $(COMMON_DEPS) $(common)/CSample.c $(common)/Parallel.c $(common)/sock.h
COMMON_TMDEPS = $(COMMON_DEPS) $(common)/MSample.c

RNG_DEPS = $(common)/Random.c $(common)/KorobovCoeff.c
R_DEPS = $(LIB)(r_helpers.o)

VEGAS_C = $(vegas)/Vegas.c
VEGAS_TM = $(vegas)/Vegas.tm
VEGAS_DEPS = $(RNG_DEPS) \
  $(vegas)/decl.h $(vegas)/common.c \
  $(vegas)/Grid.c $(vegas)/Integrate.c
VEGAS = $(VEGAS_C) $(VEGAS_TM) $(VEGAS_DEPS)

$(LIB)(Vegas.o): config.h $(VEGAS_C) $(VEGAS_DEPS) $(COMMON_CDEPS) $(R_DEPS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -I$(vegas) -DNOUNDERSCORE -c -o Vegas.o $(VEGAS_C)
	$(AR) $(ARFLAGS) $(LIB) Vegas.o
	$(RM) Vegas.o



SUAVE_C = $(suave)/Suave.c
SUAVE_TM = $(suave)/Suave.tm
SUAVE_DEPS = $(RNG_DEPS) \
  $(suave)/decl.h $(suave)/common.c \
  $(suave)/Grid.c $(suave)/Fluct.c $(suave)/Sample.c \
  $(suave)/Integrate.c
SUAVE = $(SUAVE_C) $(SUAVE_TM) $(SUAVE_DEPS)

$(LIB)(Suave.o): config.h $(SUAVE_C) $(SUAVE_DEPS) $(COMMON_CDEPS) $(R_DEPS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -I$(suave) -DNOUNDERSCORE -c -o Suave.o $(SUAVE_C)
	$(AR) $(ARFLAGS) $(LIB) Suave.o
	$(RM) Suave.o



DIVONNE_C = $(divonne)/Divonne.c
DIVONNE_TM = $(divonne)/Divonne.tm
DIVONNE_DEPS = $(RNG_DEPS) \
  $(divonne)/decl.h $(divonne)/common.c \
  $(divonne)/Rule.c $(divonne)/Sample.c \
  $(divonne)/FindMinimum.c $(divonne)/Iterate.c $(divonne)/Explore.c \
  $(divonne)/Split.c $(divonne)/Integrate.c
DIVONNE = $(DIVONNE_C) $(DIVONNE_TM) $(DIVONNE_DEPS)

$(LIB)(Divonne.o): config.h $(DIVONNE_C) $(DIVONNE_DEPS) $(COMMON_CDEPS) $(R_DEPS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -I$(divonne) -DNOUNDERSCORE -c -o Divonne.o $(DIVONNE_C)
	$(AR) $(ARFLAGS) $(LIB) Divonne.o
	$(RM) Divonne.o



CUHRE_C = $(cuhre)/Cuhre.c
CUHRE_TM = $(cuhre)/Cuhre.tm
CUHRE_DEPS = \
  $(cuhre)/decl.h $(cuhre)/common.c \
  $(cuhre)/Rule.c $(cuhre)/Integrate.c
CUHRE = $(CUHRE_C) $(CUHRE_TM) $(CUHRE_DEPS)

$(LIB)(Cuhre.o): config.h $(CUHRE_C) $(CUHRE_DEPS) $(COMMON_CDEPS) $(R_DEPS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -I$(cuhre) -DNOUNDERSCORE -c -o Cuhre.o $(CUHRE_C)
	$(AR) $(ARFLAGS) $(LIB) Cuhre.o
	$(RM) Cuhre.o

$(LIB)(r_helpers.o): $(common)/r_helpers.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DNOUNDERSCORE -c -o r_helpers.o $(common)/r_helpers.c
	$(AR) $(ARFLAGS) $(LIB) r_helpers.o
	$(RM) r_helpers.o

$(LIB)(Fork.o): $(common)/Fork.c $(common)/stddecl.h $(common)/sock.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -DNOUNDERSCORE -c -o Fork.o $(common)/Fork.c
	$(AR) $(ARFLAGS) $(LIB) Fork.o
	$(RM) Fork.o



$(LIB)(Global.o): $(common)/Global.c $(common)/stddecl.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -DNOUNDERSCORE -c -o Global.o $(common)/Global.c
	$(AR) $(ARFLAGS) $(LIB) Global.o
	$(RM) Global.o



$(LIB)(Data.o): $(common)/Data.c $(common)/stddecl.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o Data.o $(common)/Data.c
	$(AR) $(ARFLAGS) $(LIB) Data.o
	$(RM) Data.o

$(LIB): $(LIB)(Vegas.o)     \
	$(LIB)(Suave.o)     \
	$(LIB)(Divonne.o)   \
	$(LIB)(Cuhre.o)     \
	$(LIB)(Fork.o)      \
	$(LIB)(Global.o)    \
	$(LIB)(r_helpers.o) \
	$(LIB)(Data.o)
	-$(RANLIB) $(LIB)

clean:
	-$(RM) demo-c.out \
	  partview.moc partview.o makefile.partview \
	  cubal.h cubaq.h *.o */*.o */*.tm.c *.tm

distclean: clean
	-$(RM) $(LIB) demo-c demo-fortran \
	  config.h config.log config.status makefile

