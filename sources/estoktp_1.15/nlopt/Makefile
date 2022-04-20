
CPPFLAGS=-I/usr/local/include
LDFLAGS=-L/usr/local/lib64
CFLAGS=-O3 -std=c99 
F77FLAGS=-O3
F77=gfortran
LIBS=-lnlopt
MAIN = constrained_minimization_nlopt
PROBLEM ?= target_01
OBJS=$(MAIN).o $(PROBLEM).o


.PHONY : clean distclean all

all : $(MAIN)

$(PROBLEM).o : $(PROBLEM).f
	$(F77) $(F77FLAGS) -c -o $@ $<

$(MAIN).o : $(MAIN).c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

$(MAIN) : $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

clean :
	$(RM) $(OBJS)

distclean : clean
	$(RM) $(MAIN)
