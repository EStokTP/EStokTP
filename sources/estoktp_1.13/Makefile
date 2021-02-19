SHELL = /bin/sh
CC = gcc
F77 = gfortran -g  
#F77 = gfortran -g  -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Werror 
#F77 = gfortran -g  -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal  
OPTIONS_AB = -g -O3  -Wall  
OPTIONS =  -g  -Wall
CFLAGS=-O3 -std=c99
LIBS_CDF = -lnlopt -lm
LIBS1 = -lm 
GLIBS   =
# this include must specify where libraries and headers of nlopt are installed
LDFLAGS=-L/software/cinetica/local/lib -I/software/cinetica/local/include
#FLIBS = -lblas

# decomment/comment the following two lines to activate/deactivate compilation with NATST

#estoktp.x: estoktp.o g09fopt.o LineRead.o  mt.o croutines.o molprofopt.o zmat_routines.o 
#	$(F77)  -o ../../exe/estoktp.x  estoktp.o g09fopt.o LineRead.o croutines.o mt.o molprofopt.o zmat_routines.o  

estoktp.x: estoktp.o g09fopt.o LineRead.o  mt.o croutines.o molprofopt.o zmat_routines.o constrained_minimization_nlopt.o  target_01.o
	$(F77) $(LDFLAGS)  -o ../../exe/estoktp.x  estoktp.o g09fopt.o LineRead.o croutines.o mt.o molprofopt.o zmat_routines.o constrained_minimization_nlopt.o  target_01.o  $(LIBS_CDF)
estoktp.o: estoktp.f
	$(F77) -c   estoktp.f -o estoktp.o 
g09fopt.o: g09fopt.f
	$(F77)  -c   g09fopt.f -o g09fopt.o 
LineRead.o: LineRead.f
	$(F77)  -c  LineRead.f -o LineRead.o 
molprofopt.o: molprofopt.f
	$(F77) -c  molprofopt.f -o molprofopt.o 
zmat_routines.o: zmat_routines.f
	$(F77) -c zmat_routines.f -o zmat_routines.o
mt.o: mt.c
	$(CC) -c -fno-leading-underscore mt.c -o mt.o  
croutines.o: croutines.c
	$(CC) -c -fno-leading-underscore  croutines.c -o croutines.o  

# decomment/comment to activate/deactivate compilation with NATST
constrained_minimization_nlopt.o: constrained_minimization_nlopt.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c constrained_minimization_nlopt.c -o constrained_minimization_nlopt.o  $(LIBS_CDF)
target_01.o: target_01.f
	$(F77) -c -fno-leading-underscore target_01.f -o target_01.o
clean:
	rm -f *.o
	rm -f *~





