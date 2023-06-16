SHELL = /bin/sh
CC = gcc
OPTIONS_AB = -g -O3  -Wall  
OPTIONS =  -g  -Wall
CFLAGS  =
LIBS    =
LIBS1 = -lm 
GLIBS   =

PROJ_ROT:
	$(CC) -c $(OPTIONS_AB) main.c $(LIBS)
	$(CC) -c $(OPTIONS_AB) eigen2.c $(LIBS)
	$(CC) -c $(OPTIONS_AB) eigen.c $(LIBS)
	$(CC) -c $(OPTIONS_AB) nrutil.c $(LIBS)
	$(CC) -c $(OPTIONS_AB) spline.c $(LIBS)
	$(CC) -c $(OPTIONS_AB) invert_matrix3.c $(LIBS)

	$(CC)  $(OPTIONS_AB) main.o eigen2.o eigen.o nrutil.o spline.o invert_matrix3.o -o  ../../exe/RPHt.exe $(LIBS1)

clean:
	rm -f *.o
	rm -f *~


