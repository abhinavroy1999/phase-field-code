#suffixes used

.SUFFIXES: .o .h .c 

#compilation options

COMPOPS = -g -Wall 

#linking options

LINKOPS = -lgsl -lgslcblas -lfftw3 -lm

#list of object files

objects = CahnHilliard.o testcode.o evolution.o ps_file.o

#list of header files

headers = stdio.h stdlib.h math.h complex.h gsl_rng.h gsl_math.h fftw3.h headers.h 

#list of source codes

sources = CahnHilliard.c testcode.c evolution.c ps_file.c

#directory paths for the source and header files

vpath %.c source_code
vpath %.h ./header_file/
vpath %.h /usr/lib/gcc/i686-linux-gnu/5/include/ #gcc compiler directory. To be changed accordingly based on the machine.
vpath %.h /usr/include/                          #directory path for fftw3 header function. To be changed accordingly based on the directory path on the machine.
vpath %.h /usr/lib/i386-linux-gnu/               #lsgl and lsglcblas file directory. To be changed accordingly based on the machine.
vpath %.h /usr/include/gsl/                      #gsl header functions directory. To be changed accordingly based on the machine.

#actions on the source files	

CahnHilliard: $(objects) $(headers)
	gcc -o Cahn_Hilliard_evolution.out $(objects) $(LINKOPS) 

CahnHilliard.o: $(sources) $(headers)
	gcc -o $@ -c ./source_code/CahnHilliard.c $(COMPOPS)
evolution.o: evolution.c stdio.h stdlib.h math.h complex.h gsl_math.h fftw3.h headers.h
	gcc -c $(COMPOPS) $<
testcode.o: testcode.c stdlib.h stdio.h headers.h
	gcc -c $(COMPOPS) $<
ps_file.o: ps_file.c stdlib.h stdio.h
	gcc -c $(COMPOPS) $<

.PHONY : clean CLEAN
clean:
	rm -rf *.o
CLEAN:
	rm -rf Cahn_Hilliard_evolution.out


#end of the Makefile
