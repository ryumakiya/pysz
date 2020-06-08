FC	= gfortran
FFLAGS	= -O3 -fopenmp 
OBJS	= rombint1.o

default: calc_ptilde

.SUFFIXES: .f90
.f90.o:
	$(FC) $(FFLAGS) -o $@ -c $< 

default: calc_ptilde

calc_ptilde: $(OBJS) calc_ptilde.o
	$(FC) $(FFLAGS) -o $@.exe $(OBJS) $@.o

clean:
	-rm -f *.o *.mod *.exe
tidy: clean
	-rm -f calc_ptilde
