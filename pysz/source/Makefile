FC	= gfortran
FFLAGS	= -cpp -shared -fPIC -O3 -fopenmp -Wall
# FFLAGS	= -cpp -shared -fPIC -O3
OBJS	= cosmo.o global_var.o linearpk_z.o hunt.o spline.o sigma_z.o chebyshev.o dverk.o qromb.o rombint.o rombint1.o da.o input_ptilde.o mod_ptilde.o mf_module.o mf_tinker.o sz2mrs_functions.o cisi.o qgaus2.o qgaus2_n20.o

.SUFFIXES: .f90
.f90.o:
	$(FC) $(FFLAGS) -c $< 

default: calc_cl

calc_cl: $(OBJS) calc_cl.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)

clean:
	-rm -f *.o *.mod *~
tidy: clean
	-rm -f calc_cl
