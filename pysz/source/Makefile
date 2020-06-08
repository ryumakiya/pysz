FC	= gfortran
FFLAGS	= -cpp -shared -fPIC -O3 -fopenmp -Wall
OBJS	= cosmo.o global_var.o linearpk_z.o hunt.o spline.o sigma_z.o chebyshev.o dverk.o qromb.o rombint.o rombint1.o da.o input_ptilde.o mf_module.o mf_tinker.o ptilde_module.o sz_functions.o cisi.o qgaus2.o qgaus2_n20.o

default: calc_cl calc_by calc_dydzdMh

calc_cl: $(OBJS) calc_cl.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o

calc_by: $(OBJS) calc_by.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o

calc_dydzdMh: $(OBJS) calc_dydzdMh.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o

%.mod: %.o
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^

clean:
	-rm -f *.o *.mod *~
tidy: clean
	-rm -f calc_cl
	-rm -f calc_by
	-rm -f calc_dydzdMh
