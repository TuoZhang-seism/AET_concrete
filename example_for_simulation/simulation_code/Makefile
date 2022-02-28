# -xdebugformat=stabs
# -xdebugformat=dwarf


#CC = mpicc
CC = gcc
#FF = f77
#CFLAGS = -g -fast -O -xtarget=opteron -xarch=amd64  -I/usr/local/lam-7.1.2b28/include
#FFLAGS = -xarch=amd64
#CFLAGS = -ggdb -std=c99
#CFLAGS =  -O -std=c99



ALWAYS = edefinition.h scat_coeff_2D.h Makefile

ecount: ecount.o eutility.o mk_scat_pattern.o incidence.o mk_source_pattern.o scatter_2D.o surface.o layer.o $(ALWAYS)
	$(CC) $(CFLAGS) ecount.o eutility.o incidence.o mk_scat_pattern.o mk_source_pattern.o scatter_2D.o surface.o  layer.o -lm -o ecount

ecount.o: ecount.c  $(ALWAYS)
	$(CC) $(CFLAGS) -c ecount.c

eutility.o: eutility.c eutility.h $(ALWAYS)
	$(CC) $(CFLAGS) -c eutility.c

mk_scat_pattern.o: mk_scat_pattern.c mk_scat_pattern.h $(ALWAYS)
	$(CC) $(CFLAGS) -c mk_scat_pattern.c

mk_source_pattern.o: mk_source_pattern.c mk_source_pattern.h $(ALWAYS)
	$(CC) $(CFLAGS) -c mk_source_pattern.c

scatter_2D.o: scatter_2D.c scatter_2D.h $(ALWAYS)
	$(CC) $(CFLAGS) -c scatter_2D.c

surface.o: surface.c surface.h $(ALWAYS)
	$(CC) $(CFLAGS) -c surface.c

incidence.o: incidence.c incidence.h $(ALWAYS)
	$(CC) $(CFLAGS) -c incidence.c
	
old_incidence.o: incidence.f
	$(FF) $(FFLAGS) -c incidence.f

layer.o: layer.c layer.h $(ALWAYS)
	$(CC) $(CFLAGS) -c layer.c
	
exp_medium_eval: exp_medium_eval.o mk_scat_pattern.o eutility.o $(ALWAYS)
	$(CC) $(CFLAGS) exp_medium_eval.o mk_scat_pattern.o eutility.o -lm -o exp_medium_eval

exp_medium_eval.o: exp_medium_eval.c
	$(CC) $(CFLAGS) -c exp_medium_eval.c


misfit: misfit.o edefinition.h $(ALWAYS)
	$(CC) $(CFLAGS) misfit.o -lm -o misfit

misfit.o: misfit.c $(ALWAYS)
	$(CC) $(CFLAGS) -c misfit.c


clean:
	rm -f *.o
