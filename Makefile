.SUFFIXES:
.SUFFIXES: .o .f90 .f

FC := gfortran
FFLAGS := -O3 -g

ALGLIB := $(shell if [ -e $(CURDIR)/libalgencan.a ]; then echo true; fi)

ifneq ($(ALGLIB),true)

all: algencan MOPsolver

algencan: 
	$(MAKE) -C $(CURDIR)/algencan-3.1.1 
	mv -f $(CURDIR)/algencan-3.1.1/lib/libalgencan.a $(CURDIR)
endif

OBJECTS = globals.o myproblem.o quadfun.o scalefactor.o evalfuns.o checkdF.o bfgs.o  BFGSCautions.o innersolver.o armijo.o morethuente.o lsvecopt.o MOPsolverBFGS.o MOPsolverStBFGSArmijo.o MOPsolverStBFGSWolfe.o main.o

MOPsolver: $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) -L$(CURDIR) -lalgencan

globals.o: globals.f90
	$(FC) -c $(FFLAGS) globals.f90

globals.mod: globals.o

myproblem.o: globals.mod myproblem.f90
	$(FC) -c $(FFLAGS) myproblem.f90

myproblem.mod: myproblem.o

quadfun.o: globals.mod myproblem.mod quadfun.f90
	$(FC) -c $(FFLAGS) quadfun.f90

scalefactor.o: globals.mod myproblem.mod scalefactor.f90
	$(FC) -c $(FFLAGS) scalefactor.f90

evalfuns.o: globals.mod myproblem.mod evalfuns.f90
	$(FC) -c $(FFLAGS) evalfuns.f90

checkdF.o: myproblem.mod checkdF.f90
	$(FC) -c $(FFLAGS) checkdF.f90

bfgs.o: bfgs.f90
	$(FC) -c $(FFLAGS) bfgs.f90

BFGSCautions.o: BFGSCautions.f90
	$(FC) -c $(FFLAGS) BFGSCautions.f90
	
innersolver.o: globals.mod myproblem.mod innersolver.f90
	$(FC) -c $(FFLAGS) innersolver.f90

armijo.o: armijo.f90
	$(FC) -c $(FFLAGS) armijo.f90

morethuente.o: morethuente.f
	$(FC) -c $(FFLAGS) morethuente.f

lsvecopt.o: globals.mod lsvecopt.f90
	$(FC) -c $(FFLAGS) lsvecopt.f90

MOPsolverBFGS.o: globals.mod MOPsolverBFGS.f90
	$(FC) -c $(FFLAGS) MOPsolverBFGS.f90

MOPsolverStBFGSArmijo.o: globals.mod MOPsolverStBFGSArmijo.f90
	$(FC) -c $(FFLAGS) MOPsolverStBFGSArmijo.f90

MOPsolverStBFGSWolfe.o: globals.mod MOPsolverStBFGSWolfe.f90
	$(FC) -c $(FFLAGS) MOPsolverStBFGSWolfe.f90

main.o: globals.mod myproblem.mod main.f90
	$(FC) -c $(FFLAGS) main.f90

CLEAN:
	rm -f *.mod *.o MOPsolver

.PHONY: CLEAN
