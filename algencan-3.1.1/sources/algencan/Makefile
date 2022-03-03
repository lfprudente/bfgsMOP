ALGENCAN := algencan.o newtonkkt.o auglag.o gencan.o betra.o moresor.o	\
   dogleg.o tnls.o spgls.o newtd.o extrap.o backtrack.o cgm.o calc.o	\
   happrox.o evalal.o fparam.o misc.o setpoint.o drand.o lss.o		\
   memalloc.o problvls.o problvlt.o problvlu.o problvlv.o probgiven.o	\
   modsydat.o modsvdhess.o modsvdgrad.o modrspace.o modprobdata.o	\
   modouttyp.o modminsq.o modmachconst.o modhpredat.o modhappdat.o	\
   moditetyp.o modalgparam.o modalgconst.o calgencan.o accpro.o

MA57 := $(shell if [ -e $(HSLSRC)/hsl_ma57_double.mod ]; then echo true; fi)
MA86 := $(shell if [ -e $(HSLSRC)/hsl_ma86_double.mod ]; then echo true; fi)
MA97 := $(shell if [ -e $(HSLSRC)/hsl_ma97_double.mod ]; then echo true; fi)

ifeq ($(MA57),true)
  HSL := lssma57.o
else
  HSL := dummy_lssma57.o
endif

ifeq ($(MA86),true)
  HSL := $(HSL) lssma86.o
else
  HSL := $(HSL) dummy_lssma86.o
endif

ifeq ($(MA97),true)
  HSL := $(HSL) lssma97.o
else
  HSL := $(HSL) dummy_lssma97.o
endif

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

all: lib

lib: $(ALGENCAN)
	$(AR) rcs libalgencan.a $(ALGENCAN) $(HSL)

install:
	mv -f libalgencan.a $(LIB)

algencan.o: algencan.f90 modalgconst.o modalgparam.o modmachconst.o	\
	    modouttyp.o problvls.o problvlt.o problvlu.o problvlv.o	\
	    probgiven.o $(HSL)

newtonkkt.o: newtonkkt.f90 modalgparam.o modmachconst.o modouttyp.o	\
	     modalgconst.o lss.o problvls.o probgiven.o modminsq.o

accpro.o: accpro.f90 modalgparam.o modmachconst.o modouttyp.o	\
	     modalgconst.o lss.o problvls.o probgiven.o

auglag.o: auglag.f90 modalgconst.o modalgparam.o modmachconst.o	\
	  modouttyp.o problvls.o probgiven.o

gencan.o: gencan.f90 modalgconst.o modalgparam.o modmachconst.o	\
	  modouttyp.o modprobdata.o modrspace.o moditetyp.o	\
	  modsydat.o problvls.o problvlv.o

betra.o: betra.f90 modalgconst.o modmachconst.o modouttyp.o moditetyp.o	\
	 problvlv.o probgiven.o

moresor.o: moresor.f90 modmachconst.o modouttyp.o lss.o

dogleg.o: dogleg.f90 modmachconst.o modouttyp.o 

tnls.o: tnls.f90 modalgconst.o modmachconst.o modouttyp.o problvlv.o

spgls.o: spgls.f90 modouttyp.o problvlv.o

newtd.o: newtd.f90 modalgparam.o modmachconst.o modouttyp.o moditetyp.o	\
	 modsvdgrad.o modrspace.o lss.o problvls.o probgiven.o

extrap.o: extrap.f90 modouttyp.o problvlv.o

backtrack.o: backtrack.f90 modalgconst.o modmachconst.o modouttyp.o \
	problvlv.o

cgm.o: cgm.f90 modalgconst.o modalgparam.o modmachconst.o modouttyp.o

calc.o: calc.f90 modrspace.o problvls.o

happrox.o: happrox.f90 modalgconst.o modalgparam.o modmachconst.o	\
	   modsvdgrad.o modsydat.o modhappdat.o modhpredat.o	\
	   moditetyp.o

evalal.o: evalal.f90 modalgparam.o modmachconst.o modouttyp.o	\
	  modsvdgrad.o problvls.o problvlv.o probgiven.o

fparam.o: fparam.f90 modalgparam.o modouttyp.o problvls.o problvlt.o	\
	  problvlu.o probgiven.o $(HSL)

misc.o: misc.f90 modmachconst.o

setpoint.o: setpoint.f90

drand.o: drand.f90

lss.o: lss.f90 $(HSL)

dummy_lssma57.o: dummy_lssma57.f90

dummy_lssma86.o: dummy_lssma86.f90

dummy_lssma97.o: dummy_lssma97.f90

lssma57.o: lssma57.f90 modmachconst.o modouttyp.o
	$(FC) $(FFLAGS) $(OPENMPFLAG) -I$(HSLSRC) -c -o $@ $<

lssma86.o: lssma86.f90 modmachconst.o modouttyp.o
	$(FC) $(FFLAGS) $(OPENMPFLAG) -I$(HSLSRC) -c -o $@ $<

lssma97.o: lssma97.f90 modmachconst.o modouttyp.o
	$(FC) $(FFLAGS) $(OPENMPFLAG) -I$(HSLSRC) -c -o $@ $<

memalloc.o: memalloc.f90 modrspace.o modsydat.o modsvdgrad.o	\
	    modsvdhess.o modhappdat.o modhpredat.o modminsq.o

problvls.o: problvls.f90 problvlt.o problvlv.o probgiven.o	\
	    modalgparam.o modmachconst.o modouttyp.o

problvlt.o: problvlt.f90 problvlu.o problvlv.o probgiven.o modouttyp.o 

problvlu.o: problvlu.f90 problvlv.o modouttyp.o modprobdata.o

problvlv.o: problvlv.f90 probgiven.o modmachconst.o modouttyp.o	\
	    modprobdata.o modsvdhess.o modsvdgrad.o

probgiven.o: probgiven.f90

modsydat.o: modsydat.f90

modsvdhess.o: modsvdhess.f90

modsvdgrad.o: modsvdgrad.f90

modrspace.o: modrspace.f90

modprobdata.o: modprobdata.f90

modouttyp.o: modouttyp.f90

modminsq.o: modminsq.f90

modmachconst.o: modmachconst.f90

modhpredat.o: modhpredat.f90

modhappdat.o: modhappdat.f90

moditetyp.o: moditetyp.f90

modalgparam.o: modalgparam.f90

modalgconst.o: modalgconst.f90

calgencan.o: calgencan.f90

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(LIB)/libalgencan.a

.PHONY: all clean
