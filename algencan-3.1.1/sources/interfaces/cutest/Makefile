BIN ?= $(CURDIR)
DIR := $(notdir $(CURDIR))
SIF := ELFUN EXTER GROUP RANGE

# Check for the availability of the HSL library
HSL := $(shell if [ -e $(LIB)/libhsl.a ]; then echo true; fi)

# Decide which linking flags are necessary to use
LDFLAGS := -lalgencan -lcutest

ifeq ($(HSL),true)
  LDFLAGS := $(LDFLAGS) -lhsl

  ifeq ($(BLAS_LAPACK),0)
    LDFLAGS := $(LDFLAGS) -lblas -llapack
  endif
endif

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

all: sifdec algencan

sifdec:
	$(SIFDECODE)/bin/sifdecoder $(MASTSIF)/$(PROBNAME).SIF
	$(FC) -c $(FFLAGS) $(addsuffix .f,$(SIF))
	rm -f AUTOMAT.d

algencan: cutestwrapper.o $(addsuffix .o,$(SIF))
	$(FC) $^ -L$(CUTEST)/objects/$(MYARCH)/double -L$(LIB) -o $@ $(LDFLAGS) $(OPENMPFLAG)

cutestwrapper.o: cutestwrapper.f90
	$(FC) $(OPENMPFLAG) $(FFLAGS) -c -o $@ $<

install:
	mv -f OUTSDIF.d $(BIN)/$(DIR)
	mv -f algencan  $(BIN)/$(DIR)

clean:
	rm -f *.o
	rm -f $(addsuffix .f,$(SIF))

distclean: clean
	rm -f $(BIN)/$(DIR)/OUTSDIF.d
	rm -f $(BIN)/$(DIR)/algencan

.PHONY: all sifdec install clean distclean
