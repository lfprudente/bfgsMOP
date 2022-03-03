BIN ?= $(CURDIR)
DIR := $(notdir $(CURDIR))

# Check for the availability of the HSL library
HSL := $(shell if [ -e $(LIB)/libhsl.a ]; then echo true; fi)

ifeq ($(HSL),true)
  LIBHSL := $(LIB)/libhsl.a
else
  LIBHSL :=
endif

all: algencan

algencan: lib-ampl
	$(CC) $(CFLAGS) -I$(AMPL)/solvers amplwrapper.o amplsetp.o $(AMPL)/solvers/amplsolver.a libalgencan.a $(LIBHSL) -ldl -lm -lgfortran -fopenmp -o $@

lib-ampl: amplwrapper.o amplsetp.o
	cp $(LIB)/libalgencan.a ./libalgencan.a
	$(AR) d libalgencan.a setpoint.o

amplwrapper.o: amplwrapper.c
	$(CC) $(CFLAGS) -I$(AMPL)/solvers -c amplwrapper.c

amplsetp.o: amplsetp.f90
	$(FC) $(FFLAGS) -c -o $@ $<

install:
	mv -f algencan $(BIN)/$(DIR)

clean:
	rm -f amplwrapper.o amplsetp.o libalgencan.a

distclean: clean
	rm -f $(BIN)/$(DIR)/algencan

.PHONY: all install clean distclean
