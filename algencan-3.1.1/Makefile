# Set variables FC and CC with the Fortran and C/C++ compilers of your choice.
# g77 and gfortran are already tested valid options for FC. gcc and g++ are
# also already tested valid options for CC. Setting CC is not necessary to run
# the stand-alone Fortran version of ALGENCAN. Leave it blank if this is your
# case.
AR := ar
FC := gfortran
CC := gcc-4.9
CP := g++-4.9

FFLAGS   := -O3
CFLAGS   := -O3

# Set variable ALGENCAN with the absolute path of your ALGENCAN installation
# directory. The value shown expands to the path of current working directory.
ALGENCAN := $(CURDIR)

BIN      := $(ALGENCAN)/bin
LIB      := $(ALGENCAN)/lib
SRC      := $(ALGENCAN)/sources
ALGSRC   := $(SRC)/algencan
HSLSRC   := $(SRC)/hsl
INTERSRC := $(SRC)/interfaces
INTERFCS := $(notdir $(wildcard $(INTERSRC)/*))

# Set the variables below with the absolute paths of your tools. The
# values shown are mere examples. Leave them blank if you feel that
# you will not use them. None of them is mandatory to use the
# stand-alone Fortran version of ALGENCAN.
AMPL     := $(HOME)/ampl

# Set the variables below to 1 if the source files of the BLAS and
# LAPACK dependencies of the HSL subroutines were placed within folder
# $HSLSRC. Otherwise, set it to 0. In the latter case, it is assumed
# that the BLAS and LAPACK libraries exist and that they will be used
# in the linking process.
BLAS_LAPACK := 1

# If you are using HSL_MA86 and/or HSL_MA97 and you want to use
# OpenMP, set the OpenMP flag of your compiler. Otherwise, leave it
# blank.
OPENMPFLAG := -fopenmp
# OPENMPFLAG :=

# CUTEst Directories
MASTSIF  := $(ALGENCAN)/myfolder/testes
#MASTSIF  := $(HOME)/CUTEst/sif
#MASTSIF  := $(HOME)/CUTEst/sif/netlib
#MASTSIF  := $(HOME)/CUTEr/MastSIF/marosmeszaros
ARCHDEFS  := $(HOME)/CUTEst/archdefs
SIFDECODE := $(HOME)/CUTEst/sifdecode
CUTEST    := $(HOME)/CUTEst/cutest
MYARCH    := mac64.osx.gfo

# Stop your modifications here.

export

ifeq ($(BLAS_LAPACK),1)
  HSL_MAKEFILE := Makefile.1
else
  HSL_MAKEFILE := Makefile.2
endif

all: algencan

algencan: hsl
	$(MAKE) -C $(ALGSRC) all install

algencan-%: algencan
	$(MAKE) -C $(INTERSRC)/$* all install

hsl:
	$(MAKE) -f $(HSL_MAKEFILE) -C $(HSLSRC)

clean:
	$(MAKE) -C $(ALGSRC) clean
	$(MAKE) -f $(HSL_MAKEFILE) -C $(HSLSRC) clean
	$(foreach i,$(INTERFCS),$(MAKE) -C $(INTERSRC)/$(i) clean;)

distclean: clean
	$(MAKE) -f $(HSL_MAKEFILE) -C $(HSLSRC) distclean
	$(foreach i,$(INTERFCS),$(MAKE) -C $(INTERSRC)/$(i) distclean;)

.PHONY: all clean distclean algencan*
