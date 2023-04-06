### makefileBody.mk
# This file contains the static part of the Makefile. Ideally
# users will never have to make changes in this part and only
# make choices of what compiler, architecture, FFTW library,
# OpenMP and QDD options to use.
COMPILER = ifort


#####################################################################
#                     Parallelisation activation                    #
#####################################################################
ifeq "$(strip $(OMP))" "YES"
  OMPADDL = -qopenmp
endif

#####################################################################
#            Compiler- and subroutine-dependent options             #
#####################################################################
  STATIC = -static
  OPT1 = $(OMPADDL) -fpp -w -xsse4.2 -Ofast -ip -no-prec-div -align all -autodouble
#  OPT2 = $(OMPADDL) -fpp -w -xW -O3 -ip -no-prec-div -align all -autodouble
  OPT3 = $(OMPADDL) -fpp -w -g -align all -autodouble
  ifeq "$(strip $(DEBUG))" "YES"
    OPT1 = $(OMPADDL) -pg -fpp -w -g -CB -traceback -align all -autodouble
#    OPT2 =$(OMPADDL) -pg -fpp -w -g -CB -traceback -align all -autodouble
    OPT3 =$(OMPADDL) -pg -fpp -w -g -align all -autodouble
  endif

#####################################################################
#  Final pre-processor and compiler flags for OMP/FFT and MKL   #
#####################################################################
### DERIVED FROM MAIN LIBRARY LOCATIONS AT TOP OF MAKEFILE
FFTW_LIBDIR = $(FFTW_PATH)/lib
FFTW_INCLUDE = $(FFTW_PATH)/include
MKL_LIBDIR = $(MKL_PATH)/lib
MKL_INCLUDE = $(MKL_PATH)/include/fftw

FFTWADD = -Dfftw_cpu
fftw_cpu_value := 1
fftwnomkl_value := 0
omp_value := 0
dynomp_value := 0

ifeq ($(strip $(FFT_TYPE)), FFTW)
  FFTWADD += -Dfftwnomkl
  fftwnomkl_value := 1
endif

ifeq ($(strip $(OMP)), YES)
  OMPADD = -Domp
  omp_value := 1
  ifeq ($(strip $(DYNOMP)), YES) 
    OMPADD += -Ddynomp
    dynomp_value := 1
  endif
endif

ifeq ($(strip $(OMP_DEBUG)), YES)
  OMPADD += -Domp_debug
endif

# Final compiler flags
ifeq ($(strip $(LINK_STATIC)), YES) # Static linking
  COMPILERFLAGS = $(STATIC)
endif
COMPILERFLAGS += $(ARCH) $(OMPADD) $(FFTWADD) -Dnompi $(CODE_OPTIONS)
COMPILERFLAGS1 = $(strip $(OPT1) $(COMPILERFLAGS))
COMPILERFLAGS2 = $(strip $(OPT2) $(COMPILERFLAGS))
COMPILERFLAGS3 = $(strip $(OPT3) $(COMPILERFLAGS))

#####################################################################
#                Linker configuration for FFT and MKL               #
#####################################################################
LINKER = $(COMPILER)

ifneq ($(strip $(LINK_STATIC)), YES) # Dynamic linking

  ifeq ($(strip $(FFT_TYPE)), MKL) # MKL
    MKLFLAGS = -lmkl_core -lmkl_intel_lp64 -lpthread -lm -ldl
    ifeq ($(strip $(OMP)), YES)
      MKLFLAGS += -lmkl_intel_thread -liomp5
    else
      MKLFLAGS += -lmkl_sequential
    endif
    LINKERFLAGS = $(MKLFLAGS)
  endif

  ifeq ($(strip $(FFT_TYPE)), FFTW)  # FFTW
    FFTWFLAGS += -L$(FFTW_LIBDIR) -lfftw3 -lm
    ifeq ($(strip $(OMP)), YES)
      FFTWFLAGS += -lfftw3_omp -liomp5
    endif
    LINKERFLAGS = $(FFTWFLAGS)
  endif

else  # Static linking

  ifeq ($(strip $(FFT_TYPE)), MKL) # MKL

    ifneq ($(strip $(OS)), MAC) # Linux
      MKLFLAGS = -Wl,--start-group
      MKLFLAGS += $(MKL_LIBDIR)/intel64/libmkl_intel_lp64.a
      ifeq ($(strip $(OMP)), YES)
        MKLFLAGS += $(MKL_LIBDIR)/intel64/libmkl_intel_thread.a
      else
        MKLFLAGS += $(MKL_LIBDIR)/intel64/libmkl_sequential.a
      endif
      MKLFLAGS += $(MKL_LIBDIR)/intel64/libmkl_core.a
      MKLFLAGS += -Wl,--end-group
      ifeq ($(strip $(OMP)), YES)
        MKLFLAGS += $(MKL_PATH)/../compiler/lib/intel64/libiomp5.a
      endif
      MKLFLAGS += -lpthread -lm -ldl
      LINKERFLAGS = $(MKLFLAGS)

    else  # macOS
      MKLFLAGS += $(MKL_LIBDIR)/libmkl_intel_lp64.a
      ifeq ($(strip $(OMP)), YES)
        MKLFLAGS += $(MKL_LIBDIR)/libmkl_intel_thread.a
      else
        MKLFLAGS += $(MKL_LIBDIR)/libmkl_sequential.a
      endif
      MKLFLAGS += $(MKL_LIBDIR)/libmkl_core.a
      ifeq ($(strip $(OMP)), YES)
        MKLFLAGS += $(MKL_PATH)/../compiler/lib/libiomp5.a
      endif
      MKLFLAGS += -lpthread -lm -ldl
      LINKERFLAGS = $(MKLFLAGS)
    endif
  endif

  ifeq ($(strip $(FFT_TYPE)), FFTW)  # FFTW

    ifneq ($(strip $(OS)), MAC) # Linux
      FFTWFLAGS += -lm $(FFTW_LIBDIR)/libfftw3.a
      ifeq ($(strip $(OMP)), YES)
        FFTWFLAGS += $(FFTW_LIBDIR)/libfftw3_omp.a
        FFTWFLAGS += $(MKL_PATH)/../compiler/lib/intel64/libiomp5.a
      endif
      LINKERFLAGS = $(FFTWFLAGS)

    else # macOS
      FFTWFLAGS += -lm $(FFTW_LIBDIR)/libfftw3.a
      ifeq ($(strip $(OMP)), YES)
        FFTWFLAGS += $(FFTW_LIBDIR)/libfftw3_omp.a
        FFTWFLAGS += $(MKL_PATH)/../compiler/lib/libiomp5.a
      endif
      LINKERFLAGS = $(FFTWFLAGS)
    endif
  endif
endif

#####################################################################
#                          Executable name                          #
#####################################################################
EXEC = qdd

#####################################################################
#                 Phony- and default targets                        #
#####################################################################
.PHONY: all clean distclean cleanall cleanobj cleanmod # NO TRUE FILE TARGET PREREQUISITS MAY
                       # APPEAR HERE, UNLESS YOU WANT THEM TO BE
                       # REBUILT EVERY TIME!

.DEFAULT_GOAL := all # In this case the default target is already pointing to 'all'
           # because it is setup to be the first target. However, if 'all'
           # weren't to be the first target, this statement will enforce it to
           # be still the default target.

all: $(EXEC)

clean: cleanall

cleanall: cleanobj cleanmod

distclean: cleanall
	@rm -vf ../../bin/$(EXEC)

cleanobj:
	@rm -vf *.o

cleanmod:
	@rm -vf *.mod

#####################################################################
#                             Checks                                #
#####################################################################
initial_checks:
	@echo ""
	@echo "##############################################################################"
	@echo "Lists of known compilation parameters."
	@echo "##############################################################################"
	@echo "Parameters defined in Makefile (Parallelization, includes and libs):"
	@echo "------------------------------------------------------------------------------"
	@echo "COMPILER      $(COMPILER)"
	@echo "OMP           $(OMP)"
	@echo "FFT_TYPE      $(FFT_TYPE)"
	@echo "DEBUG         $(DEBUG)"
	@echo "LINK_STATIC   $(LINK_STATIC)"
ifeq ($(strip $(FFT_TYPE)), MKL)
	@echo "MKL_PATH      $(MKL_PATH)"
	@echo "MKL_LIBDIR    $(MKL_LIBDIR)"
	@echo "MKL_INCLUDE   $(MKL_INCLUDE)"
else
	@echo "FFTW_PATH     $(FFTW_PATH)"
	@echo "FFTW_LIBDIR   $(FFTW_LIBDIR)"
	@echo "FFTW_INCLUDE  $(FFTW_INCLUDE)"
endif
	@echo "##############################################################################"
	@echo "-D flags enabled (1) by Makefile rules:"
	@echo "------------------------------------------------------------------------------"
	@echo "fftw_cpu      $(fftw_cpu_value)"
	@echo "fftwnomkl     $(fftwnomkl_value)"
	@echo "omp           $(omp_value)"
	@echo "dynomp        $(dynomp_value)"
	@echo ""  
	@echo "##############################################################################"
	@echo "Done with the initial checks, starting compilation..."
	@echo "##############################################################################"

#####################################################################
#                           Gather objects                          #
#####################################################################
ifneq (,$(wildcard mfObjects.mk))
  include mfObjects.mk
else 
  include Makefiles/mfObjects.mk
endif

#####################################################################
#              Building and linking: targets and recipes            #
#####################################################################
$(EXEC): initial_checks $(OBJS)
	@echo ""  
	@echo "##############################################################################"
	@echo Linking executable \'$@\'
	@echo "##############################################################################"
	$(LINKER) -o $@ $(strip $(OBJS) $(LINKERFLAGS))
	mv -fv $(EXEC) ../../bin/

# Implicit and explicit compiler-independent build recipes for objects
ifneq (,$(wildcard mfTargets.mk))
  include mfTargets.mk
else 
  include Makefiles/mfTargets.mk
endif

