### makefileBody.mk
# This file contains the static part of the Makefile. Ideally
# users will never have to make changes in this part and only
# make choices of what compiler, architecture, FFTW library,
# OpenMP and QDD options to use.
COMPILER = gfortran

#####################################################################
#                  QDD code status info injection                   #
#####################################################################
### See getInfo.mk header for info
#ifneq (,$(wildcard getInfo.mk))
#  include getInfo.mk
#else 
#  include Makefiles/getInfo.mk
#endif
######################################################################

#####################################################################
#                         QDD option selection                      #
#####################################################################
CODE_OPTIONS =
ifeq ($(strip $(FINDIFF)),YES)
  CODE_OPTIONS += -Dfindiff
endif
ifeq ($(strip $(NUMEROV)),YES)
  CODE_OPTIONS += -Dnumerov
endif
ifeq ($(strip $(RAREGAS)),YES)
  CODE_OPTIONS += -Draregas
endif
ifeq ($(strip $(EXTENDED)),YES)
  CODE_OPTIONS += -Dextended
endif
ifeq ($(strip $(FSIC)),YES)
  CODE_OPTIONS += -Dfsic
endif

#####################################################################
#                     Parallelisation activation                    #
#####################################################################
ifeq "$(strip $(OMP))" "YES"
  OMPADDL = -fopenmp
endif

#####################################################################
#            Compiler- and subroutine-dependent options             #
#####################################################################
  STATIC = -static
  OPT1 = $(OMPADDL) -w -Ofast -mfpmath=sse -mavx -fdefault-real-8 -fdefault-double-8
  OPT2 = $(OMPADDL) -w -Ofast -mfpmath=sse -mavx -fdefault-real-8 -fdefault-double-8
  OPT3 = $(OMPADDL) -w -g -fdefault-real-8 -fdefault-double-8
  ifeq "$(DEBUG)" "YES"
    OPT1 = $(OMPADDL) -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
    OPT2 = $(OMPADDL) -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
    OPT3 = $(OMPADDL) -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
  endif

#####################################################################
#  Final pre-processor and compiler flags for OMP/FFT and MKL   #
#####################################################################
### DERIVED FROM MAIN LIBRARY LOCATIONS AT TOP OF MAKEFILE
FFTW_LIBDIR = $(FFTW_PATH)/lib
FFTW_INCLUDE = $(FFTW_PATH)
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
  ifeq ($(strip $(DYNOMP)), YES) # This option will probably vanish in the near future.
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
    MKLFLAGS = -L$(MKL_LIBDIR)/intel64 -Wl,--no-as-needed -lmkl_gf_lp64
    ifeq ($(strip $(OMP)), YES) # OpenMP
      MKLFLAGS += -lmkl_gnu_thread -lgomp
    else
      MKLFLAGS += -lmkl_sequential
    endif
    MKLFLAGS += -lmkl_core -lpthread -lm -ldl
    LINKERFLAGS = $(MKLFLAGS)
  endif

  ifeq ($(strip $(FFT_TYPE)), FFTW)  # FFTW
  #  FFTWFLAGS += -L$(FFTW_LIBDIR) -lfftw3 -lm
  #  ifeq ($(strip $(OMP)), YES) # OpenMP
  #    FFTWFLAGS += -lfftw3_omp -lgomp
  #  endif
    FFTWFLAGS += -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
    ifeq ($(strip $(OMP)), YES) # OpenMP
      FFTWFLAGS += -lfftw3_omp -lgomp
    endif
    LINKERFLAGS = $(FFTWFLAGS)
  endif

else  # Static linking
  ifeq ($(strip $(FFT_TYPE)), MKL) # MKL
    MKLFLAGS = -Wl,--start-group
    MKLFLAGS += $(MKL_LIBDIR)/intel64/libmkl_gf_lp64.a
    ifeq ($(strip $(OMP)), YES) # OpenMP
      MKLFLAGS += $(MKL_LIBDIR)/intel64/libmkl_gnu_thread.a
    else
      MKLFLAGS += $(MKL_LIBDIR)/intel64/libmkl_sequential.a  
    endif
    MKLFLAGS += $(MKL_LIBDIR)/intel64/libmkl_core.a
    MKLFLAGS += -Wl,--end-group
    ifeq ($(strip $(OMP)), YES) # OpenMP
      MKLFLAGS += -Wl,-Bstatic -lgomp -Wl,-Bdynamic
    endif
    MKLFLAGS += -lpthread -lm -ldl
    LINKERFLAGS = $(MKLFLAGS)
  endif

  ifeq ($(strip $(FFT_TYPE)), FFTW)  # FFTW
    FFTWFLAGS += -lpthread -lm $(FFTW_LIBDIR)/libfftw3.a
    ifeq ($(strip $(OMP)), YES) # OpenMP
      FFTWFLAGS += $(FFTW_LIBDIR)/libfftw3_omp.a -Wl,-Bstatic -lgomp -Wl,-Bdynamic
    endif
    LINKERFLAGS = $(FFTWFLAGS)
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
	@echo "Parameters defined in Makefile (Solvers and extensions):"
	@echo "------------------------------------------------------------------------------"
	@echo "findiff       $(FINDIFF)"
	@echo "numerov       $(NUMEROV)"
	@echo "raregas       $(RAREGAS)"
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

# Check that MAC, MKL and gfortran are not used simultaniously:
ifeq ($(strip $(OS)), MAC)
  ifeq ($(strip $(FFT_TYPE)), MKL)
  $(error ERROR: Compiling and linking on macOS with GNU Fortran against Intel's \
    Math Kernel Library is not supported by Intel. Please compile with Intel Fortran, \
    or choose FFTW for Fast Fourier Transforms by setting 'FFT_TYPE = FFTW'.)
  endif
endif

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
	$(LINKER) $(COMPILERFLAGS1) -o $@ $(strip $(OBJS) $(LINKERFLAGS))
	mv -fv $(EXEC) ../../bin/

# Implicit and explicit compiler-independent build recipes for objects
ifneq (,$(wildcard mfTargets.mk))
  include mfTargets.mk
else 
  include Makefiles/mfTargets.mk
endif

# Explicit compiler-dependent recipes for objects
#QDD_info.o: QDD_info.F90
#	$(COMPILER) -ffree-line-length-none $(VERSION_FLAGS) -c $^
