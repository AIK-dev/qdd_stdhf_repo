# Makefile: Build QDD with gfortran and Netlib only. No other libraries required.
# No MPI, no OpenMP, no Git QDD info injection

# LINK_STATIC: select static linkage of the binary
# Available options:
#   * YES
#   * NO
LINK_STATIC = NO

### DEBUG: enable debugging
# Available options:
#   * YES
#   * NO
DEBUG = NO

#####################################################################
#                         Code selection                            #
#####################################################################
# grid: FFT or finite difference (FFT options must be NO for finite differences)
FINDIFF = NO
NUMEROV = NO
#COUDOUB = YES
# Switch to extended model with polarizable raregas
RAREGAS = NO
# Include private sector
EXTENDED = NO
FSIC = NO

CODE_OPTIONS =
ifeq ($(FINDIFF),YES)
  CODE_OPTIONS += -Dfindiff
endif
ifeq ($(NUMEROV),YES)
  CODE_OPTIONS += -Dnumerov
endif
ifeq ($(COUDOUB),YES)
  CODE_OPTIONS += -Dcoudoub
endif
ifeq ($(RAREGAS),YES)
  CODE_OPTIONS += -Draregas
endif
ifeq ($(EXTENDED),YES)
  CODE_OPTIONS += -Dextended
endif
ifeq ($(FSIC),YES)
  CODE_OPTIONS += -Dfsic
endif

#####################################################################
#            Compiler- and subroutine-dependent options             #
#####################################################################
#OPT2: setting for the FFT package, needs forced double precision
#OPT3: setting for critical soubroutines which do not stand optimization
OPT1 = -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
OPT2 = -w -O3 -mfpmath=sse -fdefault-real-8 -fdefault-double-8
OPT3 = -w -g -fdefault-real-8 -fdefault-double-8
ifeq "$(DEBUG)" "YES"
	OPT1 = -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
	OPT2 = -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
	OPT3 = -pg -w -g -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcheck=bounds
endif

#####################################################################
#  Final pre-processor and compiler flags for MPI/OMP/FFT and MKL   #
#####################################################################
netlib_fft_value := 1

# Final compiler flags
COMPILERFLAGS1 = $(OPT1) $(strip -Dnetlib_fft -Dnompi) $(CODE_OPTIONS)
COMPILERFLAGS2 = $(OPT2) $(strip -Dnetlib_fft -Dnompi) $(CODE_OPTIONS)
COMPILERFLAGS3 = $(OPT3) $(strip -Dnetlib_fft -Dnompi) $(CODE_OPTIONS)

#####################################################################
#                Linker configuration for FFT and MKL               #
#####################################################################
COMPILER = gfortran
LINKER = $(COMPILER)

ifeq ($(LINK_STATIC), YES)
  LINKERFLAGS += $(STATIC)
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
	@echo "##############################################################################"
	@echo "Lists of known compilation parameters."
	@echo "##############################################################################"
	@echo "Parameters defined in Makefile (Solvers and extensions):"
	@echo "------------------------------------------------------------------------------"
	@echo "findiff       $(FINDIFF)"
	@echo "numerov       $(NUMEROV)"
	@echo "coudoub       $(COUDOUB)"
	@echo "raregas       $(RAREGAS)"
	@echo "##############################################################################"
	@echo "Parameters defined in Makefile (Parallelization, includes and libs):"
	@echo "------------------------------------------------------------------------------"
	@echo "DEBUG         $(DEBUG)"
	@echo "LINK_STATIC   $(LINK_STATIC)"
	@echo "##############################################################################"

# Check that variables which should be either YES or NO have proper values:
ifneq ($(DEBUG), $(filter $(DEBUG), YES NO))
  $(error ERROR: Unknown $(strip DEBUG value) in Makefile ($(DEBUG)))
endif

ifneq ($(LINK_STATIC), $(filter $(LINK_STATIC), YES NO))
  $(error ERROR: Unknown $(strip LINK_STATIC value) in Makefile ($(LINK_STATIC)))
endif

	@echo "Done with the initial checks, starting compilation..."
	@echo "#######################################################"

#####################################################################
#                  Targets and common dependencies                  #
#####################################################################
OBJINT = params.o main.o kinetic.o restart.o restartc.o init.o\
       static.o dynamic.o lda.o util.o abso_bc.o\
       pseudosoft.o pseudogoed.o ionmd.o forces.o\
       carlo.o localize.o localizer.o\
       sicnew.o sicnewc.o rho.o rhoc.o nonloc.o nonlocc.o\
       schmid.o zeroforce.o loc_mfield.o givens.o subgrids.o\
       parallele.o rta.o\
       HEeigensystem.o mini.o zdiag.o orthmat.o

ifeq ($(EXTENDED),YES)
  OBJINT += util_extended.o attachement.o V2ph.o
  DEPS_EXTENDED = util_extended.o attachement.o V2ph.o
else
  DEPS_EXTENDED =
endif

ifeq ($(FSIC),YES)
  OBJINT += 2stUT.o 2st_util.o 2stUTc.o
endif

ifeq ($(RAREGAS), YES)
  OBJINT +=  functions.o pot_substrate.o forces_substrate.o md_substrate.o\
    short.o image.o lattice.o util_substrate.o
endif

# List of modules:
MODLIST = params.mod kinetic.mod
OBJINT += coulsolv.o
MODLIST += coulsolv.mod

# Object file for FFT:
  FFT_OBJ = fftpack.o

ifeq ($(FSIC),YES)
  MODLIST += 2st_util.mod twostr.mod twost.mod
endif

OBJS = $(OBJINT) $(FFT_OBJ)

# Objects that depend on kinetic.mod and/or coulsolv.mod,
# with no other special requirements:
OBJS_KIN  = dynamic.o ionmd.o util.o zeroforce.o
OBJS_COUL = image.o lda.o loc_mfield.o
OBJS_KINCOUL = init.o main.o pseudosoft.o
DEPS_KINCOUL = kinetic.mod coulsolv.mod

# Dependencies for objects that depend on kinetic.o and/or coulsolv.o,
# with no other special requirements:
$(OBJS_KIN): kinetic.mod
$(OBJS_COUL): coulsolv.mod
$(OBJS_KINCOUL): $(DEPS_KINCOUL)

#####################################################################
#                   Compilation and linkage rules                   #
#####################################################################
$(EXEC): initial_checks $(OBJS)
	@echo Linking executable $@
	$(LINKER) -o $@ $(strip $(OBJS) $(LINKERFLAGS))
	mv -fv $(EXEC) ../../bin/

# Implicit rule for objects and modules together:
%.o %.mod: %.F90
	$(COMPILER) $(COMPILERFLAGS1) -c $<

# Additions for twostsic:
ifeq ($(FSIC),YES)
FSICDEPEND = 2st_util.o 2stUTc.o 2stUT.o

2st_util.o 2st_util.mod: fsic/2st_util.F90 kinetic.mod orthmat.o
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o 2st_util.o -c $<

2stUTc.o twost.mod: fsic/2stUT.F90 kinetic.mod twostr.mod
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o 2stUTc.o -c $<

2stUT.o twostr.mod: fsic/2stUT.F90 kinetic.mod 2st_util.mod
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o 2stUT.o -c $<
endif

# Explicit rules
# (complex dependencies, REAL/COMPLEX switch, non-default COMPILERFLAGS).

kinetic.o : kinetic.F90 fft.F90 findiff/findiff.F90 $(FFT_OBJ)

coulsolv.o: coulsolv.F90 falr.o coulex.o kinetic.mod $(FFT_OBJ)

static.o: static.F90 pseudosoft.F90 $(DEPS_KINCOUL) $(FSICDEPEND) $(DEPS_EXTENDED) 

# dynamic.o: dynamic.F90 pseudosoft.o twostr.mod twost.mod $(DEPS_KINCOUL) $(FSICDEPEND)
dynamic.o: dynamic.F90 pseudosoft.o $(DEPS_KINCOUL) $(FSICDEPEND) $(DEPS_EXTENDED) 

rho.o: rho.F90 $(DEPS_EXTENDED)
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o $@ -c $<

rhoc.o: rho.F90 kinetic.mod $(DEPS_EXTENDED)
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

localizer.o: localize.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o $@ -c $<

localize.o: localize.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

sicnew.o: sicnew.F90 $(DEPS_KINCOUL)
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o $@ -c $<

sicnewc.o: sicnew.F90 $(DEPS_KINCOUL)
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

nonloc.o: nonloc.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o $@ -c $<

nonlocc.o: nonloc.F90 kinetic.mod
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

fftpack.o: fftpack.F90 fftpack2.F90
	$(COMPILER) $(COMPILERFLAGS2) -c $<

main.o: main.F90 util.o orthmat.o rta.mod $(FSICDEPEND) $(DEPS_EXTENDED)
	$(COMPILER) $(COMPILERFLAGS3) -c $<

init.o: init.F90 rta.mod params.mod
	$(COMPILER) $(COMPILERFLAGS3) -c $<

givens.o: givens.F90
	$(COMPILER) $(COMPILERFLAGS3) -c $<

restart.o: restart.F90 kinetic.mod $(FSICDEPEND)
	$(COMPILER) $(COMPILERFLAGS3) -DREALSWITCH -o $@ -c $<

restartc.o: restart.F90 kinetic.mod $(FSICDEPEND)
	$(COMPILER) $(COMPILERFLAGS3) -DCOMPLEXSWITCH -o $@ -c $<

parallele.o: parallele.F90 kinetic.mod $(DEPS_EXTENDED)
	$(COMPILER) $(COMPILERFLAGS3) -c $<

rta.o: rta.F90 kinetic.mod params.mod
	$(COMPILER) $(COMPILERFLAGS3) -c $<

orthmat.o: orthmat.F90
	$(COMPILER) $(COMPILERFLAGS3) -c -cpp $<

zdiag.o: zdiag.f
	$(COMPILER) $(COMPILERFLAGS3) -c $<

mini.o: mini.f
	$(COMPILER) $(COMPILERFLAGS3) -c $<

HEeigensystem.o: HEeigensystem.f
	$(COMPILER) $(COMPILERFLAGS3) -c -cpp $<

ifeq ($(EXTENDED),YES)
attachement.o: extended/attachement.F90
	$(COMPILER) $(COMPILERFLAGS3) -c $<

util_extended.o: extended/util_extended.F90 params.o util.o
	$(COMPILER) $(COMPILERFLAGS3) -c $<

V2ph.o: extended/V2ph.f90 kinetic.mod params.mod
	$(COMPILER) $(COMPILERFLAGS1) -o $@ -c $<
endif

# QM/MM routines
ifeq ($(RAREGAS), YES)
forces_substrate.o: QMMM/forces_substrate.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

util_substrate.o: QMMM/util_substrate.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

functions.o: QMMM/functions.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

pot_substrate.o: QMMM/pot_substrate.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

md_substrate.o: QMMM/md_substrate.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

short.o: QMMM/short.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

image.o: QMMM/image.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

lattice.o: QMMM/lattice.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<

subgrids.o: subgrids.F90 params.mod
	$(COMPILER) $(COMPILERFLAGS1) -c $<
endif