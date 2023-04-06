### mfTargets.mk
# Contains all the compiler-mutual implicit and explicit
# object dependencies and build recipes.

# Implicit recipes for objects:
%.o: %.F90
	$(COMPILER) $(COMPILERFLAGS1) -c $<

kinetic.o: params.o fftw.o
coulsolv.o: params.o fftw.o kinetic.o
static.o: params.o coulsolv.o util.o $(FSICDEPEND) $(DEPS_EXTENDED)
dynamic.o: params.o kinetic.o util.o rta.o orthmat.o stochastic.o $(FSICDEPEND) $(DEPS_EXTENDED)
util.o: params.o kinetic.o


# Explicit recipes for objects:
rho.o: rho.F90 params.o kinetic.o util.o $(DEPS_EXTENDED)
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o $@ -c $<

rhoc.o: rho.F90 params.o kinetic.o util.o $(DEPS_EXTENDED)
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

localizer.o: localize.F90 kinetic.o util.o
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o $@ -c $<

localize.o: localize.F90 kinetic.o util.o
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

sicnew.o: sicnew.F90 params.o kinetic.o coulsolv.o util.o
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o $@ -c $<

sicnewc.o: sicnew.F90 params.o kinetic.o coulsolv.o util.o
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

nonloc.o: nonloc.F90 kinetic.o
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o $@ -c $<

nonlocc.o: nonloc.F90 kinetic.o
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

fftw.o: fftw.F90
ifeq ($(strip $(FFT_TYPE)), MKL)
	$(COMPILER) $(COMPILERFLAGS1) -I$(MKL_INCLUDE) -c $<
else ifeq ($(strip $(FFT_TYPE)), FFTW)
	$(COMPILER) $(COMPILERFLAGS1) -I$(FFTW_INCLUDE) -c $<
endif

ifneq ($(strip $(NO_QDD_INFO)), YES)
main.o: main.F90 params.o kinetic.o util.o coulsolv.o orthmat.o rta.o QDD_info.o $(FSICDEPEND) $(DEPS_EXTENDED)
	$(COMPILER) $(COMPILERFLAGS3) -DQDD_INFO -c $<
else
main.o: main.F90 params.o kinetic.o util.o coulsolv.o orthmat.o rta.o $(FSICDEPEND) $(DEPS_EXTENDED)
	$(COMPILER) $(COMPILERFLAGS3) -c $<
endif

givens.o: givens.F90 params.o
	$(COMPILER) $(COMPILERFLAGS3) -c $<

restart.o: restart.F90 params.o kinetic.o $(FSICDEPEND)
	$(COMPILER) $(COMPILERFLAGS3) -DREALSWITCH -o $@ -c $<

restartc.o: restart.F90 params.o kinetic.o $(FSICDEPEND)
	$(COMPILER) $(COMPILERFLAGS3) -DCOMPLEXSWITCH -o $@ -c $<

parallele.o: parallele.F90 params.o $(DEPS_EXTENDED)
	$(COMPILER) $(COMPILERFLAGS3) -c $<

rta.o: rta.F90 params.o kinetic.o util.o
	$(COMPILER) $(COMPILERFLAGS3) -c $<

orthmat.o: orthmat.F90 params.o
	$(COMPILER) $(COMPILERFLAGS3) -c -cpp $<

zdiag.o: zdiag.f
	$(COMPILER) $(COMPILERFLAGS3) -c $<

mini.o: mini.f
	$(COMPILER) $(COMPILERFLAGS3) -c $<

HEeigensystem.o: HEeigensystem.f
	$(COMPILER) $(COMPILERFLAGS3) -c -cpp $<

# Additions for twostsic:
ifeq ($(strip $(FSIC)), YES)
2st_util.o: fsic/2st_util.F90 params.o kinetic.o orthmat.o util.o
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o $@ -c $<

2stUTc.o: fsic/2stUT.F90 params.o kinetic.o orthmat.o 2stUT.o
	$(COMPILER) $(COMPILERFLAGS1) -DCOMPLEXSWITCH -o $@ -c $<

2stUT.o: fsic/2stUT.F90 params.o kinetic.o 2st_util.o
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o $@ -c $<
endif

ifeq ($(strip $(EXTENDED)),YES) # NEEDS 'FSIC = YES'
attachement.o: extended/attachement.F90 params.o util.o
	$(COMPILER) $(COMPILERFLAGS3) -c $<

util_extended.o: extended/util_extended.F90
	$(COMPILER) $(COMPILERFLAGS3) -c $<

V2ph.o: extended/V2ph.f90 params.o kinetic.o
	$(COMPILER) $(COMPILERFLAGS1) -o $@ -c $<

zeroforce.o: extended/zeroforce.F90 params.o kinetic.o
	$(COMPILER) $(COMPILERFLAGS1) -o $@ -c $<

stochastic.o: stochastic.F90 params.o kinetic.o
	$(COMPILER) $(COMPILERFLAGS1) -o $@ -c $<
endif



# QM/MM routines
ifeq ($(strip $(RAREGAS)), YES)
forces_substrate.o: QMMM/forces_substrate.F90 params.o util.o
	$(COMPILER) $(COMPILERFLAGS1) -c $<

util_substrate.o: QMMM/util_substrate.F90 params.o util.o
	$(COMPILER) $(COMPILERFLAGS1) -c $<

functions.o: QMMM/functions.F90 params.o
	$(COMPILER) $(COMPILERFLAGS1) -c $<

pot_substrate.o: QMMM/pot_substrate.F90 params.o
	$(COMPILER) $(COMPILERFLAGS1) -c $<

md_substrate.o: QMMM/md_substrate.F90 params.o util.o
	$(COMPILER) $(COMPILERFLAGS1) -c $<

short.o: QMMM/short.F90 params.o
	$(COMPILER) $(COMPILERFLAGS1) -c $<

image.o: QMMM/image.F90 params.o util.o coulsolv.o
	$(COMPILER) $(COMPILERFLAGS1) -c $<

lattice.o: QMMM/lattice.F90 params.o
	$(COMPILER) $(COMPILERFLAGS1) -c $<

subgrids.o: subgrids.F90 params.o
	$(COMPILER) $(COMPILERFLAGS1) -c $<
endif

# CUDA routines
ifeq ($(strip $(CUDA)), YES)

endif