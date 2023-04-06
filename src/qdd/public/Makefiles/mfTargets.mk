### mfTargets.mk
# Contains all the compiler-mutual implicit and explicit
# object dependencies and build recipes.

# Implicit recipes for objects:
%.o: %.F90
	$(COMPILER) $(COMPILERFLAGS1) -c $<

kinetic.o: params.o fftw.o
coulsolv.o: params.o fftw.o kinetic.o
static.o: params.o coulsolv.o util.o
dynamic.o: params.o kinetic.o util.o rta.o orthmat.o
util.o: params.o kinetic.o


# Explicit recipes for objects:
rho.o: rho.F90 params.o kinetic.o util.o
	$(COMPILER) $(COMPILERFLAGS1) -DREALSWITCH -o $@ -c $<

rhoc.o: rho.F90 params.o kinetic.o util.o
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

main.o: main.F90 params.o kinetic.o util.o coulsolv.o orthmat.o rta.o
	$(COMPILER) $(COMPILERFLAGS3) -c $<

givens.o: givens.F90 params.o
	$(COMPILER) $(COMPILERFLAGS3) -c $<

restart.o: restart.F90 params.o kinetic.o
	$(COMPILER) $(COMPILERFLAGS3) -DREALSWITCH -o $@ -c $<

restartc.o: restart.F90 params.o kinetic.o
	$(COMPILER) $(COMPILERFLAGS3) -DCOMPLEXSWITCH -o $@ -c $<

parallele.o: parallele.F90 params.o
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
