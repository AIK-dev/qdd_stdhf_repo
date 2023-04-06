### mfTargets.mk
# Contains all the compiler-mutual implicit and explicit
# object dependencies and build recipes.

%.o: %.F90
	$(COMPILER) $(COMPILERFLAGS1) -c $<

# Explicit recipes for objects:

cuparams.o :  cuda/cuparams.F90 cuda/cupseudo.F90
	$(COMPILER)  $(COMPILERFLAGS1) -Mpreprocess -o $@ -c $<

cukinetic.o :  cuda/cukinetic.F90 cuparams.o
	$(COMPILER)  $(COMPILERFLAGS1) -Mpreprocess -o $@ -c $<

cutools.o :  cuda/cutools.F90 cuparams.o
	$(COMPILER)  $(COMPILERFLAGS1) -o $@ -c $<

cudadyn.o : cuda/cudadyn.F90 cutools.o cukinetic.o cucoulsolv.o
	$(COMPILER)  $(COMPILERFLAGS1) -Mpreprocess -o $@ -c $<

cucoulsolv.o : cuda/cucoulsolv.F90 cukinetic.o cutools.o
	$(COMPILER)  $(COMPILERFLAGS1) -Mpreprocess -o $@ -c $<

dynamic.o : cuda/cudynamic.F90 cukinetic.o cutools.o cucoulsolv.o cudadyn.o
	$(COMPILER)  $(COMPILERFLAGS1) -Mpreprocess -o $@ -c $<