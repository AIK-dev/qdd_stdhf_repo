### mfObjects.mk
# Gathers all theobjects selected by the user in the 
# top-level Makefiles.
OBJS = params.o main.o kinetic.o restart.o restartc.o init.o\
       static.o dynamic.o lda.o util.o abso_bc.o\
       pseudosoft.o pseudogoed.o ionmd.o forces.o\
       carlo.o localize.o localizer.o\
       sicnew.o sicnewc.o rho.o rhoc.o nonloc.o nonlocc.o\
       schmid.o loc_mfield.o givens.o subgrids.o\
       parallele.o rta.o coulsolv.o\
       HEeigensystem.o mini.o zdiag.o orthmat.o

