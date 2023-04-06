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
       HEeigensystem.o mini.o zdiag.o orthmat.o stochastic.o

ifneq ($(strip $(NO_QDD_INFO)), YES)
  OBJS += QDD_info.o
endif

DEPS_EXTENDED =
ifeq ($(strip $(EXTENDED)),YES)
  OBJS += util_extended.o attachement.o V2ph.o zeroforce.o 
  DEPS_EXTENDED = util_extended.o attachement.o V2ph.o 
endif

ifeq ($(strip $(FSIC)),YES)
  OBJS += 2stUT.o 2st_util.o 2stUTc.o
endif

ifeq ($(strip $(RAREGAS)), YES)
  OBJS +=  functions.o pot_substrate.o forces_substrate.o md_substrate.o\
    short.o image.o lattice.o util_substrate.o
endif

# Additions for twostsic:
ifeq ($(strip $(FSIC)), YES)
  FSICDEPEND = 2st_util.o 2stUTc.o 2stUT.o
endif
