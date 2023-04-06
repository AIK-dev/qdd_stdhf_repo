### Makefile

# This file contains the basic switches which have to be set
# by the user.

# Chose the compiler:
# Presently available choices:
#   * gfortran    (the GNU Fortran compiler)
#   * ifort       (the Intel Fortran compiler)
COMPILER = gfortran

# OS: select operating system. Ignored if: FFT_TYPE = FFTW,
# otherwise ignored.
# Available options:
#   * LINUX
#   * MAC
OS = LINUX

# FFT_TYPE: FFT solvers
# Available options:
#   * FFTW
#   * MKL
FFT_TYPE = FFTW

# OMP: Use OpenMP for threading support for FFTs and wave functions
# Available options:
#   * YES
#   * NO
OMP = YES

# DYNOMP: Parallelise s.p. wfs. Is ignored if: OMP = NO
# Available options:
#   * YES
#   * NO
DYNOMP = YES

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

# Set the location to the root of Intel's MKL distribution
#MKL_PATH = $(MKLROOT)

### SET MAIN LIBRARY LOCATIONS
# Set the location to the root of the FFTW subroutine package
#FFTW_PATH = $(HOME)/Developer/fftw-3.3.8


#####################################################################
#                  Main Makefile-body inclusion                     #
#####################################################################
### See mfBody.*.mk header for info
ifeq "$(strip $(COMPILER))" "gfortran"
  include Makefiles/mfBody.gnu.mk
endif

ifeq "$(strip $(COMPILER))" "ifort"
  include Makefiles/mfBody.intel.mk
endif
######################################################################
