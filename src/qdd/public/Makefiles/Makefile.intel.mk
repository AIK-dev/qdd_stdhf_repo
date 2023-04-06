### Makefile.intel.mk

# FFT_TYPE: FFT solvers
# Available options:
#   * FFTW
#   * MKL
FFT_TYPE = MKL

# OS: select operating system. Ignored if: FFT_TYPE = FFTW,
# otherwise ignored.
# Available options:
#   * LINUX
#   * MAC
OS = LINUX

# Set the location to the root of Intel's MKL distribution
MKL_PATH = $(MKLROOT)

### SET MAIN LIBRARY LOCATIONS
# Set the location to the root of the FFTW subroutine package
FFTW_PATH = $(HOME)/Developer/fftw-3.3.8

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

### OpenMP debugging, just for me... To be deleted after testing...
OMP_DEBUG = NO



#####################################################################
#                  Main Makefile-body inclusion                     #
#####################################################################
### See mfBody.intel.mk header for info
ifneq (,$(wildcard mfBody.intel.mk))
  include mfBody.intel.mk
else 
  include Makefiles/mfBody.intel.mk
endif
######################################################################
