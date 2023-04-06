### Makefile.intel.avx2.mk

# Machine/architecture specific optimisations. Leave empty if unsure.
ARCH = -xCORE-AVX2

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
OS = MAC

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

### QDD version and source info. Flag to switch off storage of QDD version info,
# source code status and linked library versions.
# YES/NO: default should be 'NO', unless you are not working from a GIT repository,
# which is NOT recommended.
NO_QDD_INFO = YES

#####################################################################
#                         QDD option selection                      #
#####################################################################

# grid: FFT or finite difference (FFT options must be NO for finite differences)
FINDIFF = NO
NUMEROV = NO
# Switch to extended model with polarizable raregas
RAREGAS = NO
# Include private sector
EXTENDED = NO
FSIC = NO










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
