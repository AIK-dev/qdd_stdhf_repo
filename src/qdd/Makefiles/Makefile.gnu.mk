### Makefile.gnu.mk

# Machine/architecture specific optimisations. Leave empty if unsure.
#ARCH = -xCORE-AVX2

# FFT_TYPE: FFT solvers
# Available options:
#   * FFTW
#   * MKL
FFT_TYPE = FFTW

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
FFTW_PATH = $(MKLROOT)/include/fftw

# OMP: Use OpenMP for threading support for FFTs and wave functions
# Available options:
#   * YES
#   * NO
OMP = NO

# DYNOMP: Parallelise s.p. wfs. Is ignored if: OMP = NO
# Available options:
#   * YES
#   * NO
DYNOMP = NO

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
# --> make sure that EXTENDED = YES is set
RAREGAS = NO
# Switch to extended model with full SIC
FSIC = NO
# Switch to private sector with general extensions
EXTENDED = NO










#####################################################################
#                  Main Makefile-body inclusion                     #
#####################################################################
### See mfBody.gnu.mk header for info
ifneq (,$(wildcard mfBody.gnu.mk))
  include mfBody.gnu.mk
else 
  include Makefiles/mfBody.gnu.mk
endif
######################################################################
