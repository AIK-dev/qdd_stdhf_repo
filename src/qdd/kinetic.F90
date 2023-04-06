
! This file includes FFT or finite differences version of
! module kinetic

#if(netlib_fft|fftw_cpu)
#include "fft.F90"
#endif

#if(findiff|numerov)
#include "findiff/findiff.F90"
#endif

