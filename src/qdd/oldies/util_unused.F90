
! **************************

REAL(DP) FUNCTION wfnorm_r(psi1, psi2)

! Overlap of wavefunctions 'psi1' and 'psi2'. REAL version

  USE params
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN OUT) :: psi1(kdfull2)
  COMPLEX(DP), INTENT(IN OUT) :: psi2(kdfull2)

  INTEGER :: ii
  REAL(DP) :: acc

  acc = 0D0
  DO ii = 1, nxyz
    acc = REAL(psi1(ii), DP)*REAL(psi2(ii), DP) + AIMAG(psi1(ii))*AIMAG(psi2(ii)) + acc
  END DO
  wfnorm_r = acc*dvol
  RETURN
END FUNCTION wfnorm_r
! **************************

REAL(DP) FUNCTION wfnorm_i(psi1, psi2)

! Overlap of wavefunctions 'psi1' and 'psi2'. COMPLEX version

  USE params
  IMPLICIT NONE
  COMPLEX(DP), INTENT(IN OUT) :: psi1(kdfull2)
  COMPLEX(DP), INTENT(IN OUT) :: psi2(kdfull2)

  INTEGER :: ii
  REAL(DP) :: acc
  acc = 0D0
  DO ii = 1, nxyz
    acc = -REAL(psi1(ii), DP)*AIMAG(psi2(ii)) + AIMAG(psi1(ii))*REAL(psi2(ii), DP) + acc
  END DO
  wfnorm_i = acc*dvol
  RETURN
END FUNCTION wfnorm_i

!********************** from restart.F90 *******************

! **************************

SUBROUTINE restherm()

! **************************

  USE params
  USE kinetic
  IMPLICIT NONE
  INTEGER :: iact, ion

  OPEN (UNIT=20, STATUS='unknown', FORM='UNFORMATTED', FILE='therm')
! nxyz=nx2*ny2*nz2
  READ (20) iact
  IF (iact /= irest) THEN
    WRITE (7, *) 'iact=', iact
    STOP 'bad irest IN for005dyn !'
  END IF
  DO ion = 1, nion
    READ (20) cx(ion), cy(ion), cz(ion)
    READ (20) cpx(ion), cpy(ion), cpz(ion)
  END DO
  CLOSE (UNIT=20, STATUS='keep')
  RETURN
END SUBROUTINE restherm

! **************************

SUBROUTINE addcluster(psi, outna)

! **************************

  USE params
  USE kinetic
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
  CHARACTER(LEN=13), INTENT(IN OUT) :: outna

  INTEGER :: i, iact, idum, ii, ion, k, nb, nelectt, niont, nspdwt, nstatet
  OPEN (UNIT=ifile, STATUS='unknown', FORM='UNFORMATTED', FILE='SAVE.'//outna)

! READ the iteration WHERE the DATA has been saved last:

  READ (ifile) iact, nstatet, nelectt, niont, nspdwt
  DO i = 1, nstatet
    READ (ifile) occup(i + nstate)
  END DO

  irest = iact

! READ wavefunctions:
  IF (nelectt > 0) THEN
    DO nb = 1, nstatet
    DO i = 1, nxyz
      READ (ifile) psi(i, nb + nstate)
    END DO
    END DO
  END IF

  DO i = 1, ksttot
    IF (i <= nstatet) THEN
      READ (ifile) ispin(i + nstate), nrel2abs(i + nstate), nabs2rel(i + nstate)
    ELSE
      READ (ifile) idum, idum, idum
    END IF
  END DO

! READ protonic coordinates and momenta
  WRITE (6, *) nion, niont
  DO ion = 1, niont
    READ (ifile) cx(ion + nion), cy(ion + nion), cz(ion + nion), np(ion + nion)
    READ (ifile) cpx(ion + nion), cpy(ion + nion), cpz(ion + nion), np(ion + nion)
  END DO

#if(raregas)
  IF (isurf /= 0) THEN

    DO i = 1, nc
      READ (ifile) imobc(i)
      READ (ifile) xc(i), yc(i), zc(i)
      READ (ifile) pxc(i), pyc(i), pzc(i)
    END DO
    DO i = 1, NE
      READ (ifile) imobe(i)
      READ (ifile) xe(i), ye(i), ze(i)
      READ (ifile) pxe(i), pye(i), pze(i)
    END DO
    DO i = 1, nk
      READ (ifile) imobk(i)
      READ (ifile) xk(i), yk(i), zk(i)
      READ (ifile) pxk(i), pyk(i), pzk(i)
    END DO

    DO i = 1, kdfull2
      READ (ifile) potfixedion(i)
    END DO

  END IF
#endif

  IF (nelectt > 0) THEN
    DO k = 1, kmom
      READ (ifile) qe(k)
    END DO
    DO i = 1, 3
      READ (ifile) se(i)
    END DO
  END IF

  IF (nabsorb > 0) THEN
    DO ii = 1, kdfull2
      READ (ifile) rhoabso(ii)
    END DO
  END IF

  CLOSE (UNIT=ifile, STATUS='keep')

  nstate = nstate + nstatet
  nelect = nelect + nelectt
  nion = nion + niont
  nspdw = nspdw + nspdwt
  nion2 = nion

  RETURN
END SUBROUTINE addcluster

! From 'orthmat' package

!-------------------------------------------------------
! cmplxsic version
!-------------------------------------------------------

REAL(DP) FUNCTION matdorth_cmplxsic(aa, n, ndim)

  COMPLEX(DP), INTENT(IN) :: aa(n, n)
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i, j
  matdorth_cmplxsic = 0D0
  DO i = 1, ndim
    DO j = i, ndim
      IF (i == j) THEN
        matdorth_cmplxsic = matdorth_cmplxsic + SUM(aa(i, 1:ndim)*aa(j, 1:ndim)) - 1D0
      ELSE
        matdorth_cmplxsic = matdorth_cmplxsic + SUM(aa(i, 1:ndim)*aa(j, 1:ndim))
      END IF
    END DO
  END DO
END FUNCTION matdorth_cmplxsic

! From 'ionmd.F90'

! *******************

SUBROUTINE stream(rho)

! *******************

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)

  INTEGER :: i
  REAL(DP) :: query2, query3, query4
  REAL(DP), ALLOCATABLE :: rhos(:), drho(:)

  ALLOCATE (rhos(kdfull2), drho(kdfull2))

  query2 = qe(2)*qold2
  qold2 = qe(2)
  query3 = qe(3)*qold3
  qold3 = qe(3)
  query4 = qe(4)*qold4
  qold4 = qe(4)
!test
  query2 = 1.0D0
  query3 = 1.0D0
!test
  IF (tfs == 0D0) THEN
    iquery4 = 0
    WRITE (6, *) 'iquery4=', iquery4
    DO i = 1, nxyz
      rhos(i) = 0D0
    END DO
  END IF
  IF (tfs > 0D0) THEN
    IF (query2 < 0D0) THEN
      iquery4 = iquery4 + 1
      WRITE (6, *) 'iquery4=', iquery4
      IF (MOD(iquery4, 10) == 0) THEN
        DO i = 1, nxyz
          drho(i) = rho(i) - rhos(i)
        END DO
        CALL write_density(rhos, 'sdensi')
        CALL write_density(drho, 'sderiv')
! STOP 'calculated a derivative of density after x-stimulation'
      END IF
    END IF
    IF (query3 < 0D0) THEN
      iquery4 = iquery4 + 1
      WRITE (6, *) 'iquery4=', iquery4
      IF (MOD(iquery4, 10) == 0) THEN
        DO i = 1, nxyz
          drho(i) = rho(i) - rhos(i)
        END DO
        CALL write_density(rhos, 'sdensi')
        CALL write_density(drho, 'sderiv')
! STOP 'calculated a derivative of density after y-stimulation'
      END IF
    END IF
    IF (query4 < 0D0) THEN
      iquery4 = iquery4 + 1
      WRITE (6, *) 'iquery4=', iquery4
      IF (MOD(iquery4, 10) == 0) THEN
        DO i = 1, nxyz
          drho(i) = rho(i) - rhos(i)
        END DO
        CALL write_density(rhos, 'sdensi')
        CALL write_density(drho, 'sderiv')
! STOP 'calculated a derivative of density after z-stimulation'
      END IF
    END IF
  END IF
  rhos = rho

  DEALLOCATE (rhos, drho)

  RETURN
END SUBROUTINE stream

! ---eltherm----------------------------------------------eltherm-------

SUBROUTINE eltherm(rho, psi, expjx, expjy, expjz, eeth)

! Computes electronic thermal energy as variance of
! electronic s.p. velicities
! e_th = sum_i int[d^3r 0.5D0*rho_i*(v_av - v_i)**2]
!
! Input:
! rho = local electron density
! psi = set of s.p. wavefunctions
! Output:
! expjx,... = variances along x,y,z
! eeth = total thermal energy

  USE params
  USE kinetic
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
  REAL(DP), INTENT(OUT) :: expjx
  REAL(DP), INTENT(OUT) :: expjy
  REAL(DP), INTENT(OUT) :: expjz
  REAL(DP), INTENT(OUT) :: eeth

  INTEGER :: i, ind, n, nb
  REAL(DP) :: ajalpha, rhoalpha
  REAL(DP) :: dkx, dky, dkz, q1, tel
  REAL(DP) :: exjx(ksttot), exjy(ksttot), exjz(ksttot)
#if(mpi)
  REAL(DP) :: ex, ey, ez
#endif
  COMPLEX(DP) :: test

  REAL(DP), DIMENSION(:), ALLOCATABLE :: ajtx, ajty, ajtz
  COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: p, q2

#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
  REAL(DP), ALLOCATABLE :: aj(:)
  ALLOCATE (aj(kdfull2))
#endif

  IF (.NOT. ALLOCATED(akx)) STOP "ELTHERM requires FFT"
  ALLOCATE (p(kdfull2), q2(kdfull2))
  ALLOCATE (ajtx(kdfull2), ajty(kdfull2), ajtz(kdfull2))

! init derivative

  dkx = pi/(dx*REAL(nx, DP))
  dky = pi/(dy*REAL(ny, DP))
  dkz = pi/(dz*REAL(nz, DP))

! compute the currents for each direction: ajtx/y/z
! and their expectation values: exjx/y/z

  q1 = 0D0
  tel = 0D0
  DO i = 1, kdfull2
    ajtx(i) = 0D0
    ajty(i) = 0D0
    ajtz(i) = 0D0
  END DO

  DO nb = 1, nstate

    exjx(nb) = 0D0
    exjy(nb) = 0D0
    exjz(nb) = 0D0

#if(netlib_fft|fftw_cpu)
    CALL fftf(psi(1, nb), q2)

    DO ind = 1, kdfull2
      q2(ind) = q2(ind)*akx(ind)
    END DO

    CALL fftback(q2, p)
#endif

    DO ind = 1, kdfull2
      test = eye/2D0*(CONJG(psi(ind, nb))*p(ind) - psi(ind, nb)*CONJG(p(ind)))
      ajalpha = test
      exjx(nb) = exjx(nb) + ajalpha*CONJG(psi(ind, nb))*psi(ind, nb)
      ajtx(ind) = ajtx(ind) + ajalpha
    END DO

#if(netlib_fft|fftw_cpu)
    CALL fftf(psi(1, nb), q2)

    DO ind = 1, kdfull2
      q2(ind) = q2(ind)*aky(ind)
    END DO

    CALL fftback(q2, p)
#endif

    DO ind = 1, kdfull2
      test = eye/2D0*(CONJG(psi(ind, nb))*p(ind) - psi(ind, nb)*CONJG(p(ind)))
      ajalpha = test
      exjy(nb) = exjy(nb) + ajalpha*CONJG(psi(ind, nb))*psi(ind, nb)
      ajty(ind) = ajty(ind) + ajalpha
    END DO

#if(netlib_fft|fftw_cpu)
    CALL fftf(psi(1, nb), q2)

    DO ind = 1, kdfull2
      q2(ind) = q2(ind)*akz(ind)
    END DO

    CALL fftback(q2, p)
#endif

    DO ind = 1, kdfull2
      test = eye/2D0*(CONJG(psi(ind, nb))*p(ind) - psi(ind, nb)*CONJG(p(ind)))
      ajalpha = test
      exjz(nb) = exjz(nb) + ajalpha*CONJG(psi(ind, nb))*psi(ind, nb)
      ajtz(ind) = ajtz(ind) + ajalpha
    END DO

  END DO !loop over states

! compute sum of current expectation values IN either direction:

  expjx = 0D0
  expjy = 0D0
  expjz = 0D0
  DO n = 1, nstate
    expjx = expjx + exjx(n)
    expjy = expjy + exjy(n)
    expjz = expjz + exjz(n)
  END DO
  expjx = expjx*dvol
  expjy = expjy*dvol
  expjz = expjz*dvol
#if(mpi)
  DO ind = 1, kdfull2
    aj(ind) = ajtx(ind)
  END DO
  CALL pi_allreduce(aj, ajtx, kdfull2, mpi_double_precision, mpi_sum, &
                    mpi_comm_world, mpi_ierror)
  DO ind = 1, kdfull2
    aj(ind) = ajty(ind)
  END DO
  CALL pi_allreduce(aj, ajty, kdfull2, mpi_double_precision, mpi_sum, &
                    mpi_comm_world, mpi_ierror)
  DO ind = 1, kdfull2
    aj(ind) = ajtz(ind)
  END DO
  CALL pi_allreduce(aj, ajtz, kdfull2, mpi_double_precision, mpi_sum, &
                    mpi_comm_world, mpi_ierror)

  CALL pi_allreduce(expjx, ex, 1, mpi_double_precision, mpi_sum, mpi_comm_world, mpi_ierror)
  CALL pi_allreduce(expjy, ey, 1, mpi_double_precision, mpi_sum, mpi_comm_world, mpi_ierror)
  CALL pi_allreduce(expjz, ez, 1, mpi_double_precision, mpi_sum, mpi_comm_world, mpi_ierror)
  expjx = ex
  expjy = ey
  expjz = ez
#endif

! compute average velocities from currents:

  DO ind = 1, kdfull2
    ajtx(ind) = ajtx(ind)/rho(ind)
    ajty(ind) = ajty(ind)/rho(ind)
    ajtz(ind) = ajtz(ind)/rho(ind)
  END DO

! finally compute thermal energy:
! e_th = sum_i int[d^3r 0.5D0*rho_i*(v_av - v_i)**2]

  tel = 0D0
  DO nb = 1, nstate

#if(netlib_fft|fftw_cpu)
    CALL fftf(psi(1, nb), q2)

    DO ind = 1, kdfull2
      q2(ind) = q2(ind)*akx(ind)
    END DO

    CALL fftback(q2, p)
#endif

    DO ind = 1, kdfull2
      test = eye/2D0*(CONJG(psi(ind, nb))*p(ind) - psi(ind, nb)*CONJG(p(ind)))
      ajalpha = test
      rhoalpha = CONJG(psi(ind, nb))*psi(ind, nb)
      ajalpha = ajalpha/rhoalpha
      tel = tel + 0.5D0*rhoalpha*(ajalpha - ajtx(ind))**2
    END DO

#if(netlib_fft|fftw_cpu)
    CALL fftf(psi(1, nb), q2)

    DO ind = 1, kdfull2
      q2(ind) = q2(ind)*aky(ind)
    END DO

    CALL fftback(q2, p)
#endif

    DO ind = 1, kdfull2
      test = eye/2D0*(CONJG(psi(ind, nb))*p(ind) - psi(ind, nb)*CONJG(p(ind)))
      ajalpha = test
      rhoalpha = CONJG(psi(ind, nb))*psi(ind, nb)
      ajalpha = ajalpha/rhoalpha
      tel = tel + 0.5D0*rhoalpha*(ajalpha - ajty(ind))**2
    END DO

#if(netlib_fft|fftw_cpu)
    CALL fftf(psi(1, nb), q2)

    DO ind = 1, kdfull2
      q2(ind) = q2(ind)*akz(ind)
    END DO

    CALL fftback(q2, p)
#endif

    DO ind = 1, kdfull2
      test = eye/2D0*(CONJG(psi(ind, nb))*p(ind) - psi(ind, nb)*CONJG(p(ind)))
      ajalpha = test
      rhoalpha = CONJG(psi(ind, nb))*psi(ind, nb)
      ajalpha = ajalpha/rhoalpha
      tel = tel + 0.5D0*rhoalpha*(ajalpha - ajtz(ind))**2
    END DO

  END DO !loop over states

  eeth = tel*dvol*0.5D0*hbar*hbar/2D0/ame !atomic units (h/2m)**2
#if(mpi)
  tel = eeth
  CALL pi_allreduce(tel, eeth, 1, mpi_double_precision, mpi_sum, mpi_comm_world, mpi_ierror)
#endif

  WRITE (6, '(a,f12.5)') ' electronic thermal energy: ', eeth
! WRITE(17,'(a,f12.5)') 'electronic thermal energy: ',eeth

!DEALLOCATE(akkx)
!DEALLOCATE(akky)
!DEALLOCATE(akkz)
  DEALLOCATE (p, q2)
  DEALLOCATE (ajtx, ajty, ajtz)
#if(mpi)
  DEALLOCATE (aj)
#endif

  RETURN
END SUBROUTINE eltherm

! **************************

SUBROUTINE write_density(drho, filename)
! WRITE density or derivative of density
!
! **************************
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: drho(kdfull2)
  CHARACTER(LEN=6), INTENT(IN) :: filename
  CHARACTER(LEN=4) :: ext
  INTEGER :: i

  WRITE (ext, '(a,i3.3)') ".", iquery4/10 ! e.g: iquery4 = 1234 ---> ext = ".123"
  ! iquery4 = 568 ---> ext = ".056"
  OPEN (UNIT=60, STATUS='unknown', FORM='UNFORMATTED', FILE=filename//'1'//ext)

  DO i = 1, nxyz
    WRITE (60) drho(i)
  END DO

  CLOSE (UNIT=60, STATUS='keep')

  RETURN
END SUBROUTINE write_density

! from 'init.F90'

!-----init_homfield-------------------------------------------------

SUBROUTINE init_homfield()

! Initialize a homogeneous electrical field and
! adds it TO the background field 'potFixedIon'.

  USE params
  IMPLICIT NONE

  INTEGER :: ind, ix, iy, iz
  REAL(DP) :: sc, vhom, x1, y1, z1
!------------------------------------------------------------------

  sc = ehomx**2 + ehomy**2 + ehomz**2
  sc = SQRT(sc)
  ehomx = ehomx/sc
  ehomy = ehomy/sc
  ehomz = ehomz/sc
  ind = 0
  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        ind = ind + 1
        vhom = e2*(x1*ehomx + y1*ehomy + z1*ehomz)*ehom0
        potfixedion(ind) = potfixedion(ind) + vhom
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE init_homfield
