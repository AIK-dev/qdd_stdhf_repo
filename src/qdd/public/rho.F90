!-----calcrho------------------------------------------------------

#ifdef REALSWITCH
SUBROUTINE calcrhor(rho, q0)
#else
  SUBROUTINE calcrho(rho, q0)
#endif

! local density 'rho' for COMPLEX or REAL wavefunctions 'q0'
!
! Input:
! q0 = set of s.p. wavefunctions
! occup = occupation numbers, via module 'params'
! Output:
! rho = local density for spin-down and spin-up separately

    USE params
    USE util, ONLY: emoms
    IMPLICIT NONE

#ifdef REALSWITCH
    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
#else
    COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate) ! cPW
#endif
    REAL(DP), INTENT(OUT) :: rho(2*kdfull2)


    INTEGER :: ind, ishift, nb, i1, i2, i3
    REAL(DP) :: rhodif, rhotot


!-----------------------------------------------------------------



!k initialize densities:
    rho = 0D0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ind,i1,i2,i3,ishift,nb) SCHEDULE(STATIC)
    DO i3 = 1, kzbox
      DO nb = 1, nstate
        ishift = (ispin(nrel2abs(nb)) - 1)*nxyz ! store spin=2 in upper BLOCK
        ind = (i3 - 1)*(kxbox*kybox)
        DO i2 = 1, kybox; DO i1 = 1, kxbox
            ind = ind + 1
#ifdef REALSWITCH
            rho(ind + ishift) = rho(ind + ishift) + occup(nb)*q0(ind, nb)*q0(ind, nb)
#else
            rho(ind + ishift) = rho(ind + ishift) + occup(nb)*(CONJG(q0(ind, nb)))*q0(ind, nb)
#endif
          END DO; END DO
      END DO
    END DO
!$OMP END PARALLEL DO


! reorder to total density in lower BLOCK (1:nxyz)
! and difference density in upper BLOCK (nxyz+1:2nxyz)
    IF (numspin == 2) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ind,rhotot) SCHEDULE(STATIC)
      DO ind = 1, nxyz
        rhotot = rho(ind) + rho(ind + nxyz)
        rho(ind + nxyz) = (rho(ind) - rho(ind + nxyz))/MAX(rhotot, 1D-8)
        rho(ind) = rhotot
      END DO
!$OMP END PARALLEL DO
    ELSE
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(rhotot) SCHEDULE(STATIC)
      DO ind = 1, nxyz
        rho(ind) = rho(ind) + rho(ind + nxyz)
        rho(ind + nxyz) = 0D0
      END DO
!$OMP END PARALLEL DO
    END IF
    CALL emoms(rho) ! moments for the whole system (IN qe)


    RETURN
#ifdef REALSWITCH
  END SUBROUTINE calcrhor
#else
END SUBROUTINE calcrho
#endif

#ifdef COMPLEXSWITCH
!-----calc_current------------------------------------------------------
SUBROUTINE calc_current(current, q0)

! current 'current' for set of COMPLEX wavefunctions 'q0'

  USE params
  USE kinetic
  USE util, ONLY: prifld

  IMPLICIT NONE
#if(dynomp)
  REAL(DP), ALLOCATABLE :: homecurr(:, :, :)
#endif
  COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
  REAL(DP), INTENT(OUT) :: current(kdfull2, 3)

  INTEGER :: nb,ithr,ntm,nta
#if(dynomp)
    COMPLEX(DP), ALLOCATABLE :: dq0(:,:)
!    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS,OMP_GET_NUM_THREADS
#else
  COMPLEX(DP), ALLOCATABLE :: dq0(:)
#endif
!-----------------------------------------------------------------


#if(dynomp)
  ALLOCATE (dq0(kdfull2, 0:nthr))
  ALLOCATE (homecurr(kdfull2, 3, 0:nthr))
  homecurr = 0D0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ithr) SCHEDULE(STATIC)
  DO nb = 1, nstate
    ithr = OMP_GET_THREAD_NUM()
    IF(nb==1) THEN
      ntm = OMP_GET_MAX_THREADS()
      nta = OMP_GET_NUM_THREADS()
      IF(ntm.NE.nta) WRITE(*,*) &
      'CURRENT: threads not exhausted: nb,ithr,maxthreads,actthreads=',&
       nb,ithr,ntm,nta
    END IF
    CALL xgradient_rspace(q0(1, nb), dq0(:,ithr))
    homecurr(:, 1, ithr) = homecurr(:, 1, ithr) + occup(nb)*AIMAG(CONJG(q0(:, nb))*dq0(:,ithr))
!    WRITE(*,*) 'for NB=',nb,ithr,nthr,occup(nb)
!    CALL prifld(current(:, 2), 'j_1')
    CALL ygradient_rspace(q0(1, nb), dq0(:,ithr))
    homecurr(:, 2, ithr) = homecurr(:, 2, ithr) + occup(nb)*AIMAG(CONJG(q0(:, nb))*dq0(:,ithr))
    CALL zgradient_rspace(q0(1, nb), dq0(:,ithr))
    homecurr(:, 3, ithr) = homecurr(:, 3, ithr) + occup(nb)*AIMAG(CONJG(q0(:, nb))*dq0(:,ithr))
  END DO
!$OMP END PARALLEL DO
  current = 0D0
  DO ithr=0,nthr
    current(:,:) = current(:,:) + homecurr(:,:,ithr)
  END DO
  DEALLOCATE(homecurr)
#else
  ALLOCATE (dq0(kdfull2))
! reset
  current = 0D0
! accumulate
  DO nb = 1, nstate
    CALL xgradient_rspace(q0(1, nb), dq0(:))
    current(:, 1) = current(:, 1) + occup(nb)*AIMAG(CONJG(q0(:, nb))*dq0(:))
!    WRITE(*,*) 'for NB=',nb
!    CALL prifld(current(:, 2), 'j_1')
    CALL ygradient_rspace(q0(1, nb), dq0(:))
    current(:, 2) = current(:, 2) + occup(nb)*AIMAG(CONJG(q0(:, nb))*dq0(:))
    CALL zgradient_rspace(q0(1, nb), dq0(:))
    current(:, 3) = current(:, 3) + occup(nb)*AIMAG(CONJG(q0(:, nb))*dq0(:))
  END DO
#endif

  DEALLOCATE (dq0)

  RETURN

END SUBROUTINE calc_current
#endif

!-----spmoms------------------------------------------------------------

#ifdef REALSWITCH
SUBROUTINE spmomsr(wfr, iunit)
#else
  SUBROUTINE spmoms(wf, iunit)
#endif

! Multipole moments of single-particle densities from wf's.
! Input:
! wfr = set of REAL single particle wavefunctions
! wf = set of COMPLEX s.p. wavefunctions
! iunit = unit number for output
! Output:
! qeorb = array of s.p. moments, via module 'params'

    USE params
    IMPLICIT NONE

#ifdef REALSWITCH
    REAL(DP), INTENT(IN) :: wfr(kdfull2, kstate)
#else
    COMPLEX(DP), INTENT(IN) :: wf(kdfull2, kstate)
#endif
    INTEGER, INTENT(IN) :: iunit

    LOGICAL, PARAMETER :: ttest = .false.
    INTEGER :: ind, ix, iy, iz, j, k, nbe
    REAL(DP) :: r2el, s, xcmel, ycmel, zcmel, x1, y1, z1, x2, y2, z2
    REAL(DP), ALLOCATABLE :: qeorb(:, :)


!----------------------------------------------------------------------

    ALLOCATE (qeorb(kstate, 11))

      WRITE (iunit, '(a)') 'protocol of s.p. moments:', &
        ' state energy x y z variance xx yy zz xy xz yz'

    xcmel = 0D0
    ycmel = 0D0
    zcmel = 0D0
    r2el = 0D0

    DO nbe = 1, nstate
      DO k = 1, 11
        qeorb(nbe, k) = 0D0
      END DO

      ind = 0
      DO iz = minz, maxz
        z1 = (iz - nzsh)*dz
        z2 = z1*z1
        DO iy = miny, maxy
          y1 = (iy - nysh)*dy
          y2 = y1*y1
          DO ix = minx, maxx
            ind = ind + 1
            IF ((ix /= nx2) .AND. (iy /= ny2) .AND. (iz /= nz2)) THEN
              x1 = (ix - nxsh)*dx
              x2 = x1*x1
#ifdef REALSWITCH
              s = wfr(ind, nbe)*wfr(ind, nbe)
#else
              s = wf(ind, nbe)*CONJG(wf(ind, nbe))
#endif
              qeorb(nbe, 1) = amoy(nbe)
! monopole
              qeorb(nbe, 2) = qeorb(nbe, 2) + s
! dipole
              qeorb(nbe, 3) = qeorb(nbe, 3) + s*x1
              qeorb(nbe, 4) = qeorb(nbe, 4) + s*y1
              qeorb(nbe, 5) = qeorb(nbe, 5) + s*z1
! quadrupole
              qeorb(nbe, 6) = qeorb(nbe, 6) + s*x2
              qeorb(nbe, 7) = qeorb(nbe, 7) + s*y2
              qeorb(nbe, 8) = qeorb(nbe, 8) + s*z2

              qeorb(nbe, 9) = qeorb(nbe, 9) + s*x1*y1
              qeorb(nbe, 10) = qeorb(nbe, 10) + s*z1*x1
              qeorb(nbe, 11) = qeorb(nbe, 11) + s*z1*y1
            END IF
          END DO
        END DO
      END DO

      DO k = 2, 11
        qeorb(nbe, k) = qeorb(nbe, k)*dvol
      END DO

      WRITE (iunit, '(i4,f7.3,4f6.2,2x,6f7.1)') nbe, qeorb(nbe, 1), &
        (qeorb(nbe, j), j=3, 5), &
        SQRT(qeorb(nbe, 6) + qeorb(nbe, 7) + qeorb(nbe, 8)), (qeorb(nbe, j), j=6, 11)

      qeorb_all(nbe, :) = qeorb(nbe, :)
    END DO


    DEALLOCATE (qeorb)

      DO nbe = 1, nstate_all
        xcmel = xcmel + qeorb_all(nbe, 3)
        ycmel = ycmel + qeorb_all(nbe, 4)
        zcmel = zcmel + qeorb_all(nbe, 5)
        r2el = r2el + qeorb_all(nbe, 6) + qeorb_all(nbe, 7) + qeorb_all(nbe, 8)
      END DO

      xcmel = xcmel/nstate_all
      ycmel = ycmel/nstate_all
      zcmel = zcmel/nstate_all
      r2el = SQRT(r2el/nstate_all)

      WRITE (iunit, '(a11,4f6.2)') 'average: ', xcmel, ycmel, zcmel, r2el

    RETURN
#ifdef REALSWITCH
  END SUBROUTINE spmomsr
#else
END SUBROUTINE spmoms
#endif

!-----------------------------------------------------------------

!**************************from 1D TO 3D************************
#ifdef REALSWITCH
SUBROUTINE from1Dto3D(va, a, xva, yva, zva)
#else
  SUBROUTINE from1Dto3Dc(va, a, xva, yva, zva)
#endif

! Maps array 'va' given in linear storage to array 'a' in 3D storage.
! Input:
! va = array in linear storage
! xva,yva,zva = x-, y-, z-dimensions of input and output
! Output:
! a = same array 3D storage
!
! Comment: This routine may be replaced by point & TARGET construction.

    USE params
    IMPLICIT NONE

    INTEGER :: xva, yva, zva, ia, ja, ka, i0
#ifdef REALSWITCH
    REAL(DP), INTENT(IN) :: va(xva*yva*zva)
    REAL(DP), INTENT(OUT) :: a(xva, yva, zva)
#else
    COMPLEX(DP), INTENT(IN) :: va(xva*yva*zva)
    COMPLEX(DP), INTENT(OUT) :: a(xva, yva, zva)
#endif

    i0 = 0
    DO ka = 1, nz2
      DO ja = 1, ny2
        DO ia = 1, nx2
          i0 = i0 + 1
          a(ia, ja, ka) = va(i0)
        END DO
      END DO
    END DO

    RETURN
#ifdef REALSWITCH
  END SUBROUTINE from1Dto3D
#else
END SUBROUTINE from1Dto3Dc
#endif
!**************************from 3D TO 1D************************
#ifdef REALSWITCH
SUBROUTINE from3Dto1D(va, a, xva, yva, zva)
#else
  SUBROUTINE from3Dto1Dc(va, a, xva, yva, zva)
#endif

! Maps array 'a' given in 3D storage to array 'va' in linear storage.
! Input:
! a = array 3D storage
! xva,yva,zva = x-, y-, z-dimensions of input and output
! Output:
! va = same array in linear storage
!
! Comment: This routine may be replaced by point & target construction.

    USE params
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: xva
    INTEGER, INTENT(IN) :: yva
    INTEGER, INTENT(IN) :: zva
#ifdef REALSWITCH
    REAL(DP), INTENT(OUT) :: va(xva*yva*zva)
    REAL(DP), INTENT(IN) :: a(xva, yva, zva)
#else
    COMPLEX(DP), INTENT(OUT) :: va(xva*yva*zva)
    COMPLEX(DP), INTENT(IN) :: a(xva, yva, zva)
#endif

    INTEGER :: i0, ia, ja, ka

    i0 = 0
    DO ka = 1, zva
      DO ja = 1, yva
        DO ia = 1, xva
          i0 = i0 + 1
          va(i0) = a(ia, ja, ka)
        END DO
      END DO
    END DO

    RETURN
#ifdef REALSWITCH
  END SUBROUTINE from3Dto1D
#else
END SUBROUTINE from3Dto1Dc
#endif
