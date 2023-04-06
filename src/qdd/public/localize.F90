
! ---localize-----------------------------------------------------------

#ifdef REALSWITCH
SUBROUTINE localizer(rho, psi, it)
#else
  SUBROUTINE localize(rho, psi, it)
#endif

! Computes and prints localization criterion of Becke et al.
! Dynamical version possibly not up TO date -- check!
!
! Input:
! rho = local density
! psi = set of s.p.wavefunctions
! it = nr. of iteration or time step IN calling routuine
! further variables communciated through MODULE 'params'

    USE params
    USE util, ONLY: safeopen
    USE kinetic
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: it
    REAL(DP), INTENT(IN) :: rho(2*kdfull2)
#ifdef REALSWITCH
    REAL(DP), INTENT(IN) :: psi(kdfull2, kstate)
    REAL(DP) :: p(kdfull2)
#else
    COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
    COMPLEX(DP) :: p(kdfull2)
#endif

    INTEGER :: i, idirection, ind, indadd, is, ix, iy, iz, midx, midy, midz, nb
    REAL(DP) :: dkx, dky, dkz, ocfac, rp, sign, sumpart, tf_fac, x1, y1, z1
    COMPLEX(DP) :: q2(kdfull2)
    INTEGER :: i1, i2, i3
    REAL(DP) :: zkx, zky, zkz
    REAL(DP), DIMENSION(:), ALLOCATABLE :: akk
    REAL(DP) :: average_localization(2)
    REAL(DP), DIMENSION(:), ALLOCATABLE :: ajt
    REAL(DP), DIMENSION(:), ALLOCATABLE :: drho, arho
    REAL(DP), DIMENSION(:), ALLOCATABLE :: tau

    LOGICAL, PARAMETER :: tupdate = .true.

    IF (.NOT. tupdate) STOP ' LOCALIZE not up to date '
    ALLOCATE (akk(kdfull2))
    ALLOCATE (ajt(kdfull2))
    ALLOCATE (drho(kdfull2))
    ALLOCATE (arho(kdfull2))
    ALLOCATE (tau(2*kdfull2))

! init derivative

    dkx = pi/(dx*REAL(nx, DP))
    dky = pi/(dy*REAL(ny, DP))
    dkz = pi/(dz*REAL(nz, DP))

! scaling factor for Thomas-Fermi energy
    tf_fac = 0.6D0*(6.0D0*pi**2)**(2D0/3D0)

! check availability of FFT
    IF (.NOT. ALLOCATED(akv)) STOP "LOCALIZE requires FFT"

! compute the currents, gradient of density and kinetic density

    average_localization = 0D0
    DO is = 2, 1, -1
      sign = (3 - 2*is)*1D0

      DO i = 1, kdfull2
        tau(i) = 0D0
      END DO

      DO idirection = 1, 3

! prepare momentum-space factor for actual direction
        IF (idirection .EQ. 1) THEN
          ind = 0
          DO i3 = 1, nz2
            DO i2 = 1, ny2
              DO i1 = 1, nx2
                IF (i1 >= (nx + 1)) THEN
                  zkx = (i1 - nx2 - 1)*dkx
                ELSE
                  zkx = (i1 - 1)*dkx
                END IF
                ind = ind + 1
                akk(ind) = -zkx
              END DO
            END DO
          END DO
        ELSEIF (idirection .EQ. 2) THEN
          ind = 0
          DO i3 = 1, nz2
            DO i2 = 1, ny2
              IF (i2 >= (ny + 1)) THEN
                zky = (i2 - ny2 - 1)*dky
              ELSE
                zky = (i2 - 1)*dky
              END IF
              DO i1 = 1, nx2
                ind = ind + 1
                akk(ind) = -zky
              END DO
            END DO
          END DO
        ELSEIF (idirection .EQ. 3) THEN
          ind = 0
          DO i3 = 1, nz2
            IF (i3 >= (nz + 1)) THEN
              zkz = (i3 - nz2 - 1)*dkz
            ELSE
              zkz = (i3 - 1)*dkz
            END IF
            DO i2 = 1, ny2
              DO i1 = 1, nx2
                ind = ind + 1
                akk(ind) = -zkz
              END DO
            END DO
          END DO
        END IF
        DO i = 1, kdfull2
          ajt(i) = 0D0
          drho(i) = 0D0
          arho(i) = 0D0
        END DO

        DO nb = 1, nstate
          IF (ispin(nb) .NE. is) CYCLE
          ocfac = occup(nb)*2D0

#ifdef REALSWITCH
          CALL rftf(psi(1, nb), q2)
#else
          CALL fftf(psi(1, nb), q2)
#endif

          DO ind = 1, kdfull2
            q2(ind) = q2(ind)*akk(ind)*eye
          END DO
#ifdef REALSWITCH
          CALL rfftback(q2, p)
#else

          CALL fftback(q2, p)
#endif
          DO ind = 1, kdfull2
#ifdef REALSWITCH
            drho(ind) = drho(ind) + ocfac*(psi(ind, nb)*p(ind))
            tau(ind) = tau(ind) + occup(nb)*(p(ind)*p(ind))
            arho(ind) = arho(ind) + occup(nb)*(psi(ind, nb)*psi(ind, nb))
#else
            drho(ind) = drho(ind) + ocfac*REAL(CONJG(psi(ind, nb))*p(ind), DP)
            ajt(ind) = ajt(ind) - occup(nb)*AIMAG(CONJG(psi(ind, nb))*p(ind))
            tau(ind) = tau(ind) + occup(nb)*REAL(CONJG(p(ind))*p(ind), DP)
            arho(ind) = arho(ind) + occup(nb)*REAL(CONJG(psi(ind, nb))*psi(ind, nb), DP)
#endif
          END DO

        END DO !loop over states

! accumulate contribution from actual direction on 'tau'
        DO ind = 1, kdfull2
          rp = MAX(arho(ind), 1D-90)
          tau(ind) = tau(ind) - (ajt(ind)**2 + 0.25D0*drho(ind)**2)/(rp)
        END DO

      END DO ! loop over directions

! compose localization measure, result stored on 'tau'

      sumpart = 0D0
      DO ind = 1, kdfull2
        rp = MAX(arho(ind), 1D-30)
        tau(ind) = 1D0/(1D0 + (tau(ind)/(tf_fac*rp**(5D0/3D0)))**2)
        IF (tau(ind) .LT. 1D-99) tau(ind) = 0D0
        average_localization(is) = average_localization(is) + tau(ind)*rp
        sumpart = sumpart + rp
      END DO
      average_localization(is) = average_localization(is)/sumpart

      IF (is == 2) THEN
        DO i = 1, kdfull2
          tau(i + kdfull2) = tau(i) ! lower spin in upper storage
        END DO
      END IF

    END DO ! loop over spins

    DEALLOCATE (akk)
    DEALLOCATE (ajt)
    DEALLOCATE (drho)
    DEALLOCATE (arho)

! rescale to final criterion

! preliminary PRINT along axes

    midx = (minx + maxx)/2
    midy = (miny + maxy)/2
    midz = (minz + maxz)/2

#ifdef REALSWITCH
    OPEN (33, STATUS='unknown', FILE='pelfstat.'//outnam)
    WRITE (33, '(/a)') '# localization after static iteration'
#else
    OPEN (33, STATUS='unknown', FILE='pelf.'//outnam)
    WRITE (33, '(/a,1f13.5)') '# localization at time=', tfs
#endif
    WRITE (33, '(a,2(1pg13.5))') '# average localizations:', average_localization
    WRITE (33, '(a)') &
      '# coordinate localiz.-up localiz.-down rho-up rho-down'

    WRITE (33, '(/a)') '# x localization '
    indadd = (midz - 1)*nxyf + (midy - 1)*nyf
    WRITE (33, '(5(1pg13.5))') ((ix - nxsh)*dx, &
                 tau(indadd + ix), tau(indadd + ix + kdfull2), &
                 rho(indadd + ix)*(0.5D0 + rho(indadd + ix + kdfull2)), &
                 rho(indadd + ix)*(0.5D0 - rho(indadd + ix + kdfull2)), &
                                ix=minx, maxx)

    WRITE (33, '(/,/a)') '# y localization '
    indadd = (midz - 1)*nxyf + midx
    WRITE (33, '(5(1pg13.5))') ((iy - nysh)*dy, &
      tau(indadd + (iy - 1)*nyf), tau(indadd + (iy - 1)*nyf + kdfull2), &
      rho(indadd + (iy - 1)*nyf)*(0.5D0 + rho(indadd + (iy - 1)*nyf + kdfull2)), &
      rho(indadd + (iy - 1)*nyf)*(0.5D0 - rho(indadd + (iy - 1)*nyf + kdfull2)), &
                           iy=miny, maxy)

    WRITE (33, '(/,/a)') '# z localization '
    indadd = (midy - 1)*nyf + midx
    WRITE (33, '(5(1pg13.5))') ((iz - nzsh)*dz, &
       tau(indadd + (iz - 1)*nxyf), tau(indadd + (iz - 1)*nxyf + kdfull2), &
       rho(indadd + (iz - 1)*nxyf)*(0.5D0 + rho(indadd + (iz - 1)*nxyf + kdfull2)), &
       rho(indadd + (iz - 1)*nxyf)*(0.5D0 - rho(indadd + (iz - 1)*nxyf + kdfull2)), &
       iz=minz, maxz)
    WRITE (33, '(/)')
    FLUSH (33)
    CLOSE (33)

#ifdef REALSWITCH
    OPEN (33, STATUS='unknown', FILE='pelfstat2Dxy.'//outnam)
    WRITE (33, '(/a)') '# localization after static iteration'
#else
    CALL safeopen(33, it, jelf, 'pelf2Dxy')
    WRITE (33, '(/a,1f13.5)') '# localization at time=', tfs
#endif
    ind = 0
    DO iz = 1, maxz
      z1 = (iz - nzsh)*dz
      DO iy = 1, maxy
        y1 = (iy - nysh)*dy
        DO ix = 1, maxx
          x1 = (ix - nxsh)*dx
          ind = ind + 1
          temp = tau(ind)
          IF (z1 == 0D0) WRITE (33, '(3f10.3,3(1pg15.6))') &
            x1, y1, z1, tau(ind), tau(ind + kdfull2), rho(ind)
        END DO
        IF (z1 == 0D0) WRITE (33, '(1x)')
      END DO
    END DO
    FLUSH (33)
    CLOSE (33)

#ifdef REALSWITCH
    OPEN (33, STATUS='unknown', FILE='pelfstat2Dyz.'//outnam)
    WRITE (33, '(/a)') '# localization after static iteration'
#else
    CALL safeopen(33, it, jelf, 'pelf2Dyz')
    WRITE (33, '(/a,1f13.5)') '# localization at time=', tfs
#endif
    ind = 0
    DO iz = 1, maxz
      z1 = (iz - nzsh)*dz
      DO iy = 1, maxy
        y1 = (iy - nysh)*dy
        DO ix = 1, maxx
          x1 = (ix - nxsh)*dx
          ind = ind + 1
          temp = tau(ind)
          IF (x1 == 0D0) WRITE (33, '(3f10.3,3(1pg15.6))') &
            x1, y1, z1, tau(ind), tau(ind + kdfull2), rho(ind)
        END DO
        IF (x1 == 0D0) WRITE (33, '(1x)')
      END DO
      WRITE (33, '(1x)')
    END DO
    FLUSH (33)
    CLOSE (33)

#ifdef REALSWITCH
    OPEN (33, STATUS='unknown', FILE='pelfstat2Dxz.'//outnam)
    WRITE (33, '(/a)') '# localization after static iteration'
#else
    CALL safeopen(33, it, jelf, 'pelf2Dxz')
    WRITE (33, '(/a,1f13.5)') '# localization at time=', tfs
#endif
    ind = 0
    DO iz = 1, maxz
      z1 = (iz - nzsh)*dz
      DO iy = 1, maxy
        y1 = (iy - nysh)*dy
        DO ix = 1, maxx
          x1 = (ix - nxsh)*dx
          ind = ind + 1
          temp = tau(ind)
          IF (y1 == 0D0) WRITE (33, '(3f10.3,3(1pg15.6))') &
            x1, y1, z1, tau(ind), tau(ind + kdfull2), rho(ind)
        END DO
        IF (y1 == 0D0) WRITE (33, '(1x)')
      END DO
    END DO

    WRITE (33, '(1x)')
    WRITE (33, '(1x)')

    FLUSH (33)
    CLOSE (33)

#ifdef REALSWITCH
    OPEN (33, STATUS='unknown', FILE='pelfstat3D.'//outnam)
    WRITE (33, '(/a)') '# localization after static iteration'
#endif
    DO iz = minz, maxz
      DO iy = minz, maxz
        DO ix = minx, maxx
          indadd = (iy - 1)*nyf + ix + (iz - 1)*nxyf
          WRITE (33, '(1x,3f10.2,3g13.5)') &
            (ix - nxsh)*dx, (iy - nysh)*dy, (iz - nzsh)*dz, &
            tau(indadd), tau(indadd + kdfull2), rho(indadd)
        END DO
        WRITE (33, '(1x)')
      END DO
    END DO
    FLUSH (33)
    CLOSE (33)

    DEALLOCATE (tau)

    RETURN
#ifdef REALSWITCH
  END SUBROUTINE localizer
#else
END SUBROUTINE localize
#endif

