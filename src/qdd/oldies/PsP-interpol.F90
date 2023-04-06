!
! Material for computing PsP on a finer mesh.
! Yet IN testing stage.
!

!-----hpsi -------------------------------------------------------------

SUBROUTINE hpsi(qact, aloc, nbe, itpri)

! ACTION of mean-field Hamiltonian on one s.p. wavefunction:
! qact = wavefunction on which H acts and resulting w.f.
! aloc = local potential for the actual spin component
! ak = kinetic energies IN momentum space
! nbe = NUMBER of state
! itpri = switch for computing s.p. energies (for ABVS(itpri)=1)
! <0 switches TO subtract mean-value of s.p. energy
!
! Version with interpolation TO PsP on fine mesh.

  USE params
  USE util, ONLY: wfovlp
  USE kinetic
  USE twost

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: qact(kdfull2)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN) :: akv(kdfull2)
  INTEGER, INTENT(IN) :: nbe
  INTEGER, INTENT(IN) :: itpri
  COMPLEX(DP), ALLOCATABLE :: qex(:)

! workspaces
  COMPLEX(DP), ALLOCATABLE :: q1(:), q2(:), q1fine(:), q2fine(:), qactfine(:)
  COMPLEX(DP), ALLOCATABLE :: qarray(:, :, :), qarrayfine(:, :, :)
  REAL(DP) :: wnorm
  LOGICAL :: tpri
  LOGICAL, PARAMETER :: tsubmean = .TRUE.
  LOGICAL, PARAMETER :: ttest = .FALSE.
  INTEGER :: i, is, na
  COMPLEX(DP) :: cf

!----------------------------------------------------------------------

  tpri = ABS(itpri) == 1
! WRITE(*,*) 'entering HPSI, tpri=',tpri

  ALLOCATE (q1(kdfull2), q2(kdfull2))
  q1 = (0D0, 0D0)
  q2 = (0D0, 0D0)
! ACTION of kinetic energy

#if(netlib_fft|fftw_cpu)
  CALL fftf(qact, q1)
  DO i = 1, nxyz
    q1(i) = akv(i)*q1(i)
  END DO
  CALL fftback(q1, q2)
!WRITE(*,*) 'FFTW kin'
#else
  CALL ckin3d(qact, q1)
!WRITE(*,*) 'CKIN3D'
#endif
!q2 = -h2m * q1

!STOP ' HPSI not yet appropriate for finite differences'

! ACTION of potential and non-local PsP (optionally)

  IF (ipsptyp == 1) THEN
!~ IF (iswitch_interpol) THEN
!~ ALLOCATE(q1fine(kdfull2fine),q2fine(kdfull2fine),qactfine(kdfull2fine))
!~ ALLOCATE(qarray (nx2,ny2,nz2),qarrayfine (2*nx2-1,2*ny2-1,2*nz2-1))

!~ CALL from1Dto3Dc(qact,qarray,nx2,ny2,nz2) !from coarse vector TO coarse array
!~ CALL interpol3Dc(qarray,qarrayfine) !from coarse array TO fine array
!~ CALL from3Dto1Dc(qactfine,qarrayfine,2*nx2-1,2*ny2-1,2*nz2-1) !from fine array TO fine vector
!~ !
!~ CALL nonlocalcfine(qactfine,q1fine,0)

!~ !
!~ CALL from1Dto3Dc(q1fine,qarrayfine,2*nx2-1,2*ny2-1,2*nz2-1) !from fine vector TO fine array
!~ CALL smoothing3Dc(qarrayfine,qarray) !from fine array TO coarse array
!~ CALL from3Dto1Dc(q1,qarray,nx2,ny2,nz2) !from coarse array TO coarse vector
!~ !
!~ ELSE
    CALL nonlocalc(qact, q1, 0)
!~ END IF ! interpol

    IF (tpri) enonlo(nbe) = wfovlp(qact, q1)
    DO i = 1, nxyz
      q1(i) = q1(i) + qact(i)*aloc(i)
    END DO
  ELSE
    DO i = 1, nxyz
      q1(i) = qact(i)*aloc(i)
    END DO
  END IF

  IF (ifsicp == 5) THEN
    ALLOCATE (qex(kdfull2))
!~ IF(tpri) epotbefore = wfovlp(qact,q1)
    CALL exchg(qact, qex, nbe)
    q1 = q1 + qex
! IF(tpri) THEN
! epotafter = wfovlp(qact,q1)
! WRITE(6,'(a,i4,3(1pg13.5))') ' EXC: nbe,before,after,such=', &
! nbe,epotbefore,epotafter,wfovlp(qact,qex)
! CALL FLUSH(6)
! END IF
    DEALLOCATE (qex)
  END IF

!JM : subtract SIC potential for state NBE
  IF (ifsicp .GE. 8) THEN
    is = ispin(nrel2abs(nbe))
    DO na = 1, nstate
      IF (ispin(nrel2abs(na)) == is) THEN
        cf = wfovlp(psiut(:, na), qact)
        DO i = 1, nxyz
          q1(i) = q1(i) - qnewut(i, na)*cf
        END DO
      END IF
    END DO

  END IF
!JM

  IF (tpri) THEN
    epotsp(nbe) = wfovlp(qact, q1)
    amoy(nbe) = ekinsp(nbe) + epotsp(nbe)
    q2 = q1 + q2
! spvariance(nbe) = SQRT(REAL(wfovlp(q2,q2),DP)-amoy(nbe)**2)
    spvariance(nbe) = SQRT(MAX(REAL(wfovlp(q2, q2), DP) - ABS(wfovlp(qact, q2))**2, 1D-99))
    is = ispin(nrel2abs(nbe))
    IF (ttest) WRITE (*, '(a,2i4,5(1pg13.5))') &
      ' HPSI: nbe,is,esp,var=', nbe, is, amoy(nbe), spvariance(nbe), &
      ekinsp(nbe), epotsp(nbe), REAL(wfovlp(q2, q2), DP)
    CALL flush (6)
! SQRT(REAL(wfovlp(q2,q2))),ABS(wfovlp(qact,q2))
  ELSE
    q2 = q1 + q2
  END IF

  IF (itpri < 0) THEN
    qact = q2 - amoy(nbe)*qact
  ELSE
    qact = q2
  END IF

!espact = amoy(nbe)

  DEALLOCATE (q1, q2)

  RETURN
END SUBROUTINE hpsi

!----------------------------interpol3D(a)--------------------------------------
#ifdef REALSWITCH
SUBROUTINE interpol3Dn(a1, ao)
#else
  SUBROUTINE interpol3Dcn(a1, ao)
#endif

    USE params
    IMPLICIT NONE
! INTEGER :: nx2,ny2,nz2!,nxyz=nx2*ny2*nz2

#ifdef REALSWITCH
!interpolation TARGET array
! REAL(DP):: va(nxyz)
    REAL(DP), INTENT(IN) :: a1(nx2, ny2, nz2)
!interpolation variables
! REAL(DP),ALLOCATABLE ::ai(:,:,:)
    REAL(DP), INTENT(OUT) ::ao(2*nx2, 2*ny2, 2*nz2)
#else
! COMPLEX(DP):: va(nxyz)
    COMPLEX(DP), INTENT(IN) :: a1(nx2, ny2, nz2)
!interpolation variables
! REAL(DP),ALLOCATABLE ::ai(:,:,:)
    COMPLEX(DP), INTENT(OUT) ::ao(2*nx2, 2*ny2, 2*nz2)
#endif
    COMPLEX(dp) :: ai(0:2*nx2, 0:2*ny2, 0:2*nz2)
    COMPLEX(dp) :: a(0:nx2, 0:ny2, 0:nz2)
    INTEGER:: nx2i, ny2i, nz2i
    INTEGER :: ia, ja, ka
    INTEGER :: ia1, ja1, ka1

    ai(:, :, :) = 0.0D0
    a = 0.0
!WRITE (*,*)'hello1'
    nx2i = 2*nx2 - 1
    ny2i = 2*ny2 - 1
    nz2i = 2*nz2 - 1
! ALLOCATE(ai(nx2i,ny2i,nz2i))

!va is a vector no an array so the first step is TO fill the corresponding 3D array a.
!#ifdef REALSWITCH
! CALL from1Dto3D(va,a,nx2,ny2,nz2) !***********
! #else
! CALL from1Dto3D(va,a,nx2,ny2,nz2) !***********
! #endif

    DO ia = 1, nx2
      DO ja = 1, ny2
        DO ka = 1, ny2
          a(ia, ja, ka) = a1(ia, ja, ka)
        END DO
      END DO
    END DO
!-------------------------------------------------------
    DO ia = 1, nx2i
      ia1 = ia/2 + 1
      DO ja = 1, ny2i
        ja1 = ja/2 + 1
        DO ka = 1, nz2i
          ka1 = ka/2 + 1
          ai(ia, ja, ka) = (1.0D0/8.0D0)*(a(ia1 - 1, ja1 - 1, ka1 - 1) + a(ia1 - 1, ja1 - 1, ka1 + 1) + &
                                          a(ia1 - 1, ja1 + 1, ka1 - 1) + a(ia1 - 1, ja1 + 1, ka1 + 1) + &
                                          a(ia1 + 1, ja1 - 1, ka1 - 1) + a(ia1 + 1, ja1 - 1, ka1 + 1) + &
                                          a(ia1 + 1, ja1 + 1, ka1 - 1) + a(ia1 + 1, ja1 + 1, ka1 + 1))
        END DO
      END DO
    END DO
    DO ia = 1, nx2i
      DO ja = 1, ny2i
        DO ka = 1, ny2i
          ao(ia, ja, ka) = ai(ia, ja, ka)
        END DO
      END DO
    END DO
    RETURN
#ifdef REALSWITCH
  END SUBROUTINE interpol3Dn
#else
END SUBROUTINE interpol3Dcn
#endif

!----------------------------interpol3D(a)--------------------------------------
#ifdef REALSWITCH
SUBROUTINE interpol3D(a, ai)
#else
  SUBROUTINE interpol3Dc(a, ai)
#endif
! 3D interpolation by a linear method
! with 1 interpolated point between 2 REAL point of each DIMENSION
! THEN there is 1 interpolated IN a center of each elementary edge (2 REAL points),
! of each elementary face (4 REAL points) and of each elementary square (8 REAL points)

! This SUBROUTINE interpolates a nx2.ny2.nz2 DIMENSION array TO (2*nx2-1)(2*ny2-1)(2*nz2-1) array

! illustration of the method IN 1D
!
! o--.--o--.--o--.--o--.--o--.--o--.--o--.--o space discretisation
!
! f: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 interpolated vector
! F: 1 2 3 4 5 6 7 8 REAL vector
!
! f(n)=F(N) with an odd n=2N-1 (odd interpolated points correspond TO REAL points)
! f(n)=(1/2)[f(n-1)+f(n+1)] with an even n

    USE params
    IMPLICIT NONE
! INTEGER :: nx2,ny2,nz2!,nxyz=nx2*ny2*nz2

#ifdef REALSWITCH
!interpolation TARGET array
! REAL(DP):: va(nxyz)
    REAL(DP), INTENT(IN) :: a(nx2, ny2, nz2)
!interpolation variables
! REAL(DP),ALLOCATABLE ::ai(:,:,:)
    REAL(DP), INTENT(OUT) ::ai(2*nx2, 2*ny2, 2*nz2)
#else
! COMPLEX(DP):: va(nxyz)
    COMPLEX(DP), INTENT(IN) :: a(nx2, ny2, nz2)
!interpolation variables
! REAL(DP),ALLOCATABLE ::ai(:,:,:)
    COMPLEX(DP), INTENT(OUT) ::ai(2*nx2, 2*ny2, 2*nz2)
#endif
    INTEGER:: nx2i, ny2i, nz2i
    INTEGER :: ia, ja, ka

    ai(:, :, :) = 0.0D0
!WRITE (*,*)'hello1'
    nx2i = 2*nx2 - 1
    ny2i = 2*ny2 - 1
    nz2i = 2*nz2 - 1
! ALLOCATE(ai(nx2i,ny2i,nz2i))

!va is a vector no an array so the first step is TO fill the corresponding 3D array a.
!#ifdef REALSWITCH
! CALL from1Dto3D(va,a,nx2,ny2,nz2) !***********
! #else
! CALL from1Dto3D(va,a,nx2,ny2,nz2) !***********
! #endif

!-------------------------------------------------------
!filling the odd points of interpolated array with corresponding points of TARGET array
    DO ia = 1, nx2
      DO ja = 1, ny2
        DO ka = 1, nz2
          ai(2*ia - 1, 2*ja - 1, 2*ka - 1) = a(ia, ja, ka)
        END DO
      END DO
    END DO

!*******************************
!filling the even interpolated points IN the center of elementary squares
    DO ia = 2, nx2i, 2
      DO ja = 2, ny2i, 2
        DO ka = 2, nz2i, 2
          ai(ia, ja, ka) = (1.0D0/8.0D0)*(ai(ia - 1, ja - 1, ka - 1) + ai(ia - 1, ja - 1, ka + 1) + &
                                          ai(ia - 1, ja + 1, ka - 1) + ai(ia - 1, ja + 1, ka + 1) + &
                                          ai(ia + 1, ja - 1, ka - 1) + ai(ia + 1, ja - 1, ka + 1) + &
                                          ai(ia + 1, ja + 1, ka - 1) + ai(ia + 1, ja + 1, ka + 1))
        END DO
      END DO
    END DO

!*******************************
!filling the even interpolated points IN the center of elementary faces
    DO ia = 2, nx2i, 2
      DO ja = 2, ny2i, 2
        DO ka = 2, nz2i, 2
          ai(ia, ja, ka) = (1.0D0/4.0D0)*(ai(ia, ja - 1, ka - 1) + ai(ia, ja - 1, ka + 1) + &
                                          ai(ia, ja + 1, ka - 1) + ai(ia, ja + 1, ka + 1))
        END DO
      END DO
    END DO

    DO ja = 2, ny2i, 2
      DO ka = 2, nz2i, 2
        DO ia = 2, nx2i, 2
          ai(ia, ja, ka) = (1.0D0/4.0D0)*(ai(ia - 1, ja, ka - 1) + ai(ia - 1, ja, ka + 1) + &
                                          ai(ia + 1, ja, ka - 1) + ai(ia + 1, ja, ka + 1))
        END DO
      END DO
    END DO

    DO ka = 2, nz2i, 2
      DO ia = 2, nx2i, 2
        DO ja = 2, ny2i, 2
          ai(ia, ja, ka) = (1.0D0/4.0D0)*(ai(ia - 1, ja - 1, ka) + ai(ia - 1, ja + 1, ka) + &
                                          ai(ia + 1, ja - 1, ka) + ai(ia + 1, ja + 1, ka))
        END DO
      END DO
    END DO

!*******************************
!filling the even interpolated points IN the center of elementary edges
    DO ja = 2, ny2i, 2
      DO ka = 2, nz2i, 2
        DO ia = 2, nx2i, 2
          ai(ia, ja, ka) = (1.0D0/2.0D0)*(ai(ia - 1, ja, ka) + ai(ia + 1, ja, ka))
        END DO
      END DO
    END DO

    DO ka = 2, nz2i, 2
      DO ia = 2, nx2i, 2
        DO ja = 2, ny2i, 2
          ai(ia, ja, ka) = (1.0D0/2.0D0)*(ai(ia, ja - 1, ka) + ai(ia, ja + 1, ka))
        END DO
      END DO
    END DO

    DO ia = 2, nx2i, 2
      DO ja = 2, ny2i, 2
        DO ka = 2, nz2i, 2
          ai(ia, ja, ka) = (1.0D0/2.0D0)*(ai(ia, ja, ka - 1) + ai(ia, ja, ka + 1))
        END DO
      END DO
    END DO

!-------------------------------------------------------

!WRITE (*,*)'END interpol'

    RETURN
#ifdef REALSWITCH
  END SUBROUTINE interpol3D
#else
END SUBROUTINE interpol3Dc
#endif

!----------------------------smoothing3D(a)--------------------------------------
#ifdef REALSWITCH
SUBROUTINE smoothing3D(a, as)
#else
  SUBROUTINE smoothing3Dc(a, as)
#endif

!the smoothing method consists on the average of the 27 points. The 26 points around the REAL point count for the half of the contribution and the main point for the other half.
!
!how the 27 points are counted? e.g. IN 2D, at z(k)
! y(j+1) . . .
!
! y(j) . o .
!
! y(j-1) . . .
!
! x(i+1) x(i) x(i-1)

! o is the REAL main point at x(i),y(j),z(k)
! . are the interpolated points so there are 8 at the z(k) plan, there are also 9 on each z(k-1) and z(k+1) plans

    USE params
    IMPLICIT NONE
#ifdef REALSWITCH
!smothing array
    REAL(DP), INTENT(OUT) :: as(nx2, ny2, nz2)
!smoothing TARGET array
    REAL(DP), INTENT(IN) ::a(2*nx2, 2*ny2, 2*nz2)

#else
    COMPLEX(DP), INTENT(OUT) :: as(nx2, ny2, nz2)
!smoothing TARGET array
    COMPLEX(DP), INTENT(IN) ::a(2*nx2, 2*ny2, 2*nz2)
#endif
    INTEGER :: ia, ja, ka, i1, i2, i3
    REAL(DP) :: fac
    as = 0.0

!*******************************
!DO ia=1,nx2
! IF ((ia==1).or.(ia==nx2))THEN ! for the orthogonal border plans TO the x axis
! DO ja=1,ny2
! DO ka=1,nz2
! as(ia,ja,ka)=a(2*ia-1,2*ja-1,2*ka-1)
! END DO
! END DO
! END IF
!END DO

    as(1, 1:ny2, 1:nz2) = a(1, 1:2*ny2:2, 1:2*nz2:2)
    as(nx2, 1:ny2, 1:nz2) = a(nx2, 1:2*ny2:2, 1:2*nz2:2)

!DO ja=1,ny2
! IF ((ja==1).or.(ja==ny2))THEN ! for the orthogonal border plans TO the y axis
! DO ka=1,nz2
! DO ia=1,nx2
! ! IF((ia.ne.1).and.(ia.ne.nx2))THEN
! as(ia,ja,ka)=a(2*ia-1,2*ja-1,2*ka-1)
! ! END IF
! END DO
! END DO
! END IF
!END DO

    as(1:nx2, 1, 1:nz2) = a(1:2*nx2:2, 1, 1:2*nz2:2)
    as(1:ny2, ny2, 1:nz2) = a(1:2*ny2:2, ny2, 1:2*nz2:2)

!DO ka=1,nz2
! IF ((ka==1).or.(ka==nz2))THEN ! for the orthogonal border plans TO the z axis
! DO ia=1,nx2
! !IF((ia.ne.1).and.(ia.ne.nx2))THEN
! DO ja=1,ny2
! IF((ja.ne.1).and.(ja.ne.ny2)) as(ia,ja,ka)=a(2*ia-1,2*ja-1,2*ka-1)
! END DO
! !END IF
! END DO
! END IF
!END DO

    as(1:nx2, 1:ny2, 1) = a(1:2*nx2:2, 1:2*ny2:2, nz2)
    as(1:nx2, 1:ny2, nz2) = a(1:nx2:2, 1:2*ny2:2, nz2)

!*******************************
!for the inner points
    DO ia = 2, nx2 - 1
      DO ja = 2, ny2 - 1
        DO ka = 2, nz2 - 1
          DO i1 = -2, 0
            DO i2 = -2, 0
              DO i3 = -2, 0
                fac = 1.0D0/27.D0
                IF (i1 .eq. -1) THEN
                  IF (i2 .eq. -1) THEN
                    IF (i3 .eq. -1) THEN
                      fac = 26.0D0/27.0D0
                    END IF
                  END IF
                END IF
                as(ia, ja, ka) = as(ia, ja, ka) + fac*a(2*ia + i1, 2*ja + i2, 2*ka + i3)
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

! todo : Compare speed with this FORM :
!~ !REAL(DP):: fac(-2:1, -2:1, -2:1)

!~ fac = 1.0D0/27.D0
!~ fac(-1,-1,-1)= 26.0D0/27.0D0

!~ DO ka=2,nz2-1
!~ DO ja=2,ny2-1
!~ DO ia=2,nx2-1
!~ DO i3=-2,0
!~ DO i2=-2,0
!~ DO i1=-2,0
!~ as(ia,ja,ka)=as(ia,ja,ka)+ fac(i1,i2,i3) * a(2*ia+i1,2*ja+i2,2*ka+i3)
!~ END DO
!~ END DO
!~ END DO
!~ END DO
!~ END DO
!~ END DO

    RETURN
#ifdef REALSWITCH
  END SUBROUTINE smoothing3D
#else
END SUBROUTINE smoothing3Dc
#endif
