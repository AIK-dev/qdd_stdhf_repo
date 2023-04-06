#if(raregas)
!-----analyze_surf------------------------------------------------

SUBROUTINE analyze_surf(it)

! OPEN protocol files and PRINT properties of substrate

  USE params
  USE util, ONLY: gettemperature, safeopen

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: it

  INTEGER :: i, ii, ion, nimobc, nimobk
  REAL(DP) :: sumcx, sumcy, sumcz, sumkx, sumky, sumkz
  REAL(DP) :: ekx, eky, ekz, ekcx, ekcy, ekcz, ekvx, ekvy, ekvz
  REAL(DP) :: amfac1, amfac2
  REAL(DP) :: rho(2*kdfull2)

!----------------------------------------------------------------

! tfs=it*dt1*0.0484
  IF (nc + NE + nk > 0) THEN
    CALL safeopen(24, it, jpos, 'pposcore')
    CALL safeopen(157, it, jpos, 'ptempsurf')
    CALL safeopen(25, it, jpos, 'pposdip')
    CALL safeopen(125, it, jpos, 'pposkat')
    CALL safeopen(424, it, jpos, 'penerSurf')
    CALL safeopen(26, it, jvel, 'pvelcore')
    CALL safeopen(126, it, jvel, 'pvelkat')
    IF (ifadiadip /= 1) CALL safeopen(27, it, jvel, 'pveldip')
    CALL safeopen(152, it, jvel, 'pkinencore')
    CALL safeopen(153, it, jvel, 'pkinenkat')
    IF (((jpos > 0 .AND. MOD(it, jpos) == 0) &
         .OR. (jvel > 0 .AND. MOD(it, jvel) == 0))) THEN

      IF (MOD(it, jpos) == 0) THEN
        CALL getsurfprops
        WRITE (424, '(f12.4,4e15.5)') tfs, ekincsurf, ekinesurf, ekinksurf, ekinsurf
      END IF

      IF (MOD(it, jsavesurf) == 0) THEN
        OPEN (464, STATUS='unknown', FILE='for006surf.'//outnam)
        DO i = 1, nc
          WRITE (464, '(6e17.7,i6)') xc(i), yc(i), zc(i), &
            xe(i), ye(i), ze(i), imobc(i), imobe(i)
        END DO
        DO i = 1, nk
          WRITE (464, '(3e17.7,i6)') xk(i), yk(i), zk(i), imobk(i)
        END DO
        CLOSE (464)
      END IF

      sumcx = 0D0
      sumcy = 0D0
      sumcz = 0D0
      nimobc = 0
      DO ion = 1, nc
        IF (MOD(it, jpos) == 0) THEN
          IF (iprintonlyifmob == 0 .OR. imobc(ion) /= 0) &
            WRITE (24, '(1f13.5,3e17.8)') tfs, xc(ion), yc(ion), zc(ion)
          IF (iprintonlyifmob == 0 .OR. imobe(ion) /= 0) &
            WRITE (25, '(1f13.5,3e17.8)') tfs, xe(ion), ye(ion), ze(ion)
        END IF

        IF (MOD(it, jvel) == 0) THEN
          IF (iprintonlyifmob == 0 .OR. imobc(ion) /= 0) &
            WRITE (26, '(1f13.5,3e17.8)') tfs, pxc(ion), pyc(ion), pzc(ion)
          IF (ifadiadip /= 1 .AND. &
              (iprintonlyifmob == 0 .OR. imobe(ion) /= 0)) THEN
            WRITE (27, '(1f13.5,3e17.8)') tfs, pxe(ion), pye(ion), pze(ion)
            nimobc = nimobc + 1
          END IF
        END IF

        IF (imobc(ion) /= 0) THEN
          sumcx = sumcx + (pxc(ion)**2)/mion/1836D0
          sumcy = sumcy + (pyc(ion)**2)/mion/1836D0
          sumcz = sumcx + (pzc(ion)**2)/mion/1836D0
        END IF

      END DO

      WRITE (6, *) mion, mkat
      IF (MOD(it, jvel) == 0) THEN
        WRITE (152, '(1f13.5,4e17.8,i6)') &
          tfs, sumcx, sumcy, sumcz, sumcx + sumcy + sumcz, nimobc
        CALL gettemperature(1)
      END IF

      sumkx = 0D0
      sumky = 0D0
      sumkz = 0D0
      nimobk = 0
      DO ion = 1, nk
        IF (MOD(it, jpos) == 0 .AND. &
            (iprintonlyifmob == 0 .OR. imobk(ion) /= 0)) &
          WRITE (125, '(1f13.5,3e17.8)') tfs, xk(ion), yk(ion), zk(ion)
        IF (MOD(it, jvel) == 0 .AND. &
            (iprintonlyifmob == 0 .OR. imobk(ion) /= 0)) THEN
          WRITE (126, '(1f13.5,3e17.8)') tfs, pxk(ion), pyk(ion), pzk(ion)
          nimobk = nimobk + 1
        END IF

        IF (imobk(ion) /= 0) THEN
          sumkx = sumkx + (pxk(ion)**2)/mkat/1836D0
          sumky = sumky + (pyk(ion)**2)/mkat/1836D0
          sumkz = sumkx + (pzk(ion)**2)/mkat/1836D0
        END IF
      END DO

      IF (MOD(it, jvel) == 0) WRITE (153, '(1f13.5,4e17.8,i6)') &
        tfs, sumkx, sumky, sumkz, nimobk

      CALL flush (24)
      CALL flush (25)
      CALL flush (26)
      IF (ifadiadip /= 1) CALL flush (27)
    END IF
  END IF

  IF (myn == 0 .AND. jener > 0 .AND. MOD(it, jener) == 0 .AND. nrare > 0) THEN
! matrix
    ekcx = 0D0
    ekcy = 0D0
    ekcz = 0D0
    ekvx = 0D0
    ekvy = 0D0
    ekvz = 0D0
    amfac1 = amu(np(1))*1836D0*ame*2D0
    amfac2 = amu(np(nrare + 1))*1836D0*ame*2D0
    DO ion = 1, nrare
      ekcx = ekcx + pxc(ion)*pxc(ion)
      ekcy = ekcy + pyc(ion)*pyc(ion)
      ekcz = ekcz + pzc(ion)*pzc(ion)
      ekvx = ekvx + pxe(ion)*pxe(ion)
      ekvy = ekvy + pye(ion)*pye(ion)
      ekvz = ekvz + pze(ion)*pze(ion)
    END DO
    ekx = ekcx/amfac1 + ekvx/amfac2
    eky = ekcy/amfac1 + ekvy/amfac2
    ekz = ekcz/amfac1 + ekvz/amfac2
    CALL safeopen(35, it, jener, 'penermat')
    WRITE (35, '(4f13.5)') tfs, ekx, eky, ekz, etot
    CALL flush (35)
  END IF

  IF (isurf == 1 .AND. iforcecl2co /= 0 .AND. MOD(it, iforcecl2co) == 0) THEN

    CALL getforces_clust2cores(rho, 0)

    OPEN (256, STATUS='unknown', FILE='force_clust2cores.'//outnam)
    WRITE (256, '(a,f10.5)') 'core pos val pos force on cores at t=', tfs

    DO ii = 1, nc

      WRITE (256, '(9e20.10)') xc(ii), yc(ii), zc(ii), xe(ii), ye(ii), &
        ze(ii), fxc(ii), fyc(ii), fzc(ii)

    END DO

    CLOSE (256)
  END IF

  RETURN
END SUBROUTINE analyze_surf
#endif

#if(raregas)
!-----vstep--------------------------------------------------------

SUBROUTINE vstep(rho, psi, it, dt)

! propagation of rare gas valence clouds by leapfrog method

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
  INTEGER, INTENT(IN) :: it
  REAL(DP), INTENT(IN) :: dt

  INTEGER :: i
  REAL(DP), ALLOCATABLE :: xm(:)

!------------------------------------------------------------------

! optionally freeze valence momenta

  DO i = 1, NE
    IF (imobe(i) == 0) THEN
      pxe(i) = 0D0
      pye(i) = 0D0
      pze(i) = 0D0
    END IF
  END DO

! propagation of positions first

! xm=amu(np(nrare+1))*1836.0*ame
! CALL leapfr(cx(nrare+1),cy(nrare+1),cz(nrare+1),
! & cpx(nrare+1),cpy(nrare+1),cpz(nrare+1),dt,xm,nrare)

  ALLOCATE (xm(1:ne))
  xm = amu(-18)*1836.0D0*ame
  CALL leapfr(xe(1), ye(1), ze(1), pxe(1), pye(1), pze(1), dt, xm, ne, 2)

! update subgrids IN CASE of pseudo-densities

  IF (ipseudo) THEN
    CALL updatesubgrids
  END IF

! THEN propagation of momenta

  CALL getforces(rho, psi, it, 0) ! forces on valences with new positions

  xm = 1D0
  CALL leapfr(pxe(1), pye(1), pze(1), fxe(1), fye(1), fze(1), dt, xm, ne, 2)
  DEALLOCATE (xm)

  RETURN
END SUBROUTINE vstep
!-----vstep--------------------------------------------------------

SUBROUTINE vstepv(rho, psi, it, dt)

! propagation of rare gas valence clouds by velocity Verlet

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
  INTEGER, INTENT(IN) :: it
  REAL(DP), INTENT(IN) :: dt

  INTEGER :: i
  REAL(DP), ALLOCATABLE :: xm(:)

!------------------------------------------------------------------

! optionally freeze valence momenta

  DO i = 1, NE
    IF (imobe(i) == 0) THEN
      pxe(i) = 0D0
      pye(i) = 0D0
      pze(i) = 0D0
    END IF
  END DO

! propagation of positions first

! xm=amu(np(nrare+1))*1836.0*ame

  ALLOCATE (xm(1:ne))
  xm = amu(-18)*1836D0*ame
  CALL velverlet1(xe(1), ye(1), ze(1), pxe(1), pye(1), pze(1), &
                  fxe(1), fye(1), fze(1), dt, xm, ne, 2)
  DEALLOCATE (xm)

! update subgrids IN CASE of pseudo-densities

  IF (ipseudo) THEN
    CALL updatesubgrids
  END IF

! forces on valences with new positions

  CALL getforces(rho, psi, it, 0)

! THEN propagation of momenta

  CALL velverlet2(pxe(1), pye(1), pze(1), fxe(1), fye(1), fze(1), dt, ne, 2)

  RETURN
END SUBROUTINE vstepv
#endif

! From former 'analse.F90'

SUBROUTINE getsurfprops
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! calculate kinetic energy of surface ions
#if(raregas)
  INTEGER :: i

  ekincsurf = 0D0
  ekinesurf = 0D0
  ekinksurf = 0D0
  ekinsurf = 0D0

  DO i = 1, nc
    ekincsurf = ekincsurf + (pxc(i)*pxc(i) + pyc(i)*pyc(i) + pzc(i)*pzc(i))/ &
                2D0/mion/1836D0/ame
  END DO
  DO i = 1, NE
    ekinesurf = ekinesurf + (pxe(i)*pxe(i) + pye(i)*pye(i) + pze(i)*pze(i))/ &
                2D0/me/1836D0/ame
  END DO
  DO i = 1, NE
    ekinksurf = ekinksurf + (pxk(i)*pxk(i) + pyk(i)*pyk(i) + pzk(i)*pzk(i))/ &
                2D0/mkat/1836D0/ame
  END DO

  ekinsurf = ekincsurf + ekinesurf + ekinksurf
#endif

  RETURN
END SUBROUTINE getsurfprops
!------------------------------------------------------------
