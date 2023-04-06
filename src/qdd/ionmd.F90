
!------------------------------------------------------------
#if(extended)
SUBROUTINE heatcluster(totemp)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
  REAL(DP), INTENT(IN):: totemp
! This routine heats up the cluster.
! The momenta of the cluster ions
! are reset to values such that the kinetic energy
! corresponds to an instantaneous "thermal" energy 'toTemp'.
! 'toTemp' is in Kelvin!!

  RETURN
END SUBROUTINE heatcluster
!------------------------------------------------------------
#endif

!-----reset_ions--------------------------------------------

SUBROUTINE reset_ions()

! Resets all ionic momenta to zero, is used in cooling.

  USE params
  IMPLICIT NONE

  INTEGER :: ion

!-----------------------------------------------------------

  IF ((ekionold - ekion) > 0D0) THEN
    jekion = jekion + 1
    DO ion = 1, nion
      cpx(ion) = 0D0
      cpy(ion) = 0D0
      cpz(ion) = 0D0
    END DO
#if(raregas)
    DO ion = 1, nc
      pxc(ion) = 0D0
      pyc(ion) = 0D0
      pzc(ion) = 0D0
      pxe(ion) = 0D0
      pye(ion) = 0D0
      pze(ion) = 0D0
    END DO
    DO ion = 1, nk
      pxk(ion) = 0D0
      pyk(ion) = 0D0
      pzk(ion) = 0D0
    END DO
#endif
    icooltimes = icooltimes + 1
    WRITE (6, *) 'jekion=', jekion
    ekionold = 0D0
  END IF
  IF (ekionold < ekion) ekionold = ekion

  RETURN
END SUBROUTINE reset_ions

!-----energkin_ions--------------------------------------------

REAL(DP) FUNCTION enerkin_ions()

! Kinetic energy of ions

  USE params
  IMPLICIT NONE

  INTEGER :: ion
  REAL(DP) :: ek, sumekinion

! calculate kinetic energy of ions

  sumekinion = 0D0
  DO ion = 1, nion
    ek = cpx(ion)*cpx(ion) + cpy(ion)*cpy(ion) + cpz(ion)*cpz(ion)
    ek = ek/(2D0*1836D0*amu(np(ion))*ame)
    sumekinion = sumekinion + ek
  END DO
#if(raregas)
  DO ion = 1, nc
    ek = pxc(ion)*pxc(ion)
    ek = ek + pyc(ion)*pyc(ion)
    ek = ek + pzc(ion)*pzc(ion)
    ek = ek/(2D0*1836D0*mion*ame)
    sumekinion = sumekinion + ek
  END DO

  DO ion = 1, NE
    ek = pxe(ion)*pxe(ion)
    ek = ek + pye(ion)*pye(ion)
    ek = ek + pze(ion)*pze(ion)
    ek = ek/(2D0*1836D0*me*ame)
    sumekinion = sumekinion + ek
  END DO

  DO ion = 1, nk
    ek = pxk(ion)*pxk(ion)
    ek = ek + pyk(ion)*pyk(ion)
    ek = ek + pzk(ion)*pzk(ion)
    ek = ek/(2D0*1836D0*mkat*ame)
    sumekinion = sumekinion + ek
  END DO
#endif
  enerkin_ions = sumekinion

  RETURN
END FUNCTION enerkin_ions

!-----energy_ions--------------------------------------------

REAL(DP) FUNCTION energ_ions()

! Potential energy of ions

  USE params
  IMPLICIT NONE

  INTEGER :: ion, ion1, ion2
  REAL(DP):: dist, dist2, sumion
  INTEGER, EXTERNAL :: iptyp
  REAL(DP), EXTERNAL :: v_ion_ion

#if(raregas)
  INTEGER :: ii, ii2, ion11, jj, jj2
  REAL(DP) :: ccc1, ccc2, sss, sss1, sss2, xtmp, ytmp, ztmp
  INTEGER, EXTERNAL :: iconvlongtoshort
  REAL(DP), EXTERNAL :: sigsig
  REAL(DP), EXTERNAL :: v_soft
  REAL(DP), EXTERNAL :: getcharge
  REAL(DP), EXTERNAL :: getdistance2
  REAL(DP), EXTERNAL :: getmixedwidth
#endif

  IF (nion2 == 2) THEN
    energ_ions = 0D0
    RETURN
  END IF

  sumion = 0.0D0
  enii = 0.0D0
  enig = 0.0D0
  engg = 0.0D0

  IF (ipsptyp == 0) THEN
    IF (nion2 /= 0) THEN
      sumion = 0D0
      DO ion = 2, nion
        DO ion1 = 1, ion - 1
          dist2 = (cx(ion) - cx(ion1))**2 + (cy(ion) - cy(ion1))**2 &
                  + (cz(ion) - cz(ion1))**2
          dist = SQRT(dist2)
          sumion = sumion + v_ion_ion(dist, ion, ion1)
        END DO
      END DO
    END IF
  ELSE IF (ipsptyp >= 1) THEN
    IF (nion2 /= 0) THEN
      sumion = 0D0
      DO ion1 = 2, nion
        DO ion2 = 1, ion1 - 1
          dist2 = (cx(ion1) - cx(ion2))**2 + (cy(ion1) - cy(ion2))**2 &
                  + (cz(ion1) - cz(ion2))**2
          dist = SQRT(dist2)
          sumion = sumion + e2*ch(np(ion1))*ch(np(ion2))/dist
        END DO
      END DO

    END IF
  ELSE
    STOP 'this sort of PsP not yet implemented'
  END IF

#if(raregas)
  IF (idielec /= 0) THEN

! each ion interacts with each image ion!!

    IF (ipsptyp == 0) THEN
      IF (nion2 /= 0) THEN
        DO ion = 1, nion
          DO ion1 = 1, nion
            dist2 = (cx(ion) + cx(ion1) - 2D0*xdielec)**2 + (cy(ion) - cy(ion1))**2 &
                    + (cz(ion) - cz(ion1))**2
            dist = SQRT(dist2)
            sumion = sumion - (epsdi - 1)/(epsdi + 1)*v_ion_ion(dist, ion, ion1)
          END DO
        END DO
      END IF
    ELSE IF (ipsptyp >= 1) THEN
      IF (nion2 /= 0) THEN
        DO ion1 = 1, nion
          DO ion2 = 1, nion
            dist2 = (cx(ion1) + cx(ion2) - 2D0*xdielec)**2 + (cy(ion1) - cy(ion2))**2 &
                    + (cz(ion1) - cz(ion2))**2
            dist = SQRT(dist2)
            sumion = sumion - (epsdi - 1)/(epsdi + 1)* &
                     e2*ch(np(ion1))*ch(np(ion2))/dist
          END DO
        END DO

      END IF
    ELSE
      STOP 'this sort of PsP not yet implemented'
    END IF

  END IF
#endif

  enii = sumion

#if(raregas)
! now the ion-GSM energy
  IF (isurf /= 0) THEN
    IF (ipsptyp == 0) THEN

      DO ion = 1, nion

        DO ion1 = 1, nc + NE + nk

          CALL getparas(ion1)

          dist2 = (rvectmp(1) - cx(ion))**2 + (rvectmp(2) - cy(ion))**2 &
                  + (rvectmp(3) - cz(ion))**2
          dist = SQRT(dist2)
          dist = MAX(dist, small)

          sss1 = sigsig(sigtmp, sgm1(np(ion)))
          sss2 = sigsig(sigtmp, sgm2(np(ion)))

          !WRITE(*,*) ' ion1 etc:',ion,ion1,dist, &
          ! e2*chgtmp*chg1(np(ion))* v_soft(dist,sss1*sq2), &
          ! e2*chgtmp*chg2(np(ion))* v_soft(dist,sss2*sq2),sumion

          sumion = sumion + e2*chgtmp*chg1(np(ion))*v_soft(dist, sss1*sq2)
          sumion = sumion + e2*chgtmp*chg2(np(ion))*v_soft(dist, sss2*sq2)

          !WRITE(*,*) ' ion1 etc:',ion,ion1,dist, &
          ! e2*chgtmp*chg1(np(ion))* v_soft(dist,sss1*sq2), &
          ! e2*chgtmp*chg2(np(ion))* v_soft(dist,sss2*sq2),sumion

          IF (idielec /= 0) THEN

            dist2 = (-rvectmp(1) + 2D0*xdielec - cx(ion))**2 + &
                    (rvectmp(2) - cy(ion))**2 + (rvectmp(3) - cz(ion))**2
            dist = SQRT(dist2)
            dist = MAX(dist, small)

            sumion = sumion - (epsdi - 1)/(epsdi + 1)*e2*chgtmp*chg1(np(ion))* &
                     v_soft(dist, sss1*sq2)
            sumion = sumion - (epsdi - 1)/(epsdi + 1)*e2*chgtmp*chg2(np(ion))* &
                     v_soft(dist, sss2*sq2)

          END IF

          ion11 = iconvlongtoshort(ion1)
          CALL getshortenergy(iptyp(ion1), 4, ion11, ion)

          !WRITE(*,*) ' ion11:',ion11,rscaltmp
          sumion = sumion + rscaltmp

        END DO

      END DO

    ELSE

      STOP 'Goedecker still to be re-implemented with isurf.ne.1'

    END IF

    enig = sumion - enii

! now the GSM-GSM

    DO ii = 1, nc + NE + nk - 1

      ccc1 = getcharge(ii)

      DO jj = ii + 1, nc + NE + nk

        dist2 = getdistance2(ii, jj)

        IF (ii <= nc .AND. jj == ii + nc) THEN ! spring force
          sumion = sumion + 0.5D0*cspr*dist2
          ! WRITE(*,*) ' cspr:',ii,jj,0.5D0*cspr*dist2
        ELSE ! Coulombian erf force

          dist = SQRT(dist2)

          ccc2 = getcharge(jj)
          sss = getmixedwidth(ii, jj)

! electrostatic part
          sumion = sumion + ccc1*ccc2*e2*v_soft(dist, sss*sq2)

          ii2 = iconvlongtoshort(ii)
          jj2 = iconvlongtoshort(jj)

          CALL getshortenergy(iptyp(ii), iptyp(jj), ii2, jj2)

          sumion = sumion + rscaltmp
          !WRITE(*,*) ' iijj:',ii,jj,ccc1,ccc2,dist, &
          ! ccc1*ccc2*e2*v_soft(dist,sss*sq2),rscaltmp

        END IF

! WRITE(6,'(2a,3f15.8)') 'Achtung############sumion=',
! & 'rScalTmp= dist^2=',
! & sumion,rScalTmp,dist2
      END DO
    END DO

    WRITE (*, *) ' GSM: sumion=', sumion, idielec

    IF (idielec /= 0) THEN

      DO ii = 1, nc + NE + nk

        CALL getparas(ii)

        ccc1 = chgtmp
        xtmp = rvectmp(1)
        ytmp = rvectmp(2)
        ztmp = rvectmp(3)

        DO jj = 1, nc + NE + nk

          CALL getparas(jj)

          dist2 = (xtmp + rvectmp(1) - 2D0*xdielec)**2 + (ytmp - rvectmp(2))**2 + &
                  (ztmp - rvectmp(3))**2

          dist = dist2**0.5D0

          ccc2 = chgtmp
          sss = getmixedwidth(ii, jj)

! electrostatic part
          sumion = sumion - (epsdi - 1D0)/(epsdi + 1D0)* &
                   ccc1*ccc2*e2*v_soft(dist, sss*sq2)

        END DO
      END DO
    END IF

  END IF ! isurf eq 0
#endif

  engg = sumion - enii - enig
  energ_ions = sumion

  RETURN
END FUNCTION energ_ions
!------------------------------------------------------------

!-----itstep-------------------------------------------------

SUBROUTINE itstep(rho, it, psi)

! Ionic time step, propagates ionic cores.
!
! Input:
! rho = local electron density
! psi = set of s.p. wavefunctions
! it = time step number in calling routine
! ionic positions and velocities communicated via module 'params'

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: it
  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)

  REAL(DP), ALLOCATABLE :: xm(:) ! array for actual particle mass
  REAL(DP), ALLOCATABLE :: tfac(:) ! array acceleration factors

  LOGICAL, PARAMETER :: tsmooth = .true. ! switch to smooth termination of acceleration
  INTEGER :: i
  REAL(DP) :: dtaccel, xcmion, ycmion, zcmion
  REAL(DP) :: fxcm, fycm, fzcm
#if(raregas)
  REAL(DP) :: pxcmion, pycmion, pzcmion
#endif

#if(extended)
! what time is it ?
  tfs = it*dt1*0.0484D0/(2D0*ame)
#endif

#if(raregas)
  IF (isurf /= 0) THEN
    DO i = 1, nc
      IF (imobc(i) == 0) THEN
        pxc(i) = 0D0
        pyc(i) = 0D0
        pzc(i) = 0D0
      END IF
    END DO
    DO i = 1, nk
      IF (imobk(i) == 0) THEN
        pxk(i) = 0D0
        pyk(i) = 0D0
        pzk(i) = 0D0
      END IF
    END DO
  END IF
#endif

! propagation of positions

#if(raregas)
  IF (isurf /= 0 .AND. nc > 0) THEN
    IF (ALLOCATED(xm)) DEALLOCATE (xm)
    ALLOCATE (xm(1:nc))
    xm = mion*1836D0*ame
    CALL leapfr(xc(1:nc), yc(1:nc), zc(1:nc), pxc(1:nc), pyc(1:nc), pzc(1:nc), dt1, xm, nc, 1)
    DEALLOCATE (xm)
  END IF

  IF (nk > 0) THEN
! correction of velocities to fox c.m.
    IF (tfixcmion) THEN
      pxcmion = SUM(pxk(1:nion))/nion
      pycmion = SUM(pyk(1:nion))/nion
      pzcmion = SUM(pzk(1:nion))/nion
      pxk = pxk - pxcmion
      pyk = pyk - pycmion
      pzk = pzk - pzcmion
    END IF
    IF (ALLOCATED(xm)) DEALLOCATE (xm)
    ALLOCATE (xm(1:nk))
    xm = mkat*1836D0*ame
    CALL leapfr(xk(1:nk), yk(1:nk), zk(1:nk), pxk(1:nk), pyk(1:nk), pzk(1:nk), dt1, xm, nk, 3)
    DEALLOCATE (xm)
  END IF
#endif

  IF (nion > 0 .AND. ionmdtyp /= 0) THEN

! compute c.m. of ions
    IF (tfixcmion) THEN
      xcmion = SUM(cx(1:nion))/nion
      ycmion = SUM(cy(1:nion))/nion
      zcmion = SUM(cz(1:nion))/nion
    END IF

! propagation of molecule ions
    IF (ALLOCATED(xm)) DEALLOCATE (xm)
    ALLOCATE (xm(1:nion))
    DO i = 1, nion
      xm(i) = amu(np(i))*1836D0*ame
    END DO
! xm(:)=amu(np(:))*1836D0*ame
    CALL leapfr(cx(1:nion), cy(1:nion), cz(1:nion), &
                cpx(1:nion), cpy(1:nion), cpz(1:nion), &
                dt1, xm, nion, 4)
    DEALLOCATE (xm)

! correct ionic c.m. to restore value before step
    IF (tfixcmion) THEN
      xcmion = xcmion - SUM(cx(1:nion))/nion
      ycmion = ycmion - SUM(cy(1:nion))/nion
      zcmion = zcmion - SUM(cz(1:nion))/nion
      cx = cx + xcmion
      cy = cy + ycmion
      cz = cz + zcmion
    END IF

  END IF

! PROPAGATION OF MOMENTA

#if(extended)
  IF (ipseudo) THEN
    CALL updatesubgrids
  END IF
#endif

! compute forces
  CALL getforces(rho, psi, it, 0)
! correct forces to have no effect on c.m.
  IF (tfixcmion) THEN
    fxcm = SUM(fx(1:nion))/nion
    fycm = SUM(fy(1:nion))/nion
    fzcm = SUM(fz(1:nion))/nion
    fx(1:nion) = fx(1:nion) - fxcm
    fy(1:nion) = fy(1:nion) - fycm
    fz(1:nion) = fz(1:nion) - fzcm
  END IF
#if(extended)
  IF (tfs < taccel - 1D-6) THEN
    ALLOCATE (tfac(1:nion))
    tfac(:) = ame*amu(np(:))*1836D0*0.048D0/taccel
    IF (tsmooth) THEN
      dtaccel = taccel/10D0
      tfac = tfac*taccel/(taccel - dtaccel/2D0)
      IF (tfs > taccel - dtaccel) tfac = tfac*sin(PI*(taccel - tfs)/(dtaccel*2D0))**2
    END IF
    WRITE (*, *) ' tfs,tfac=', tfs, tfac
    fx(nion) = fx(nion) + vpx*tfac(nion)
    fy(nion) = fy(nion) + vpy*tfac(nion)
    fz(nion) = fz(nion) + vpz*tfac(nion)
    IF (nion > 1) THEN
      fx(1:nion - 1) = fx(1:nion - 1) - vpx*tfac(1:nion - 1)/(nion - 1)
      fy(1:nion - 1) = fy(1:nion - 1) - vpy*tfac(1:nion - 1)/(nion - 1)
      fz(1:nion - 1) = fz(1:nion - 1) - vpz*tfac(1:nion - 1)/(nion - 1)
    END IF
    DEALLOCATE (tfac)
  END IF
#endif

! compute forces on molecule ions with new positions

#if(raregas)
  IF (isurf /= 0) THEN
! propagation of rare gas cores
    IF (nc > 0) THEN
      IF (ALLOCATED(xm)) DEALLOCATE (xm)
      ALLOCATE (xm(1:nc))
      xm = 1D0 ! setting for propagation of momenta
      CALL leapfr(pxc(1:nc), pyc(1:nc), pzc(1:nc), &
                  fxc(1:nc), fyc(1:nc), fzc(1:nc), dt1, xm, nc, 1)
      DEALLOCATE (xm)
    END IF

    IF (nk > 0) THEN
      IF (ALLOCATED(xm)) DEALLOCATE (xm)
      ALLOCATE (xm(1:nk))
      xm = 1D0 ! setting for propagation of momenta
      CALL leapfr(pxk(1:nk), pyk(1:nk), pzk(1:nk), &
                  fxk(1:nk), fyk(1:nk), fzk(1:nk), dt1, xm, nk, 3)
      DEALLOCATE (xm)
    END IF
  END IF
#endif

! propagation of molecule ions
  IF (nion > 0 .AND. ionmdtyp /= 0) THEN
    IF (ALLOCATED(xm)) DEALLOCATE (xm)
    ALLOCATE (xm(1:nion))
    xm = 1D0 ! setting for propagation of momenta
    CALL leapfr(cpx(1:nion), cpy(1:nion), cpz(1:nion), &
                fx(1:nion), fy(1:nion), fz(1:nion), dt1, xm, nion, 4)
    DEALLOCATE (xm)
  END IF

  RETURN
END SUBROUTINE itstep

! ****************************************************

SUBROUTINE leapfr(x, y, z, xprop, yprop, zprop, ddt, xm, n, ityp)

! Computes one leap frog step for 'n' particles of masses 'xm'
! Propagation of positions and momenta are separate
! if x,y,z are positions, then xprop,... are momenta : ddt=dt, xm=mass
! if x,y,z are momenta, then xprop,... are forces : ddt=dt, xm=1.
! ityp = type of ions/atoms (here only active for 'ityp==4')

  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: x(n)
  REAL(DP), INTENT(OUT) :: y(n)
  REAL(DP), INTENT(OUT) :: z(n)
  REAL(DP), INTENT(IN) :: xprop(n)
  REAL(DP), INTENT(IN) :: yprop(n)
  REAL(DP), INTENT(IN) :: zprop(n)
  REAL(DP), INTENT(IN) :: ddt
  REAL(DP), INTENT(IN) :: xm(1:n)
  INTEGER, INTENT(IN) :: ityp

  INTEGER :: i

#if(raregas)
  DO i = 1, n
    IF ((ityp == 4) .OR. (ityp == 1 .AND. imobc(i) == 1) .OR. &
        (ityp == 2 .AND. imobe(i) == 1) .OR. &
        (ityp == 3 .AND. imobk(i) == 1)) THEN
      x(i) = x(i) + ddt*xprop(i)/xm(i)
      y(i) = y(i) + ddt*yprop(i)/xm(i)
      z(i) = z(i) + ddt*zprop(i)/xm(i)

    END IF
  END DO
#else
  IF (ityp == 4) THEN
    DO i = 1, n
      x(i) = x(i) + xprop(i)*(ddt/xm(i))
      y(i) = y(i) + yprop(i)*(ddt/xm(i))
      z(i) = z(i) + zprop(i)*(ddt/xm(i))
    END DO
  END IF
#endif

  RETURN
END SUBROUTINE leapfr

!-----itstepv-------------------------------------------------

SUBROUTINE itstepv(rho, it, psi)

! Ionic time step with velocity Verlet.
! Input:
! rho = local electron density
! psi = set of s.p. wavefunctions
! it = time step number in calling routine
! ionic positions and velocities communicated via module 'params'

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: it
  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)

  REAL(DP), ALLOCATABLE :: xm(:) ! array for raregas particle masses
  REAL(DP), ALLOCATABLE, SAVE :: xmion(:) ! array for ion masses
  REAL(DP), ALLOCATABLE :: tfac(:) ! acceleration factors
  LOGICAL, PARAMETER :: tsmooth = .true. ! switch to smooth termination of acceleration
  LOGICAL, PARAMETER :: ttestprint = .false.

  REAL(DP) :: dtaccel, xcmion, ycmion, zcmion
  REAL(DP) :: fxcm, fycm, fzcm
  INTEGER :: i
#if(raregas)
  REAL(DP) :: pxcmion, pycmion, pzcmion
#endif

#if(extended)
! what time is it ?
  tfs = it*dt1*0.0484D0/(2D0*ame)
#endif
  IF (ttestprint) WRITE (*, *) 'enter ITSTEPV: time[fs],cpx,cpy,cpz=', &
    it*dt1*0.0484D0, cpx(nion), cpy(nion), cpz(nion)

#if(raregas)
! zero in mobile atoms
  IF (isurf /= 0) THEN
    DO i = 1, nc
      IF (imobc(i) == 0) THEN
        pxc(i) = 0D0
        pyc(i) = 0D0
        pzc(i) = 0D0
      END IF
    END DO
    DO i = 1, nk
      IF (imobk(i) == 0) THEN
        pxk(i) = 0D0
        pyk(i) = 0D0
        pzk(i) = 0D0
      END IF
    END DO
  END IF
#endif

! PROPAGATION OF POSITIONS

#if(raregas)
  IF (isurf /= 0 .AND. nc > 0) THEN
    IF (ALLOCATED(xm)) DEALLOCATE (xm)
    ALLOCATE (xm(1:nc))
    xm = mion*1836D0*ame
    CALL velverlet1(xc(1:nc), yc(1:nc), zc(1:nc), pxc(1:nc), pyc(1:nc), pzc(1:nc), &
                    fxc(1:nc), fyc(1:nc), fzc(1:nc), dt1, xm, nc, 1)
    DEALLOCATE (xm)
  END IF
  IF (nk > 0) THEN
! correction of velocities to fix c.m.
    IF (tfixcmion) THEN
      pxcmion = SUM(pxk(1:nk))/nk
      pycmion = SUM(pyk(1:nk))/nk
      pzcmion = SUM(pzk(1:nk))/nk
      pxk = pxk - pxcmion
      pyk = pyk - pycmion
      pzk = pzk - pzcmion
    END IF
    IF (ALLOCATED(xm)) DEALLOCATE (xm)
    ALLOCATE (xm(1:nk))
    xm = mkat*1836D0*ame
    CALL velverlet1(xk(1:nk), yk(1:nk), zk(1:nk), pxk(1:nk), pyk(1:nk), pzk(1:nk), &
                    fxk(1:nk), fyk(1:nk), fzk(1:nk), dt1*modionstep, xm, nk, 3)
    DEALLOCATE (xm)
  END IF
#endif

  IF (nion > 0 .AND. ionmdtyp /= 0) THEN

! compute c.m. of ions
    IF (tfixcmion) THEN
      xcmion = SUM(cx(1:nion))/nion
      ycmion = SUM(cy(1:nion))/nion
      zcmion = SUM(cz(1:nion))/nion
    END IF

! propagation of molecule ions, first half of momenta
    IF (.NOT. ALLOCATED(xmion)) ALLOCATE (xmion(1:nion))
    DO i = 1, nion
      xmion(i) = amu(np(i))*1836D0*ame
    END DO
    CALL velverlet1(cx(1:nion), cy(1:nion), cz(1:nion), &
                    cpx(1:nion), cpy(1:nion), cpz(1:nion), &
                    fx(1:nion), fy(1:nion), fz(1:nion), dt1*modionstep, xmion, nion, 4)

! correct ionic c.m. to restore value before step
    IF (tfixcmion) THEN
      xcmion = xcmion - SUM(cx(1:nion))/nion
      ycmion = ycmion - SUM(cy(1:nion))/nion
      zcmion = zcmion - SUM(cz(1:nion))/nion
      cx = cx + xcmion
      cy = cy + ycmion
      cz = cz + zcmion
    END IF

  END IF

! update forces

#if(extended)
  IF (ipseudo) THEN
    CALL updatesubgrids
  END IF
#endif

! WRITE(*,*) 'before GETFORCES'

! compute forces with new positions
  CALL getforces(rho, psi, it, 0)
! correct forces to have no effect on c.m.
  fxcm = SUM(fx(1:nion))/nion
  fycm = SUM(fy(1:nion))/nion
  fzcm = SUM(fz(1:nion))/nion
  WRITE (*, '(a,3(1pg13.5))') ' c.m. forces=', fxcm, fycm, fzcm
  IF (tfixcmion) THEN
    fx(1:nion) = fx(1:nion) - fxcm
    fy(1:nion) = fy(1:nion) - fycm
    fz(1:nion) = fz(1:nion) - fzcm
  END IF

! WRITE(*,*) 'after GETFORCES'

#if(extended)
  IF (tfs < taccel - 1D-6) THEN
    ALLOCATE (tfac(1:nion))
    tfac(:) = ame*amu(np(:))*1836D0*0.048D0/taccel
    IF (tsmooth) THEN
      dtaccel = taccel/10D0
      tfac = tfac*taccel/(taccel - dtaccel/2D0)
      IF (tfs > taccel - dtaccel) tfac = tfac*sin(PI*(taccel - tfs)/(dtaccel*2D0))**2
    END IF
    WRITE (*, *) ' tfs,tfac=', tfs, tfac
    fx(nion) = vpx*tfac(nion)
    fy(nion) = vpy*tfac(nion)
    fz(nion) = vpz*tfac(nion)
    IF (nion > 1) THEN
      fx(1:nion - 1) = -vpx*tfac(1:nion - 1)/(nion - 1)
      fy(1:nion - 1) = -vpy*tfac(1:nion - 1)/(nion - 1)
      fz(1:nion - 1) = -vpz*tfac(1:nion - 1)/(nion - 1)
    END IF
    DEALLOCATE (tfac)
  END IF
#endif

! second half of propagation of momenta

#if(raregas)
  IF (isurf /= 0) THEN
! propagation of rare gas cores
    IF (nc > 0) THEN
      CALL velverlet2(pxc(1:nc), pyc(1:nc), pzc(1:nc), &
                      fxc(1:nc), fyc(1:nc), fzc(1:nc), dt1*modionstep, nc, 1)
    END IF

    IF (nk > 0) THEN
      CALL velverlet2(pxk(1:nk), pyk(1:nk), pzk(1:nk), &
                      fxk(1:nk), fyk(1:nk), fzk(1:nk), dt1*modionstep, nk, 3)
    END IF
  END IF
#endif

! propagation of molecule ions
  IF (nion > 0 .AND. ionmdtyp /= 0) &
    CALL velverlet2(cpx(1:nion), cpy(1:nion), cpz(1:nion), &
                    fx(1:nion), fy(1:nion), fz(1:nion), dt1*modionstep, nion, 4)

  IF (ttestprint) THEN
    WRITE (6, '(a,3f15.6)') 'fionx,fiony,fionz', fx(1), fy(1), fz(1)
    WRITE (6, '(a,3f15.6)') 'pix,piy,piz', cpx(1), cpy(1), cpz(1)
    WRITE (6, '(a,3f15.6)') 'c.m.: fionx,fiony,fionz', &
      SUM(fx(1:nion)), SUM(fy(1:nion)), SUM(fz(1:nion))
    WRITE (6, '(a,3f15.6)') 'c.m.: pionx,piony,pionz', &
      SUM(cpx(1:nion)), SUM(cpy(1:nion)), SUM(cpz(1:nion))
  END IF

  RETURN
END SUBROUTINE itstepv

SUBROUTINE velverlet1(x, y, z, px, py, pz, fox, foy, foz, ddt, xm, n, ityp)

! ****************************************************

! first part of ionic time step with velocity-Verlet
! x,y,z = arrays of positions
! px,py,pz = arrays of momenta
! fox,foy,foz = arrays of forces, before the step
! ddt = size of time step
! xm = ionic mass
! n = number of ions
! ityp = type of ions/atoms (here ONLY for 'ityp==4')
!
! note that momenta are ONLY stepped by a half step,
! waiting for a second half step after computing the new
! forces

  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: x(n), y(n), z(n)
  REAL(DP), INTENT(IN OUT) :: px(n), py(n), pz(n)
  REAL(DP), INTENT(IN) :: fox(n), foy(n), foz(n)
  REAL(DP), INTENT(IN) :: ddt
  REAL(DP), INTENT(IN) :: xm(n)
  INTEGER, INTENT(IN) :: ityp

  INTEGER :: i
  REAL(DP) :: ddth

  ddth = 0.5D0*ddt
#if(raregas)
  DO i = 1, n

    IF ((ityp == 1 .AND. imobc(i) == 1) .OR. &
        (ityp == 2 .AND. imobe(i) == 1) .OR. &
        (ityp == 3 .AND. imobk(i) == 1) .OR. (ityp == 4)) THEN

      x(i) = x(i) + (px(i) + ddth*fox(i))*ddt/xm(i)
      y(i) = y(i) + (py(i) + ddth*foy(i))*ddt/xm(i)
      z(i) = z(i) + (pz(i) + ddth*foz(i))*ddt/xm(i)
      px(i) = px(i) + ddth*fox(i)
      py(i) = py(i) + ddth*foy(i)
      pz(i) = pz(i) + ddth*foz(i)

    END IF
  END DO
#else
  IF (ityp == 4) THEN
    DO i = 1, n
      x(i) = x(i) + (px(i) + ddth*fox(i))*(ddt/xm(i))
      y(i) = y(i) + (py(i) + ddth*foy(i))*(ddt/xm(i))
      z(i) = z(i) + (pz(i) + ddth*foz(i))*(ddt/xm(i))
      px(i) = px(i) + ddth*fox(i)
      py(i) = py(i) + ddth*foy(i)
      pz(i) = pz(i) + ddth*foz(i)
    END DO
  END IF
#endif

  RETURN
END SUBROUTINE velverlet1

SUBROUTINE velverlet2(px, py, pz, fox, foy, foz, ddt, n, ityp)

! ****************************************************

! second part of ionic time step with velocity-Verlet
! px,py,pz = arrays of momenta
! fox,foy,foz = arrays of forces, after the step
! ddt = size of time step
! xm = ionic mass
! n = number of ions
! ityp = type of ions/atoms

  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: px(n), py(n), pz(n)
  REAL(DP), INTENT(IN) :: fox(n), foy(n), foz(n)
  REAL(DP), INTENT(IN) :: ddt
  INTEGER, INTENT(IN) :: ityp

  INTEGER :: i
  REAL(DP) :: ddth

  ddth = 0.5D0*ddt

#if(raregas)
  DO i = 1, n

    IF ((ityp == 1 .AND. imobc(i) == 1) .OR. &
        (ityp == 2 .AND. imobe(i) == 1) .OR. &
        (ityp == 3 .AND. imobk(i) == 1) .OR. (ityp == 4)) THEN

      px(i) = px(i) + ddth*fox(i)
      py(i) = py(i) + ddth*foy(i)
      pz(i) = pz(i) + ddth*foz(i)

    END IF
  END DO
#else
  IF (ityp == 4) THEN
    DO i = 1, n
      px(i) = px(i) + ddth*fox(i)
      py(i) = py(i) + ddth*foy(i)
      pz(i) = pz(i) + ddth*foz(i)
    END DO
  END IF
#endif

  RETURN
END SUBROUTINE velverlet2

! *******************************

SUBROUTINE lffirststep(rho, psi)

! Preparation of leap-frog steps by a few short Euler steps.
! Input:
! rho = local electron density
! psi = set of s.p. wavefunctions

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)

  REAL(DP), ALLOCATABLE :: xm(:)

  REAL(DP), EXTERNAL :: enerkin_ions

  tfs = 0D0

#if(extended)
  IF (ipseudo) THEN
    CALL updatesubgrids
  END IF
#endif

  CALL getforces(rho, psi, -1, 0)

#if(raregas)
  IF (isurf /= 0) THEN
    IF (nc + NE > 0) THEN
! propagation of cores
      IF (ALLOCATED(xm)) DEALLOCATE (xm)
      ALLOCATE (xm(1:nc))
      xm = 1D0
      CALL leapfr(pxc(1:nc), pyc(1:nc), pzc(1:nc), &
                  fxc(1:nc), fyc(1:nc), fzc(1:nc), dt1/2D0, xm, nc, 1)
      DEALLOCATE (xm)
! propagation of clouds
      IF (ifadiadip /= 1) THEN
        ALLOCATE (xm(1:ne)) ! possible buf with gfortran = debian64
        xm = 1D0
        CALL leapfr(pxe(1:ne), pye(1:ne), pze(1:ne), &
                    fxe(1:ne), fye(1:ne), fze(1:ne), dt1/4D0, xm, ne, 2)
        DEALLOCATE (xm)
      END IF
    END IF
! propagation of cation
    IF (nk > 0) THEN
      IF (ALLOCATED(xm)) DEALLOCATE (xm)
      ALLOCATE (xm(1:nk))
      xm = 1D0
      CALL leapfr(pxk(1:nk), pyk(1:nk), pzk(1:nk), &
                  fxk(1:nk), fyk(1:nk), fzk(1:nk), dt1/2D0, xm, nk, 3)
      DEALLOCATE (xm)
    END IF

  END IF
#endif

! propagation of molecule ions
  IF (nion > 0) THEN
    IF (ALLOCATED(xm)) DEALLOCATE (xm)
    ALLOCATE (xm(1:nion))
    xm = 1D0
    CALL leapfr(cpx(1:nion), cpy(1:nion), cpz(1:nion), &
                fx(1:nion), fy(1:nion), fz(1:nion), dt1/2D0, xm, nion, 4)
    DEALLOCATE (xm)
  END IF

  WRITE (6, *) 'initial momenta and kinetic energy'
  ekion = enerkin_ions()

  RETURN
END SUBROUTINE lffirststep

! ******************************

SUBROUTINE spheric(r, x, y, z, n)

! Samples stichastically a spherical distribution of n particles on
! a sphere of radius r.
!
! Input:
! r = radius of sampling sphere
! n = number of particles
! Output:
! x,y,z = arrays of ionic coordinates

  USE params, ONLY: DP, Pi
  IMPLICIT NONE

! ******************************

!INTEGER, PARAMETER :: mm=50000
  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(OUT) :: x(n)
  REAL(DP), INTENT(OUT) :: y(n)
  REAL(DP), INTENT(OUT) :: z(n)

  INTEGER :: i, ia, ic, im, jran
  REAL(DP) :: cth, sth, cosp, sinp, phi, rr, xx

  DATA im, ia, ic/259200, 7141, 54773/
  jran = 12345

  WRITE (6, *) r, n
  DO i = 1, n
    rr = r
    jran = MOD(jran*ia + ic, im)
    xx = REAL(jran, DP)/REAL(im, DP)
    cth = 1D0 - 2D0*xx
    sth = SQRT(1D0 - cth*cth)
    jran = MOD(jran*ia + ic, im)
    xx = REAL(jran, DP)/REAL(im, DP)
    phi = 2D0*Pi*xx
    cosp = COS(phi)
    sinp = SIN(phi)
    x(i) = rr*sth*cosp
    y(i) = rr*sth*sinp
    z(i) = rr*cth
  END DO
  RETURN
END SUBROUTINE spheric

! ************************************

SUBROUTINE conslw(x, y, z, px, py, pz, n)
  USE params, ONLY: DP
  IMPLICIT NONE

! Renormalization of c.o.m, of momentum and of angular momentum
!
! Input:
! n = number of particles
! Input/Output:
! x,y,z = arrays of ionic coordinates
! px,py,pz = arrays of ionic momenta

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(IN OUT) :: x(n)
  REAL(DP), INTENT(IN OUT) :: y(n)
  REAL(DP), INTENT(IN OUT) :: z(n)
  REAL(DP), INTENT(IN OUT) :: px(n)
  REAL(DP), INTENT(IN OUT) :: py(n)
  REAL(DP), INTENT(IN OUT) :: pz(n)

  INTEGER :: i, npart
  REAL(DP) :: rlx, rly, rlz, reno
  REAL(DP) :: rixx, riyy, rizz, rixy, rixz, riyz
  REAL(DP) :: rkx, rky, rkz, rx, ry, rz, rx2, ry2, rz2
  REAL(DP) :: rlxcm, rlycm, rlzcm, xcm, ycm, zcm, xkcm, ykcm, zkcm, &
              xxkcm, yykcm, zzkcm, xxcm, yycm, zzcm
  REAL(DP) :: det, wx, wy, wz

  xcm = 0D0
  ycm = 0D0
  zcm = 0D0
  xkcm = 0D0
  ykcm = 0D0
  zkcm = 0D0
  rlx = 0D0
  rly = 0D0
  rlz = 0D0
  rixx = 0D0
  riyy = 0D0
  rizz = 0D0
  rixy = 0D0
  rixz = 0D0
  riyz = 0D0
  npart = n
  reno = 1D0/REAL(n, DP)

  DO i = 1, n
    rx = x(i)
    ry = y(i)
    rz = z(i)
    rkx = px(i)
    rky = py(i)
    rkz = pz(i)
    rx2 = rx*rx
    ry2 = ry*ry
    rz2 = rz*rz

! c.o.m. in r space

    xcm = xcm + rx
    ycm = ycm + ry
    zcm = zcm + rz

! c.o.m. in p space

    xkcm = xkcm + rkx
    ykcm = ykcm + rky
    zkcm = zkcm + rkz

! angular momentum

    rlx = rlx + ry*rkz - rz*rky
    rly = rly + rz*rkx - rx*rkz
    rlz = rlz + rx*rky - ry*rkx

! tensor of inertia

    rixx = rixx + ry2 + rz2
    riyy = riyy + rx2 + rz2
    rizz = rizz + rx2 + ry2
    rixy = rixy - rx*ry
    rixz = rixz - rx*rz
    riyz = riyz - ry*rz
  END DO

! renormalization of c.o.m. in r and p space
! and of the total angular momentum

  xcm = xcm*reno
  ycm = ycm*reno
  zcm = zcm*reno
  xkcm = xkcm*reno
  ykcm = ykcm*reno
  zkcm = zkcm*reno

  rixx = rixx*reno - ycm*ycm - zcm*zcm
  riyy = riyy*reno - xcm*xcm - zcm*zcm
  rizz = rizz*reno - xcm*xcm - ycm*ycm
  rixy = rixy*reno + xcm*ycm
  rixz = rixz*reno + xcm*zcm
  riyz = riyz*reno + ycm*zcm

  rlxcm = ycm*zkcm - zcm*ykcm
  rlycm = zcm*xkcm - xcm*zkcm
  rlzcm = xcm*ykcm - ycm*xkcm

  rlx = rlx*reno - rlxcm
  rly = rly*reno - rlycm
  rlz = rlz*reno - rlzcm

  det = rixx*riyy*rizz + 2.*rixy*rixz*riyz &
        - rixx*riyz**2 - riyy*rixz**2 - rizz*rixy**2
  wx = (rlx*(riyy*rizz - riyz**2) + rly*(rixz*riyz - rizz*rixy) &
        + rlz*(rixy*riyz - riyy*rixz))/det
  wy = (rlx*(rixz*riyz - rixy*rizz) + rly*(rixx*rizz - rixz**2) &
        + rlz*(rixy*rixz - rixx*riyz))/det
  wz = (rlx*(rixy*riyz - riyy*rixz) &
        + rly*(rixz*rixy - rixx*riyz) + rlz*(rixx*riyy - rixy**2))/det

  WRITE (6, *) ' before...'
  WRITE (6, *) ' r c.m....', xcm, ycm, zcm
  WRITE (6, *) ' k c.m....', xkcm, ykcm, zkcm
  WRITE (6, *) ' l........', rlx, rly, rlz

  xxcm = 0D0
  yycm = 0D0
  zzcm = 0D0
  xxkcm = 0D0
  yykcm = 0D0
  zzkcm = 0D0
  rlx = 0D0
  rly = 0D0
  rlz = 0D0

  DO i = 1, n
    x(i) = x(i) - xcm
    y(i) = y(i) - ycm
    z(i) = z(i) - zcm
    xxcm = xxcm + x(i)
    yycm = yycm + y(i)
    zzcm = zzcm + z(i)
    px(i) = px(i) - xkcm - (wy*z(i) - wz*y(i))
    py(i) = py(i) - ykcm - (wz*x(i) - wx*z(i))
    pz(i) = pz(i) - zkcm - (wx*y(i) - wy*x(i))
    xxkcm = xxkcm + px(i)
    yykcm = yykcm + py(i)
    zzkcm = zzkcm + pz(i)
    rx = x(i)
    ry = y(i)
    rz = z(i)
    rkx = px(i)
    rky = py(i)
    rkz = pz(i)
    rlx = rlx + ry*rkz - rz*rky
    rly = rly + rz*rkx - rx*rkz
    rlz = rlz + rx*rky - ry*rkx
    rixx = rixx + ry*ry + rz*rz
    rixy = rixy + rx*ry
    rixz = rixz + rx*rz
    riyy = riyy + rx*rx + rz*rz
    riyz = riyz + ry*rz
    rizz = rizz + rx*rx + ry*ry
  END DO

  xxcm = xxcm*reno
  yycm = yycm*reno
  zzcm = zzcm*reno
  xxkcm = xxkcm*reno
  yykcm = yykcm*reno
  zzkcm = zzkcm*reno
  rlx = rlx*reno
  rly = rly*reno
  rlz = rlz*reno

  WRITE (6, *) ' after....'
  WRITE (6, *) ' r c.m....', xxcm, yycm, zzcm
  WRITE (6, *) ' k c.m....', xxkcm, yykcm, zzkcm
  WRITE (6, *) ' l .......', rlx, rly, rlz

! END of the renormalizations ( r , p , l )

  RETURN
END SUBROUTINE conslw

