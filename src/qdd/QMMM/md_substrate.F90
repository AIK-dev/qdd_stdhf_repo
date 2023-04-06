
#if(raregas)
!-----valence_step---------------------------------------------

SUBROUTINE valence_step(rho, dt, it, tdyn)

! Propagates the substrates valence electron cloud.
! Input:
! rho = electron density
! dt = stepsize (IN CASE of dynamics)
! it = current iteration
! tdyn = switch TO dynamic CASE
! Output:
! valence positions via COMMON

! This routine can still be used IN statics as well as
! IN dynamics (see the switch 'tdyn'). It is presently
! ONLY called IN dynamics.

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  REAL(DP), INTENT(IN) :: dt
  INTEGER, INTENT(IN) :: it
  LOGICAL, INTENT(IN) :: tdyn

#if(raregas)
  COMPLEX(DP) :: psidummy(1)
#endif
  LOGICAL, PARAMETER :: tspinprint = .true.

!---------------------------------------------------------------

  IF (isurf /= 0 .AND. nc + NE + nk > 0) THEN ! check condition ??

    IF (ifadiadip == 1) THEN
! instantaneous adjustment of substrate dipoles
      CALL adjustdip(rho, it)
    ELSE
! dynamic propagation of substrate dipoles
      IF (tdyn) THEN
        IF (ipsptyp == 1) STOP ' VSTEP must not be used with Goedecker PsP' ! ???
        IF (ionmdtyp == 1) CALL vstep(rho, psidummy, it, dt)
        IF (ionmdtyp == 2) CALL vstepv(rho, psidummy, it, dt)
      ELSE
        CALL adjustdip(rho, it)
      END IF
    END IF

! update of the subgrids ? here or outside substrate loop ?

    IF (.NOT. tdyn .AND. ipseudo) CALL updatesubgrids ! probably done IN vstep
    CALL calcpseudo() !! ????????

    IF (nc > 0 .AND. ivdw == 1) CALL calc_frho(rho)

  END IF

  RETURN
END SUBROUTINE valence_step
#endif

!------------------------------------------------------------

REAL(DP) FUNCTION sigsig(s1, s2)
!------------------------------------------------------------
  USE params, ONLY: DP
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: s1, s2

  sigsig = SQRT(s1**2 + s2**2)

  RETURN
END FUNCTION sigsig
!------------------------------------------------------------

#if(raregas)
!------------------------------------------------------------

SUBROUTINE initsurface
!------------------------------------------------------------
  USE params
  USE util, ONLY: printfield, printfield2
  IMPLICIT NONE

  INTEGER :: i, icheckflag, ico1, ico2, ii, ind, ion, ix, iy, iz
  INTEGER :: icountc, icounte, icountk, iefixed, icfixed, ikfixed
  INTEGER :: jj, jx, jy, jz, ntrans1, ntrans2
  REAL(DP) :: ala, be, dum, cf, cqq, rc, rr, rarmass, sigsig, z0, z1, zmin
  REAL(DP) :: xcm, ycm, zcm
  REAL(DP) :: hfield(kdfull2), hfield2(kdfull2)

  INTEGER, EXTERNAL :: isoutofbox
#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
  CALL mpi_comm_rank(mpi_comm_world, myn, mpi_ierror)
#else
  myn = 0
#endif

  NE = nc

  IF (myn == 0) THEN
    OPEN (120, STATUS='old', FILE='for005surf.'//outnam)
    READ (120, *) ! comment line

    IF (iusecell == 0) THEN

      DO i = 1, nc
        READ (120, *) xc(i), yc(i), zc(i), xe(i), ye(i), ze(i), imobc(i), imobe(i)
      END DO

      DO i = 1, nk
        READ (120, *) xk(i), yk(i), zk(i), imobk(i)
      END DO

    ELSE IF (iusecell == 1) THEN ! build from UNIT cell

      READ (120, *) icheckflag

      IF (icheckflag /= iusecell) &
        STOP 'for005surf incompatible for this iUseCell'

      READ (120, *) rlattvec(1), rlattvec(2), rlattvec(3)
      READ (120, *) ncells1x, ncells1y, ncells1z
      READ (120, *) ncells2x, ncells2y, ncells2z
      READ (120, *) natoms

      DO i = 1, natoms
        READ (120, *) xatom(i, 1), xatom(i, 2), xatom(i, 3), ipsort(i)
        WRITE (6, '(3f20.8,i6)') xatom(i, 1), xatom(i, 2), xatom(i, 3), ipsort(i)

      END DO

! build region 1

      icountc = 0
      icounte = 0
      icountk = 0

      WRITE (6, *) 'a', zc(1), ze(1)

      DO jx = -ncells1x, ncells1x
        DO jy = -ncells1y, ncells1y
          DO jz = 1, ncells1z

            DO i = 1, natoms

              rvectmp(1) = xatom(i, 1) + rlattvec(1)*jx + shiftx
              rvectmp(2) = xatom(i, 2) + rlattvec(2)*jy + shifty
              rvectmp(3) = xatom(i, 3) - rlattvec(3)*(jz - 1) + shiftz

! WRITE(6,*) rVecTmp(1), rVecTmp(2), rVecTmp(3)

              IF (ipsort(i) == 1) THEN
                icountc = icountc + 1
                xc(icountc) = rvectmp(1)
                yc(icountc) = rvectmp(2)
                zc(icountc) = rvectmp(3)
                chgc(icountc) = chgc0
                imobc(icountc) = 0
              ELSE IF (ipsort(i) == 2) THEN
                icounte = icounte + 1
                xe(icounte) = rvectmp(1)
                ye(icounte) = rvectmp(2)
                ze(icounte) = rvectmp(3)
                chge(icounte) = chge0
                imobe(icounte) = 0
              ELSE IF (ipsort(i) == 3) THEN
                icountk = icountk + 1
                xk(icountk) = rvectmp(1)
                yk(icountk) = rvectmp(2)
                zk(icountk) = rvectmp(3)
                chgk(icountk) = chgk0
                imobk(icountk) = 0
              ELSE
                STOP 'md.F: invalid particle sort'
              END IF

            END DO

          END DO
        END DO
      END DO

      nc = icountc
      NE = icounte
      nk = icountk

      WRITE (6, *) 'nc, nk = ', nc, nk

      IF (nc /= NE) STOP 'md.F: nc must equal ne!!'

! calculate cell multipole moments
      WRITE (6, *) 'calculate cell multipole moments'
      DO ico1 = 1, 3
        celldipole(ico1) = 0D0
        DO ico2 = 1, 3
          cellquad(ico1, ico2) = 0D0
        END DO
      END DO

!$$$ DO ilm=-2,2
!$$$ cellMult(2,ilm) = CMPLX(0D0,0D0,DP)
!$$$ END DO

      DO i = 1, natoms
        IF (ipsort(i) == 1) THEN
          cqq = chgc0
        ELSE IF (ipsort(i) == 2) THEN
          cqq = chge0
        ELSE
          cqq = chgk0
        END IF

        DO ico1 = 1, 3
          celldipole(ico1) = celldipole(ico1) + xatom(i, ico1)*cqq
          DO ico2 = 1, 3
            IF (ico1 == ico2) THEN
              rr = xatom(i, 1)**2 + xatom(i, 2)**2 + xatom(i, 3)**2
              cellquad(ico1, ico2) = cellquad(ico1, ico2) + &
                                     cqq*(3*xatom(i, ico1)*xatom(i, ico2) - rr)
            ELSE
              cellquad(ico1, ico2) = cellquad(ico1, ico2) + &
                                     cqq*3*xatom(i, ico1)*xatom(i, ico2)
            END IF
          END DO
        END DO

      END DO

      WRITE (6, '(a,3f17.7)') 'Cell dipole moment: ', celldipole(1), &
        celldipole(2), celldipole(3)
      WRITE (6, *) 'Cell quadrupole moment: '
      WRITE (6, '(3f17.7)') cellquad(1, 1), cellquad(1, 2), cellquad(1, 3)
      WRITE (6, '(3f17.7)') cellquad(2, 1), cellquad(2, 2), cellquad(2, 3)
      WRITE (6, '(3f17.7)') cellquad(3, 1), cellquad(3, 2), cellquad(3, 3)

!$$$ CALL createOuterPot
!get electrostatic potential
! from outer region by multipole expansion

    ELSE IF (iusecell == 2) THEN ! construct lattice by translational vectors

      READ (120, *) icheckflag

      IF (icheckflag /= iusecell) &
        STOP 'for005surf incompatible for this iUseCell'

      READ (120, *) rtransvec(1, 1), rtransvec(1, 2)
      READ (120, *) rtransvec(2, 1), rtransvec(2, 2)
      READ (120, *) ntrans1, ntrans2
      READ (120, *) ncells2x, ncells2y, ncells2z
      READ (120, *) natoms

      DO ii = 1, natoms
        READ (120, *) xatom(ii, 1), xatom(ii, 2), xatom(ii, 3), ipsort(ii)
        xatom(ii, 3) = xatom(ii, 3) + shiftz
      END DO

      icountc = 0
      icountk = 0
      icounte = 0

      DO ii = -ntrans1, ntrans1
        DO jj = -ntrans2, ntrans2

          DO ind = 1, natoms
            rvectmp(1) = xatom(ind, 1) + ii*rtransvec(1, 1) + jj*rtransvec(2, 1)
            rvectmp(2) = xatom(ind, 2) + ii*rtransvec(1, 2) + jj*rtransvec(2, 2)
            rvectmp(3) = xatom(ind, 3) + ii*rtransvec(1, 3) + jj*rtransvec(2, 3)

            IF (ipsort(ind) == 1) THEN
              icountc = icountc + 1
              xc(icountc) = rvectmp(1)
              yc(icountc) = rvectmp(2)
              zc(icountc) = rvectmp(3)
              chgc(icountc) = chgc0
              imobc(icountc) = 0
            ELSE IF (ipsort(ind) == 2) THEN
              icounte = icounte + 1
              xe(icounte) = rvectmp(1)
              ye(icounte) = rvectmp(2)
              ze(icounte) = rvectmp(3)
              chge(icounte) = chge0
              imobe(icounte) = 0
            ELSE IF (ipsort(ind) == 3) THEN
              icountk = icountk + 1
              xk(icountk) = rvectmp(1)
              yk(icountk) = rvectmp(2)
              zk(icountk) = rvectmp(3)
              chgk(icountk) = chgk0
              imobk(icountk) = 0
            ELSE
              STOP 'md.F: invalid particle sort'
            END IF

          END DO

        END DO
      END DO

      nc = icountc
      nk = icountk
      NE = nc

      WRITE (6, '(a,2i6)') 'nc, nk: ', nc, nk

    ELSE
      STOP 'Error IN md.F: iUseCell'
    END IF

    CLOSE (120)

  END IF ! myn?

! CALL debugMultExp(natoms,ipsort)
! CALL debugFME(natoms,ipsort)

! reshuffle initial DATA TO array

  DO i = 1, nc
    chgc(i) = chgc0
    chge(i) = chge0
  END DO
  DO i = 1, nk
    chgk(i) = chgk0
  END DO

  ndfull2 = kdfull2

! shift inner region of surface IN space IF desired

  IF (shiftx**2 + shifty**2 + shiftz**2 > small .AND. iusecell == 0) THEN
    DO i = 1, nc
      xc(i) = xc(i) + shiftx
      yc(i) = yc(i) + shifty
      zc(i) = zc(i) + shiftz
      xe(i) = xe(i) + shiftx
      ye(i) = ye(i) + shifty
      ze(i) = ze(i) + shiftz
    END DO
    DO i = 1, nk
      xk(i) = xk(i) + shiftx
      yk(i) = yk(i) + shifty
      zk(i) = zk(i) + shiftz
    END DO
  END IF

! shift done.

! re-scale distances IN inner region IF desired

  IF (ABS(scaledist - 1D0) > small .AND. iusecell == 0) THEN

    xcm = 0D0
    ycm = 0D0
    zcm = 0D0

    DO i = 1, nc
      xcm = xcm + xc(i)
      ycm = ycm + yc(i)
      zcm = zcm + zc(i)
      xcm = xcm + xe(i)
      ycm = ycm + ye(i)
      zcm = zcm + ze(i)
    END DO
    DO i = 1, nk
      xcm = xcm + xk(i)
      ycm = ycm + yk(i)
      zcm = zcm + zk(i)
    END DO

    xcm = xcm/(nc + nk + NE)
    ycm = ycm/(nc + nk + NE)
    zcm = zcm/(nc + nk + NE)

    DO i = 1, nc
      xc(i) = (xc(i) - xcm)*scaledist + xcm
      yc(i) = (yc(i) - ycm)*scaledist + ycm
      zc(i) = (zc(i) - zcm)*scaledist + zcm
      xe(i) = (xe(i) - xcm)*scaledist + xcm
      ye(i) = (ye(i) - ycm)*scaledist + ycm
      ze(i) = (ze(i) - zcm)*scaledist + zcm
    END DO
    DO i = 1, nk
      xk(i) = (xk(i) - xcm)*scaledist + xcm
      yk(i) = (yk(i) - ycm)*scaledist + ycm
      zk(i) = (zk(i) - zcm)*scaledist + zcm
    END DO

  END IF

! re-scaling done.

! create subgrids for fast integration and evaluation

  WRITE (6, *) 'Init subgrids'
  CALL initsubgrids

  told = 0D0
  ekinel = 0D0
  ekinion = 0D0
  epots = 0D0
  energyg = 0D0

  CALL initimpulses(0)

  IF (ekmat > 0D0) THEN
    rarmass = (REAL(nc, DP))*mion*1836.D0*ame
!c ekmat IN eV
    ekmat = ekmat/13.6D0
!c ekmat IN Ry
    DO ion = 1, nc
      pxc(ion) = (1836.D0*mion*ame)*SQRT(2D0*ekmat/rarmass)
      pxe(ion) = (1836.D0*me*ame)*SQRT(2D0*ekmat/rarmass)
    END DO
    ekion = ekmat
  END IF

  WRITE (6, '(a,1f12.5,2e15.5)') 'ekmat, mta core & cld :', ekmat, pxc(1), pxe(1)

  DO i = 1, nc
    xcold(i) = xc(i)
    ycold(i) = yc(i)
    zcold(i) = zc(i)
    xeold(i) = xe(i)
    yeold(i) = ye(i)
    zeold(i) = ze(i)

    xcinit(i) = xc(i)
    ycinit(i) = yc(i)
    zcinit(i) = zc(i)
    xeinit(i) = xe(i)
    yeinit(i) = ye(i)
    zeinit(i) = ze(i)

  END DO

  DO i = 1, nk
    xkold(i) = xk(i)
    ykold(i) = yk(i)
    zkold(i) = zk(i)

    xkinit(i) = xk(i)
    ykinit(i) = yk(i)
    zkinit(i) = zk(i)

  END DO

!29 IF (iuselast == 1) THEN ! USE ion positions from last run
  IF (iuselast == 1) THEN ! USE ion positions from last run
    OPEN (308, STATUS='old', FILE='for005surf.init')
    READ (308, *) nc, nk
    NE = nc
    DO ii = 1, nc
      READ (308, *) xc(ii), yc(ii), zc(ii), xe(ii), ye(ii), ze(ii), imobc(ii), imobe(ii)
    END DO
    DO ii = 1, nk
      READ (308, *) xk(ii), yk(ii), zk(ii), imobk(ii)
    END DO

    CLOSE (308)
  END IF

! set mixed Gaussian widths:
  sigmacc = sigsig(sigmac, sigmac)
  sigmacv = sigsig(sigmac, sigmav)
  sigmack = sigsig(sigmac, sigmak)
  sigmavv = sigsig(sigmav, sigmav)
  sigmakv = sigsig(sigmak, sigmav)
  sigmakk = sigsig(sigmak, sigmak)

!$$$ IF (iPotFixed.eq.0) THEN
!$$$
!$$$ CALL setzero(potion,kdfull2)
!$$$ CALL setzero(potfixedion,kdfull2)
!$$$ CALL setzero(potStat,kdfull2)
!$$$
!$$$
!$$$
!$$$ IF (nelect.gt.0) CALL addGSMPot(potFixedIon,0) ! add fixed particles
!$$$ IF (iPotFixed.eq.0) CALL addGSMPot(potFixedIon,-1) ! add pure pseudo-particles
!$$$ ! remember that potFixedIon as potion ONLY
!$$$ ! acts on the DFT electrons
!$$$ IF (iPotFixed.eq.0) CALL addGSMPot(potStat,-1)
!$$$
!$$$
!$$$ END IF

!$$$ IF (iPriFixed.ne.0) THEN
!$$$
!$$$ OPEN(543,STATUS='unknown',FILE='potFixed.'//outnam)
!$$$ CALL printField(543,potFixedIon)
!$$$ CLOSE(543)
!$$$c OPEN(544,STATUS='unknown',FILE='potMob.'//outnam)
!$$$c CALL printField(544,potion)
!$$$c CLOSE(544)
!$$$
!$$$
!$$$ IF (iPriFixed.eq.2) STOP
!$$$ & 'STOP after calculating fixed potential'
!$$$
!$$$ END IF

  WRITE (6, *) 'Done.'

! check which remaining ions are OUT of box TO exclude them from
! pseudodensity description

!$$$ DO ii=1,nc
!$$$ CALL getParas(ii)
!$$$ WRITE(6,'(3i,3f)') ii,imobc(ii),iOutOfBox(ii),
!$$$ & rVecTmp(1),rVecTmp(2),rVecTmp(3)
!$$$ END DO
!$$$
!$$$ DO ii=1,ne
!$$$ CALL getParas(ii+nc)
!$$$ WRITE(6,'(3i,3f)') ii,imobe(ii),iOutOfBox(ii+nc),
!$$$ & rVecTmp(1),rVecTmp(2),rVecTmp(3)
!$$$ END DO
!$$$
!$$$ DO ii=1,nk
!$$$ CALL getParas(ii+nc+ne)
!$$$ WRITE(6,'(3i,3f)') ii,imobk(ii),iOutOfBox(ii+nc+ne),
!$$$ & rVecTmp(1),rVecTmp(2),rVecTmp(3)
!$$$ END DO
!$$$
!$$$ STOP

! FIXING AND UNFIXING STUFF...

  IF (iunfixall /= 0) THEN
    DO i = 1, nc
      imobc(i) = 1
      imobe(i) = 1
    END DO
    DO i = 1, nk
      imobk(i) = 1
    END DO
  END IF
  IF (ifixall /= 0) THEN
    DO i = 1, nc
      imobc(i) = 0
      imobe(i) = 0
    END DO
    DO i = 1, nk
      imobk(i) = 0
    END DO
  END IF

  IF (iunfixall /= 0 .AND. ifixall /= 0) &
    STOP 'iUnFixAll != 0 and iFixAll !=0 incompatible'

  IF (unfixcrad > 0) THEN
    DO i = 1, nc
      IF (xc(i)**2 + yc(i)**2 + zc(i)**2 < unfixcrad**2) THEN
        imobc(i) = 1
      END IF
    END DO
  END IF

  IF (unfixerad > 0) THEN
    DO i = 1, NE
      IF (xe(i)**2 + ye(i)**2 + ze(i)**2 < unfixerad**2) THEN
        imobe(i) = 1
      END IF
    END DO
  END IF

  IF (unfixkrad > 0) THEN
    DO i = 1, nk
      IF (xk(i)**2 + yk(i)**2 + zk(i)**2 < unfixkrad**2) THEN
        imobk(i) = 1
      END IF
    END DO
  END IF

  IF (unfixclateralrad > 0) THEN
    DO i = 1, nc
      IF (xc(i)**2 + yc(i)**2 < unfixclateralrad**2) THEN
        imobc(i) = 1
      END IF
    END DO
  END IF

  IF (unfixelateralrad > 0) THEN
    DO i = 1, NE
      IF (xe(i)**2 + ye(i)**2 < unfixelateralrad**2) THEN
        imobe(i) = 1
      END IF
    END DO
  END IF

  IF (unfixklateralrad > 0) THEN
    DO i = 1, nk
      IF (xk(i)**2 + yk(i)**2 < unfixklateralrad**2) THEN
        imobk(i) = 1
      END IF
    END DO
  END IF
! y-direction

  IF (unfixclateralrady > 0) THEN
    DO i = 1, nc
      IF (xc(i)**2 + zc(i)**2 < unfixclateralrady**2) THEN
        imobc(i) = 1
      END IF
    END DO
  END IF

  IF (unfixelateralrady > 0) THEN
    DO i = 1, NE
      IF (xe(i)**2 + ze(i)**2 < unfixelateralrady**2) THEN
        imobe(i) = 1
      END IF
    END DO
  END IF

  IF (unfixklateralrady > 0) THEN
    DO i = 1, nk
      IF (xk(i)**2 + zk(i)**2 < unfixklateralrady**2) THEN
        imobk(i) = 1
      END IF
    END DO
  END IF

! x-direction

  IF (unfixclateralradx > 0) THEN
    DO i = 1, nc
      IF (zc(i)**2 + yc(i)**2 < unfixclateralradx**2) THEN
        imobc(i) = 1
      END IF
    END DO
  END IF

  IF (unfixelateralradx > 0) THEN
    DO i = 1, NE
      IF (ze(i)**2 + ye(i)**2 < unfixelateralradx**2) THEN
        imobe(i) = 1
      END IF
    END DO
  END IF

  IF (unfixklateralradx > 0) THEN
    DO i = 1, nk
      IF (zk(i)**2 + yk(i)**2 < unfixklateralradx**2) THEN
        imobk(i) = 1
      END IF
    END DO
  END IF

! z-direction

  DO i = 1, nc
    IF (zc(i) < fixcbelow) THEN
      imobc(i) = 0
    END IF
  END DO
  DO i = 1, NE
    IF (ze(i) < fixebelow) THEN
      imobe(i) = 0
    END IF
  END DO
  DO i = 1, nk
    IF (zk(i) < fixkbelow) THEN
      imobk(i) = 0
    END IF
  END DO

! y-direction

  DO i = 1, nc
    IF (yc(i) < fixcbelowy) THEN
      imobc(i) = 0
    END IF
  END DO
  DO i = 1, NE
    IF (ye(i) < fixebelowy) THEN
      imobe(i) = 0
    END IF
  END DO
  DO i = 1, nk
    IF (yk(i) < fixkbelowy) THEN
      imobk(i) = 0
    END IF
  END DO

! x-direction

  DO i = 1, nc
    IF (xc(i) < fixcbelowx) THEN
      imobc(i) = 0
    END IF
  END DO
  DO i = 1, NE
    IF (xe(i) < fixebelowx) THEN
      imobe(i) = 0
    END IF
  END DO
  DO i = 1, nk
    IF (xk(i) < fixkbelowx) THEN
      imobk(i) = 0
    END IF
  END DO

  icfixed = 0
  iefixed = 0
  ikfixed = 0

  DO i = 1, nc
    IF (imobc(i) == 0) icfixed = icfixed + 1
    IF (imobe(i) == 0) iefixed = iefixed + 1
  END DO
  DO i = 1, nk
    IF (imobk(i) == 0) ikfixed = ikfixed + 1
  END DO

  WRITE (6, '(i6,a,i6,a)') icfixed, ' of ', nc, ' cores fixed.'
  WRITE (6, '(i6,a,i6,a)') iefixed, ' of ', NE, ' shells fixed.'
  WRITE (6, '(i6,a,i6,a)') ikfixed, ' of ', nk, ' cations fixed.'

  IF (ishiftcmtoorigin) THEN

    WRITE (6, '(a,3f12.3)') 'Shifting surface by : ', rvectmp2(1), &
      rvectmp2(2), rvectmp2(3)

    DO ii = 1, nc
      xc(ii) = xc(ii) - rvectmp2(1)
      yc(ii) = yc(ii) - rvectmp2(2)
      zc(ii) = zc(ii) - rvectmp2(3)
    END DO
    DO ii = 1, NE
      xe(ii) = xe(ii) - rvectmp2(1)
      ye(ii) = ye(ii) - rvectmp2(2)
      ze(ii) = ze(ii) - rvectmp2(3)
    END DO
    DO ii = 1, nk
      xk(ii) = xk(ii) - rvectmp2(1)
      yk(ii) = yk(ii) - rvectmp2(2)
      zk(ii) = zk(ii) - rvectmp2(3)
    END DO

  END IF

  DO i = 1, nc + NE + nk

    CALL getparas(i)
    ioutofbox(i) = isoutofbox(rvectmp(1), rvectmp(2), rvectmp(3))
!$$$ WRITE(6,'(3i,3f12.2)') i,iOutOfBox(i),imobTmp,
!$$$ & rVecTmp(1),rVecTmp(2),rVecTmp(3)
  END DO

! calculate potential of the fixed ions

! iPotFixed allows TO READ/WRITE the electrostatic potential caused by particles
! with imob=0, so that the run-time calculation can be skipped
!
! 0 --> DO not READ IN; calculate electrostatic potential of all (even fixed)
! GSM particles at each iteration
! 1 --> READ IN potFixedIon() from FILE
! -1 --> calculate potFixedIon() once and WRITE the RESULT TO a FILE which can
! be later READ IN by option 1, STOP after that
! 2 --> calculate potFixedIon() once at the beginning and make your run; does
! not WRITE potFixedIon() TO FILE

  IF (ipotfixed == 1) THEN
    WRITE (6, *) 'Reading electrostatic potential of fixed ions from FILE.'
    OPEN (121, STATUS='old', FILE='potFixed.'//outnam)
    READ (121, *)
    DO ind = 1, kdfull2
      READ (121, *) dum, dum, dum, potfixedion(ind)
    END DO
    CLOSE (121)
  ELSE IF (ipotfixed == 2) THEN
    DO i = 1, kdfull2
      potfixedion(i) = 0D0
      hfield(i) = 0D0
    END DO
    WRITE (6, *) 'Calculating potential of fixed ions. '
    WRITE (6, *) 'This may take a WHILE.'
    CALL addgsmpot(potfixedion, 0)
    IF (idielec /= 0) THEN
      CALL addgsmpot2(hfield, 0)
      DO ind = 1, kdfull2
        CALL conv1to3(ind)
        IF (iindtmp(1) <= (nint(xdielec/dx) + nx)) potfixedion(ind) = hfield(ind)
      END DO
    END IF
    CALL addshortrepulsivepot(potfixedion, 0)
    OPEN (121, STATUS='unknown', FILE='potFixed.'//outnam)
    WRITE (121, *) '# nc, nk = ', nc, nk
    WRITE (6, *) 'printing field...'
    CALL printfield(121, potfixedion, 'tp.potFix')
    CLOSE (121)
  ELSE IF (ipotfixed == -1) THEN
    DO i = 1, kdfull2
      potfixedion(i) = 0D0
      hfield(i) = 0D0
      hfield2(i) = 0D0
    END DO
    WRITE (6, *) 'Calculating potential of fixed ions. '
    WRITE (6, *) 'This may take a WHILE.'
    CALL addgsmpot(hfield, 1)
    CALL addgsmpot(potfixedion, 0)
    IF (idielec /= 0) THEN
      CALL addgsmpot2(hfield2, 0)
      DO ind = 1, kdfull2
        CALL conv1to3(ind)
        IF (iindtmp(1) <= (nint(xdielec/dx) + nx)) potfixedion(ind) = hfield2(ind)
      END DO
    END IF
    CALL addshortrepulsivepotonsubgrid(potfixedion, 0)
    WRITE (6, *) 'Done. Writing TO FILE...'
    OPEN (121, STATUS='unknown', FILE='potFixed.'//outnam)
    WRITE (121, *) '# nc, nk = ', nc, nk
    WRITE (6, *) 'printing field...'
    CALL printfield(121, potfixedion, 'tp.potFix')
    CALL printfield2(924, potfixedion, hfield)
    CLOSE (121)
    STOP 'STOP after calculating potential of fixed ions'
  ELSE IF (ipotfixed == 0) THEN
    DO i = 1, kdfull2
      potfixedion(i) = 0D0
    END DO
  ELSE
    STOP 'IN initSurface: invalid value for iPotFixed'
  END IF

  IF (isystem == 1) THEN ! the MgO CASE

! first determine lowest surface slab

    zmin = 1D9

    DO i = 1, nc
      IF (zc(i) < zmin) zmin = zc(i)
    END DO

    DO i = 1, nk
      IF (zk(i) < zmin) zmin = zk(i)
    END DO

    WRITE (6, *) 'minimaler z-Wert: ', zmin

    ind = 0

    z0 = (minz - nzsh)*dz
    be = 9.7D0
    ala = 3.98D-2
    rc = zmin - ala - 2/be - z0
    WRITE (6, *) 'rc = ', rc
    cf = 1.7D0/13.6D0

    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        DO ix = minx, maxx

          ind = ind + 1

          potfixedion(ind) = potfixedion(ind) - cf/(1 + EXP(be*(z1 - z0 - rc)))

        END DO
      END DO
    END DO

  END IF

  IF (ifredmas) THEN
    WRITE (6, *) 'REDUCING SURFACE MASS'
    mion = mion/23*0.5D0
    mkat = mkat/23*0.5D0
    me = me/23*0.5D0
  END IF

! define protocol files

  IF (myn == 0) THEN
    OPEN (121, STATUS='unknown', FILE='pposel.'//outnam)
    WRITE (121, *) '# t xe ye ze '
    OPEN (122, STATUS='unknown', FILE='pposani.'//outnam)
    WRITE (122, *) '# t x y z '
    OPEN (123, STATUS='unknown', FILE='penergy.'//outnam)
    WRITE (123, *) '# t ekinion ekinel ekinkat epots ecoul &
      eshort totalenergy'
    OPEN (124, STATUS='unknown', FILE='pforces.'//outnam)
    WRITE (124, *) '# t forx fory forz '
! OPEN(125,STATUS='unknown',FILE='pposkat.'//outnam)

    OPEN (226, STATUS='unknown', FILE='pdisplace.'//outnam)
! OPEN(127,STATUS='unknown',FILE='pgeomel.'//outnam)
! OPEN(147,STATUS='unknown',FILE='pgeomion.'//outnam)
! OPEN(149,STATUS='unknown',FILE='ptempsurf.'//outnam)
! OPEN(150,STATUS='unknown',FILE='ptempion.'//outnam)
! OPEN(151,STATUS='unknown',FILE='pkinenion.'//outnam)
! OPEN(152,STATUS='unknown',FILE='pkinencore.'//outnam)
! OPEN(153,STATUS='unknown',FILE='pkinenkat.'//outnam)
! OPEN(128,STATUS='unknown',FILE='ptherm.'//outnam)

    OPEN (130, STATUS='unknown', FILE='pstatus.'//outnam)
  END IF

  RETURN

END SUBROUTINE initsurface
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE iperiogsm
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! the electrostatic properties for the GSM particles
! allow for different types of anions (and cations) at the same
! time, but so far it makes no sense, because IN that CASE
! we would need more parameters for the short range interaction
  INTEGER :: ii

  DO ii = 1, nc
    IF (itypc(ii) == 1) THEN ! MgO CASE

      mion = 16D0 ! maybe skip the term IN brackets
! sigmac = 1.0
! chgc(ii) = 0.9345
! me = mrel * mion
! sigmav = 1.0
! chge(ii) = -2.9345
! cspr = 1.0

    ELSE IF (itypc(ii) == 2) THEN ! Ar CASE

      mion = 39.95D0
      sigmac = 1.43D0/sq2
      chgc(ii) = 6.119D0
      me = 4.38D0/1836D0
      sigmav = sigmac
      chge(ii) = -6.119D0
      cspr = -6.119D0**2*e2/11.08D0

    ELSE
      STOP 'This KIND of GSM ion is not yet implemented'
    END IF
  END DO

  DO ii = 1, nk
    IF (itypc(ii) == 1) THEN ! MgO CASE
      mkat = 24.31D0

    ELSE
      STOP 'This KIND of GSM cations is not yet implemented'
    END IF
  END DO

  WRITE (6, *) 'ALMOST END OF SURFACE INIT'

  IF (ifredmas) THEN
    WRITE (6, *) 'REDUCING SURFACE MASS'
    mion = mion/23*0.5D0
    mkat = mkat/23*0.5D0
    me = me/23*0.5D0
  END IF

  RETURN
END SUBROUTINE iperiogsm
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION getdistance2(ii, jj)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! returns distances**2 of particles with indices ii and jj
! WHERE ii,jj is from 1 TO nc+ne+nk
  INTEGER, INTENT(IN) :: ii
  INTEGER, INTENT(IN) :: jj
  INTEGER:: dis2

  IF (ii <= nc) THEN ! ii is core
    IF (jj <= nc) THEN ! jj is core
      dis2 = (xc(ii) - xc(jj))**2 + (yc(ii) - yc(jj))**2 + (zc(ii) - zc(jj))**2
    ELSE IF (jj <= nc + NE) THEN ! jj is v-cloud
      dis2 = (xc(ii) - xe(jj - nc))**2 + (yc(ii) - ye(jj - nc))**2 &
             + (zc(ii) - ze(jj - nc))**2
    ELSE ! jj is cation
      dis2 = (xc(ii) - xk(jj - nc - NE))**2 + (yc(ii) - yk(jj - nc - NE))**2 &
             + (zc(ii) - zk(jj - nc - NE))**2
    END IF
  ELSE IF (ii <= nc + NE) THEN ! ii is valence cloud
    IF (jj <= nc + NE) THEN ! jj is v-cloud
      dis2 = (xe(ii - nc) - xe(jj - nc))**2 + (ye(ii - nc) - ye(jj - nc))**2 &
             + (ze(ii - nc) - ze(jj - nc))**2
    ELSE ! jj is cation
      dis2 = (xe(ii - nc) - xk(jj - nc - NE))**2 + (ye(ii - nc) - yk(jj - nc - NE))**2 &
             + (ze(ii - nc) - zk(jj - nc - NE))**2
    END IF
  ELSE ! ii and jj are cations
    dis2 = (xk(ii - nc - NE) - xk(jj - nc - NE))**2 + (yk(ii - nc - NE) - yk(jj - nc - NE))**2 &
           + (zk(ii - nc - NE) - zk(jj - nc - NE))**2
  END IF

  getdistance2 = dis2

  RETURN
END FUNCTION getdistance2
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION getmixedwidth(ii, jj)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! returns mixed Gaussian widths between particles with indices ii and jj
! WHERE ii,jj are "long" indices from 1 TO nc+ne+nk
! some of the sigsig functions can be replaced by simple
! multiplication with SQ2 !!!
! actually all sigsigs can be avoided by replacing by
! stored variables
  INTEGER :: ii, jj
  REAL(DP) :: sigsig, sss

  IF (ii <= nc) THEN ! ii is core
    IF (jj <= nc) THEN ! jj is core

      sss = sigsig(sigmac, sigmac)

    ELSE IF (jj <= nc + NE) THEN ! jj is v-cloud

      sss = sigsig(sigmac, sigmav)

    ELSE ! jj is cation

      sss = sigsig(sigmac, sigmak)

    END IF
  ELSE IF (ii <= nc + NE) THEN ! ii is valence cloud
    IF (jj <= nc + NE) THEN ! jj is v-cloud

      sss = sigsig(sigmav, sigmav)

    ELSE ! jj is cation

      sss = sigsig(sigmav, sigmak)

    END IF
  ELSE ! ii and jj are cations

    sss = sigsig(sigmak, sigmak)

  END IF

  getmixedwidth = sss

  RETURN
END FUNCTION getmixedwidth
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getparas(ind)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! stores POSITION of substrate particule NUMBER "ind" IN COMMON vector rvectmp(1:3)
  INTEGER, INTENT(IN) :: ind

  IF (ind <= nc) THEN
    chgtmp = chgc(ind)
    sigtmp = sigmac
    rvectmp(1) = xc(ind)
    rvectmp(2) = yc(ind)
    rvectmp(3) = zc(ind)
    imobtmp = imobc(ind)
  ELSE IF (ind <= nc + NE) THEN
    chgtmp = chge(ind - nc)
    sigtmp = sigmav
    rvectmp(1) = xe(ind - nc)
    rvectmp(2) = ye(ind - nc)
    rvectmp(3) = ze(ind - nc)
    imobtmp = imobe(ind - nc)
  ELSE IF (ind <= nc + NE + nk) THEN
    chgtmp = chgk(ind - nc - NE)
    sigtmp = sigmak
    rvectmp(1) = xk(ind - nc - NE)
    rvectmp(2) = yk(ind - nc - NE)
    rvectmp(3) = zk(ind - nc - NE)
    imobtmp = imobk(ind - nc - NE)
  ELSE
    STOP 'TO be implemented IN getParas'
  END IF

  RETURN
END SUBROUTINE getparas
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION getcharge(ii)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
  INTEGER ::ii
  IF (ii <= nc) THEN
    getcharge = chgc(ii)
  ELSE IF (ii <= nc + NE) THEN
    getcharge = chge(ii - nc)
  ELSE
    getcharge = chgk(ii - nc - NE)
  END IF

  RETURN
END FUNCTION getcharge
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE addforce(i, forcx, forcy, forcz)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! adds force forcx,forcy,forcz TO correct array depending
! on "long" index i
  INTEGER :: i
  REAL(DP) :: forcx, forcy, forcz
  IF (i <= nc) THEN
    fxc(i) = fxc(i) + forcx
    fyc(i) = fyc(i) + forcy
    fzc(i) = fzc(i) + forcz
  ELSE IF (i <= nc + NE) THEN
    fxe(i - nc) = fxe(i - nc) + forcx
    fye(i - nc) = fye(i - nc) + forcy
    fze(i - nc) = fze(i - nc) + forcz
  ELSE
    fxk(i - nc - NE) = fxk(i - nc - NE) + forcx
    fyk(i - nc - NE) = fyk(i - nc - NE) + forcy
    fzk(i - nc - NE) = fzk(i - nc - NE) + forcz
  END IF

  RETURN
END SUBROUTINE addforce
!------------------------------------------------------------
!#endif

!------------------------------------------------------------

INTEGER FUNCTION iconvlongtoshort(ind1)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind1

  INTEGER :: ireturn

  IF (ind1 <= nc) THEN
    ireturn = ind1
  ELSE IF (ind1 <= nc + NE) THEN
    ireturn = ind1 - nc
  ELSE
    ireturn = ind1 - nc - NE
  END IF

  iconvlongtoshort = ireturn

  RETURN
END FUNCTION iconvlongtoshort
!------------------------------------------------------------

!------------------------------------------------------------

INTEGER FUNCTION iconvshorttolong(ityp1, ind1)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! converts "short" index ind1 of given particle TYPE ityp1
! TO "long" index
  INTEGER, INTENT(IN) :: ityp1
  INTEGER, INTENT(IN) :: ind1

  INTEGER :: ireturn

  IF (ityp1 == 1) THEN
    ireturn = ind1
  ELSE IF (ityp1 == 2) THEN
    ireturn = ind1 + nc
  ELSE IF (ityp1 == 3) THEN
    ireturn = ind1 + nc + NE
  ELSE
    STOP 'Error IN iconvShortToLong(...)'
  END IF

  iconvshorttolong = ireturn

  RETURN
END FUNCTION iconvshorttolong
!------------------------------------------------------------

! ------------------------------------------------------------

SUBROUTINE propagate
! ------------------------------------------------------------
  USE params
  IMPLICIT NONE
! propagation using velocity Verlet step
  INTEGER :: i, npar
  REAL(DP) :: eth, ek, resc
! at first propagate locations with old force
! DO ii = 1,nmovc
! i = ixyzcm(ii)
! rescale impulses IF cooling TO given temperature
  IF (icool == 1 .AND. itindex <= icoolsteps .AND. &
      MOD(itindex, icoolevry) == 0) THEN
    npar = nc + NE + nk
    eth = 1.5D0*npar*temperature
    ek = 0.0D0
    DO i = 1, nc
      ek = ek + (pxc(i)**2 + pyc(i)**2 + pzc(i)**2)/2D0/mion
    END DO
    DO i = 1, NE
      ek = ek + (pxe(i)**2 + pye(i)**2 + pze(i)**2)/2D0/me
    END DO
    DO i = 1, nc
      ek = ek + (pxk(i)**2 + pyk(i)**2 + pzk(i)**2)/2D0/mkat
    END DO
    resc = SQRT(eth/ek)
    WRITE (6, *) 'resc=', resc

    DO i = 1, nc
! WRITE(6,*) pxc(i)
      pxc(i) = pxc(i)*resc
      pyc(i) = pyc(i)*resc
      pzc(i) = pzc(i)*resc
      pxcold(i) = pxcold(i)*resc
      pycold(i) = pycold(i)*resc
      pzcold(i) = pzcold(i)*resc
    END DO
    DO i = 1, NE
      pxe(i) = pxe(i)*resc
      pye(i) = pye(i)*resc
      pze(i) = pze(i)*resc
      pxeold(i) = pxeold(i)*resc
      pyeold(i) = pyeold(i)*resc
      pzeold(i) = pzeold(i)*resc
    END DO
    DO i = 1, nk
      pxk(i) = pxk(i)*resc
      pyk(i) = pyk(i)*resc
      pzk(i) = pzk(i)*resc
      pxkold(i) = pxkold(i)*resc
      pykold(i) = pykold(i)*resc
      pzkold(i) = pzkold(i)*resc
    END DO

  END IF

  DO i = 1, nc

    xc(i) = xcold(i) + dte*pxcold(i)/mion + 0.5D0*dte**2*fxcold(i)/mion
    yc(i) = ycold(i) + dte*pycold(i)/mion + 0.5D0*dte**2*fycold(i)/mion
    zc(i) = zcold(i) + dte*pzcold(i)/mion + 0.5D0*dte**2*fzcold(i)/mion

  END DO

! DO ii = 1,nmove

! i = ixyzem(ii)

  DO i = 1, NE

    xe(i) = xeold(i) + dte*pxeold(i)/me + 0.5D0*dte**2*fxeold(i)/me
    ye(i) = yeold(i) + dte*pyeold(i)/me + 0.5D0*dte**2*fyeold(i)/me
    ze(i) = zeold(i) + dte*pzeold(i)/me + 0.5D0*dte**2*fzeold(i)/me

  END DO

! DO ii=1,nmovk
! i = ixyzkm(ii)

  DO i = 1, nk

    xk(i) = xkold(i) + dte*pxkold(i)/mkat + 0.5D0*dte**2*fxkold(i)/mkat
    yk(i) = ykold(i) + dte*pykold(i)/mkat + 0.5D0*dte**2*fykold(i)/mkat
    zk(i) = zkold(i) + dte*pzkold(i)/mkat + 0.5D0*dte**2*fzkold(i)/mkat
  END DO

! locations propagated, now calculate new force

! CALL getForces(itime,rho,psi,iflag)

! now propagate impulses with new force

! DO ii=1,nmovc

! i = ixyzcm(ii)

  DO i = 1, nc

    pxc(i) = pxcold(i) + 0.5D0*dte*(fxcold(i) + fxc(i))
    pyc(i) = pycold(i) + 0.5D0*dte*(fycold(i) + fyc(i))
    pzc(i) = pzcold(i) + 0.5D0*dte*(fzcold(i) + fzc(i))

    fxcold(i) = fxc(i)
    fycold(i) = fyc(i)
    fzcold(i) = fzc(i)
    pxcold(i) = pxc(i)
    pycold(i) = pyc(i)
    pzcold(i) = pzc(i)
    xcold(i) = xc(i)
    ycold(i) = yc(i)
    zcold(i) = zc(i)

  END DO

! DO ii=1,nmove

! i = ixyzem(ii)

  DO i = 1, NE

    pxe(i) = pxeold(i) + 0.5D0*dte*(fxeold(i) + fxe(i))
    pye(i) = pyeold(i) + 0.5D0*dte*(fyeold(i) + fye(i))
    pze(i) = pzeold(i) + 0.5D0*dte*(fzeold(i) + fze(i))

    fxeold(i) = fxe(i)
    fyeold(i) = fye(i)
    fzeold(i) = fze(i)
    pxeold(i) = pxe(i)
    pyeold(i) = pye(i)
    pzeold(i) = pze(i)
    xeold(i) = xe(i)
    yeold(i) = ye(i)
    zeold(i) = ze(i)

  END DO ! END of prop. for ions loop

! DO ii=1,nmovk

! i = ixyzkm(ii)

  DO i = 1, nk

    pxk(i) = pxkold(i) + 0.5D0*dte*(fxkold(i) + fxk(i))
    pyk(i) = pykold(i) + 0.5D0*dte*(fykold(i) + fyk(i))
    pzk(i) = pzkold(i) + 0.5D0*dte*(fzkold(i) + fzk(i))
    fxkold(i) = fxk(i)
    fykold(i) = fyk(i)
    fzkold(i) = fzk(i)
    pxkold(i) = pxk(i)
    pykold(i) = pyk(i)
    pzkold(i) = pzk(i)
    xkold(i) = xk(i)
    ykold(i) = yk(i)
    zkold(i) = zk(i)

  END DO

! update subgrids IF necessary

  CALL updatesubgrids

  CALL calctemp

END SUBROUTINE propagate
! ------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE printpos
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER :: i
  REAL(DP) :: t

  t = dte*itindex

  DO i = 1, nc
    WRITE (121, '(1f10.3,6e17.5)') t, xe(i), ye(i), ze(i), pxe(i), pye(i), pze(i)
    WRITE (124, '(1f10.3,3e17.5)') t, fxcold(i), fycold(i), fzcold(i)
    WRITE (122, '(1f10.3,6e17.5)') t, xc(i), yc(i), zc(i), pxc(i), pyc(i), pzc(i)
  END DO

  DO i = 1, nk
    WRITE (125, '(1f10.3,6e17.5)') t, xk(i), yk(i), zk(i), pxk(i), pyk(i), pzk(i)
  END DO

END SUBROUTINE printpos
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE printdipoles(ishift)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! PRINT shell displacements
  INTEGER, INTENT(IN) :: ishift
  INTEGER :: i
  REAL(DP) :: r, t, xi, yi, zi
  t = dte*itindex

  DO i = 1, NE
    r = SQRT((xc(i) - xe(i))**2 + (yc(i) - ye(i))**2 + (zc(i) - ze(i))**2)
    xi = xc(i) - xe(i)
    yi = yc(i) - ye(i)
    zi = zc(i) - ze(i)

    WRITE (226 + ishift, '(1f10.4,4e15.5)') t, r, xi, yi, zi
  END DO

END SUBROUTINE printdipoles
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE cool
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER :: i

  DO i = 1, nc
    pxc(i) = pxc(i)
    pyc(i) = pyc(i)
    pzc(i) = pzc(i)
    pxcold(i) = pxcold(i)
    pycold(i) = pycold(i)
    pzcold(i) = pzcold(i)
  END DO
  DO i = 1, NE
    pxe(i) = pxe(i)*0.01D0
    pye(i) = pye(i)*0.01D0
    pze(i) = pze(i)*0.01D0
    pxeold(i) = pxeold(i)*0.01D0
    pyeold(i) = pyeold(i)*0.01D0
    pzeold(i) = pzeold(i)*0.01D0

  END DO
  DO i = 1, nk
    pxk(i) = pxk(i)*0.01D0
    pyk(i) = pyk(i)*0.01D0
    pzk(i) = pzk(i)*0.01D0
    pxkold(i) = pxkold(i)*0.01D0
    pykold(i) = pykold(i)*0.01D0
    pzkold(i) = pzkold(i)*0.01D0

  END DO

END SUBROUTINE cool
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE calctemp
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! calculates thermal energies
  INTEGER :: i
  REAL(DP) :: sumc, sume, sumk, t, tempc, tempe, tempk, temper

  sumc = 0D0
  sume = 0D0
  sumk = 0D0

  DO i = 1, nc
    sumc = sumc + 0.5D0*(pxc(i)*pxc(i) + pyc(i)*pyc(i) + pzc(i)*pzc(i))/mion
  END DO
  DO i = 1, NE
    sume = sume + 0.5D0*(pxe(i)*pxe(i) + pye(i)*pye(i) + pze(i)*pze(i))/me
  END DO
  DO i = 1, nk
    sumk = sumk + 0.5D0*(pxk(i)*pxk(i) + pyk(i)*pyk(i) + pzk(i)*pzk(i))/mkat
  END DO

! temp is not T but rather kT

  ttempc = ttempc + 2D0/3D0*sumc/nc**2
  ttempe = ttempe + 2D0/3D0*sume/NE**2
  ttempk = ttempk + 2D0/3D0*sumk/nc**2
  ttemp = ttemp + 2D0/3D0*(sumc + sume + sumk)/(nc + NE + nk)**2

  IF (MOD(itindex, 10) == 0) THEN
    tempc = ttempc/10
    tempe = ttempe/10
    tempk = ttempk/10
    temper = ttemp/10
    ttempc = 0D0
    ttempe = 0D0
    ttempk = 0D0
    ttemp = 0D0
  END IF

  IF (MOD(itindex, ithermevry) == 0) THEN
    WRITE (128, '(1f15.5,5e15.5)') t, tempc, tempe, tempk, temper, &
      temp*1.5789D5*(nc + NE + nk)
  END IF

END SUBROUTINE calctemp
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE calcenergyg
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER :: i, j, n
  REAL(DP) :: r, t, ek, epotel, epotion, esh

  ekinion = 0D0
  epotion = 0D0
  ekinel = 0D0
  epotel = 0D0
  ekinkat = 0D0
  eshort = 0D0

! WRITE(6,*) n,ne,nk

  t = dte*itindex

  DO i = 1, nc
    ekinion = ekinion + (pxc(i)*pxc(i) + pyc(i)*pyc(i) + pzc(i)*pzc(i))/2/mion
  END DO

  DO i = 1, NE
    ekinel = ekinel + (pxe(i)*pxe(i) + pye(i)*pye(i) + pze(i)*pze(i))/2/me
  END DO

  DO i = 1, nk
    ekinkat = ekinkat + (pxk(i)*pxk(i) + pyk(i)*pyk(i) + pzk(i)*pzk(i))/2/mkat
  END DO

! potential energies:

  epots = 0D0
! springs
  DO i = 1, NE
    epots = epots + 0.5D0*cspr*((xc(i) - xe(i))**2 + &
                                (yc(i) - ye(i))**2 + (zc(i) - ze(i))**2)
  END DO

! Coulomb interaction energy
  ecoul = 0D0
! ion-ion
  DO i = 1, nc
    DO j = 1, nc
      IF (j /= i) THEN
        r = SQRT((xc(i) - xc(j))**2 + (yc(i) - yc(j))**2 + (zc(i) - zc(j))**2)
        IF (r < crcut) THEN
          r = crcut

        END IF

        ecoul = ecoul + 0.5D0*chgc(i)*chgc(j)*e2*erf(r/sigmacc)/r

      END IF
    END DO
  END DO

! ion-cation
! something's wrong here
  DO i = 1, nc

    DO j = 1, nk
      r = SQRT((xc(i) - xk(j))**2 + (yc(i) - yk(j))**2 + (zc(i) - zk(j))**2)
      IF (r < crcut) THEN
        r = crcut

      END IF

      ecoul = ecoul + chgc(i)*chgk(j)*e2*erf(r/sigmack)/r

    END DO
! i-k done
  END DO

! ion-(other)electron
  DO i = 1, n
    DO j = 1, NE
      IF (j /= i) THEN
        r = SQRT((xc(i) - xe(j))**2 + (yc(i) - ye(j))**2 + (zc(i) - ze(j))**2)
        IF (r < crcut) THEN
          r = crcut

        END IF

! WRITE(6,*) ecoul
! IF (itindex.gt.0) STOP

        ecoul = ecoul + chgc(i)*chge(j)*e2*erf(r/sigmacv)/r

      END IF
    END DO
! i-e done
  END DO

! electron-electron
  DO i = 1, NE
    DO j = 1, NE
      IF (i /= j) THEN
        r = SQRT((xe(i) - xe(j))**2 + (ye(i) - ye(j))**2 + (ze(i) - ze(j))**2)
        IF (r < crcut) THEN
          r = crcut

        END IF

        ecoul = ecoul + 0.5D0*chge(i)*chge(j)*e2*erf(r/sigmavv)/r
      END IF
    END DO
! e-e done
  END DO

! electron-cation
  DO i = 1, NE
    DO j = 1, nk
      r = SQRT((xe(i) - xk(j))**2 + (ye(i) - yk(j))**2 + (ze(i) - zk(j))**2)
      IF (r < crcut) THEN
        r = crcut

      END IF

      ecoul = ecoul + chge(i)*chgk(j)*e2*erf(r/sigmakv)/r
    END DO
! e-k done
  END DO

! cation-cation
  DO i = 1, nk
    DO j = 1, nk
      IF (i /= j) THEN
        r = SQRT((xk(i) - xk(j))**2 + (yk(i) - yk(j))**2 + (zk(i) - zk(j))**2)
        IF (r < crcut) THEN
          r = crcut

        END IF

        ecoul = ecoul + 0.5D0*chgk(i)*chgk(j)*e2*erf(r/sigmakk)/r

      END IF
    END DO
  END DO

  IF (ifmdshort == 1) THEN
! now short range interactions

    DO i = 1, nk
      DO j = 1, NE
        r = SQRT((xe(j) - xk(i))**2 + (ye(j) - yk(i))**2 + (ze(j) - zk(i))**2)
        esh = bkv*EXP(-r/sigkv) + ckv/r**6
        eshort = eshort + esh
      END DO
    END DO

    DO i = 1, NE
      DO j = 1, NE
        IF (i /= j) THEN
          r = SQRT((xe(j) - xe(i))**2 + (ye(j) - ye(i))**2 + (ze(j) - ze(i))**2)
          esh = bvv*EXP(-r/sigvv) + cvv/r**6
          eshort = eshort + 0.5D0*esh
        END IF
      END DO
    END DO

    DO i = 1, nk
      DO j = 1, nk
        IF (i /= j) THEN
          r = SQRT((xk(j) - xk(i))**2 + (yk(j) - yk(i))**2 + (zk(j) - zk(i))**2)
          esh = bkk*EXP(-r/sigkk) + ckk/r**6
          eshort = eshort + 0.5D0*esh
        END IF
      END DO
    END DO

  END IF

  IF (eshort <= 1D-40) THEN
    eshort = 0D0
  END IF
  IF (ekinel <= 1D-40) THEN
    ekinel = 0D0
  END IF
  IF (ekinkat <= 1D-40) THEN
    ekinkat = 0D0
  END IF
  IF (ekinion <= 1D-40) THEN
    ekinion = 0D0
  END IF

  energyg = ekinion + ekinel + ekinkat + epots + ecoul + eshort

  ek = ekinel + ekinion + ekinkat

  IF (myn == 0) THEN
    WRITE (123, '(1f10.2,7e20.8)') t, ekinion, ekinel, ekinkat &
      , epots, ecoul, eshort, energyg
  END IF

END SUBROUTINE calcenergyg
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE initimpulses(INDEX)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: INDEX
  INTEGER :: i

  IF (INDEX == 0) THEN

    DO i = 1, nc
      pxcold(i) = 0D0
      pycold(i) = 0D0
      pzcold(i) = 0D0
      pxc(i) = 0D0
      pyc(i) = 0D0
      pzc(i) = 0D0
      pxeold(i) = 0D0
      pyeold(i) = 0D0
      pzeold(i) = 0D0
      pxe(i) = 0D0
      pye(i) = 0D0
      pze(i) = 0D0
    END DO

    DO i = 1, nk
      pxkold(i) = 0D0
      pykold(i) = 0D0
      pzkold(i) = 0D0
      pxk(i) = 0D0
      pyk(i) = 0D0
      pzk(i) = 0D0
    END DO

  ELSE IF (INDEX > 0) THEN

    ! nothing ?? F.L.

  ELSE
    STOP 'Could not set initial impulses'

  END IF

END SUBROUTINE initimpulses
!------------------------------------------------------------

!~ !------------------------------------------------------------

!~ SUBROUTINE disturblattice
!~ !------------------------------------------------------------
!~ USE params
!~ IMPLICIT NONE

!~ LOGICAL :: tterm
!~ INTEGER :: i
!~ DATA tterm/.true./

!~ IF(tterm) RETURN

!~ DO i=1,nc
!~ ! dist = (drand(0)-0.5D0)*maxDist
!~ xc(i) = xc(i) + dist
!~ xcold(i) = xcold(i) + dist

!~ ! dist = (drand(0)-0.5D0)*maxDist
!~ yc(i) = yc(i) + dist
!~ ycold(i) = ycold(i) + dist

!~ ! dist = (drand(0)-0.5D0)*maxDist
!~ zc(i) = zc(i) + dist
!~ zcold(i) = zcold(i) + dist
!~ END DO

!~ DO i=1,NE
!~ ! dist = (drand(0)-0.5D0)*maxDist
!~ xe(i) = xe(i) + dist
!~ xeold(i) = xeold(i) + dist

!~ ! dist = (drand(0)-0.5D0)*maxDist
!~ ye(i) = ye(i) + dist
!~ yeold(i) = yeold(i) + dist

!~ ! dist = (drand(0)-0.5D0)*maxDist
!~ ze(i) = ze(i) + dist
!~ zeold(i) = zeold(i) + dist
!~ END DO

!~ DO i=1,nk
!~ ! dist = (drand(0)-0.5D0)*maxDist
!~ xk(i) = xk(i) + dist
!~ xkold(i) = xkold(i) + dist

!~ ! dist = (drand(0)-0.5D0)*maxDist
!~ yk(i) = yk(i) + dist
!~ ykold(i) = ykold(i) + dist

!~ ! dist = (drand(0)-0.5D0)*maxDist
!~ zk(i) = zk(i) + dist
!~ zkold(i) = zkold(i) + dist
!~ END DO

!~ RETURN

!~ END SUBROUTINE disturblattice
!~ !------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE printoutparameters
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER :: i

  WRITE (6, *) '*************************************************'

  WRITE (6, '(a,3i6)') 'nc, ne, nk: ', nc, NE, nk
  WRITE (6, '(a,2i6)') 'nion, nelect: ', nion, nelect

  WRITE (6, *) '*************************************************'

  WRITE (6, *) '*'
  WRITE (6, *) 'cspr = ', cspr
  WRITE (6, *) '*'

  WRITE (6, *) '************************************************'

  WRITE (6, *) 'Gaussian widths:'

  WRITE (6, '(a,3e17.7)') 'sigmac, sigmav, sigmak: ', sigmac, sigmav, sigmak
  WRITE (6, '(a,3e17.7)') 'sigmacc, sigmacv, sigmack: ', &
    sigmacc, sigmacv, sigmack
  WRITE (6, '(a,3e17.7)') 'sigmavv, sigmakv, sigmakk: ', &
    sigmavv, sigmakv, sigmakk

  WRITE (6, *) '*************************************************'

  WRITE (6, *) '*'
  WRITE (6, *) '*'

  WRITE (6, *) 'POSITION of Na: '

  DO i = 1, nion
    WRITE (6, '(1i6,3f15.4)') i, cx(i), cy(i), cz(i)
  END DO

  WRITE (6, *) '*************************************************'

  WRITE (6, *) '*'
  WRITE (6, *) '*'

  WRITE (6, *) 'POSITION of Cores: '

  DO i = 1, nc
    WRITE (6, '(1i6,3f15.4)') i, xc(i), yc(i), zc(i)
  END DO

  WRITE (6, *) '*************************************************'

  WRITE (6, *) '*'
  WRITE (6, *) '*'

  WRITE (6, *) 'POSITION of Valences: '

  DO i = 1, NE
    WRITE (6, '(1i6,3f15.4)') i, xe(i), ye(i), ze(i)
  END DO

  WRITE (6, *) '*************************************************'

  WRITE (6, *) '*'
  WRITE (6, *) '*'

  WRITE (6, *) 'POSITION of cations: '

  DO i = 1, nk
    WRITE (6, '(1i6,3f15.4)') i, xk(i), yk(i), zk(i)
  END DO

  WRITE (6, *) '*************************************************'

  WRITE (6, *) '*'
  WRITE (6, *) '*'

END SUBROUTINE printoutparameters
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION rgetmaxpol()
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
  INTEGER :: i
  REAL(DP) :: rr
  rgetmaxpol = -1D0

  DO i = 1, nc
    rr = (xc(i) - xe(i))**2 + (yc(i) - ye(i))**2 + (zc(i) - ze(i))**2
    rr = SQRT(rr)
    IF (rr > rgetmaxpol) rgetmaxpol = rr
  END DO

  RETURN
END FUNCTION rgetmaxpol
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE rgetmaxforce(iflag)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iflag

  INTEGER :: i
  REAL(DP) :: rgetmaxx, rgetmaxy, rgetmaxz

  rgetmaxx = -1D60
  rgetmaxy = -1D60
  rgetmaxz = -1D60

  IF (iflag == 2) THEN
    DO i = 1, NE
      IF (imobe(i) == 1) THEN
        IF (ABS(fxe(i)) > rgetmaxx) rgetmaxx = fxe(i)
        IF (ABS(fye(i)) > rgetmaxy) rgetmaxy = fye(i)
        IF (ABS(fze(i)) > rgetmaxz) rgetmaxz = fze(i)
      END IF
    END DO
  END IF

  rvectmp(1) = rgetmaxx
  rvectmp(2) = rgetmaxy
  rvectmp(3) = rgetmaxz

  RETURN
END SUBROUTINE rgetmaxforce
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION rgetmeanzpol()
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER :: i, ic
  REAL(DP) :: rr

  rr = 0D0
  ic = 0

  DO i = 1, nc
    IF (imobe(i) == 1) THEN

      rr = zc(i) - ze(i) + rr
      ic = ic + 1
    END IF
  END DO

  rr = rr/ic ! Is it possible that ic = 0 ??

  rgetmeanzpol = rr

  RETURN
END FUNCTION rgetmeanzpol
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE pripolariz
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! prints information on polarization of substrate

  INTEGER :: i, icou, icurlay, ii, iter1
  REAL(DP) :: dmspolx, dmspoly, dmspolz
  REAL(DP) :: polav, polmin, polmax, polstddev, polz
  REAL(DP) :: r, rr, rrr, rrre
  REAL(DP) :: polavl(maxnlayers)
  REAL(DP) :: polstddevl(maxnlayers)
  REAL(DP) :: polminl(maxnlayers)
  REAL(DP) :: polmaxl(maxnlayers)

! CALL printoutparameters

  IF (iaxis == 0) THEN

    WRITE (6, *) nc, NE, nk

    polav = 0D0
    polstddev = 0D0
    polmin = 1D6
    polmax = 0D0

    icou = 0

    dmspolx = 0D0
    dmspoly = 0D0
    dmspolz = 0D0

    DO i = 1, nc

      rr = (xc(i) - xe(i))**2 + (yc(i) - ye(i))**2 + (zc(i) - ze(i))**2
      r = SQRT(rr)

      IF (imobe(i) == 1) THEN

        icou = icou + 1

        dmspolx = (xc(i) - xe(i))**2 + dmspolx
        dmspoly = (yc(i) - ye(i))**2 + dmspoly
        dmspolz = (zc(i) - ze(i))**2 + dmspolz

        polav = polav + r
        polstddev = polstddev + rr

        IF (r < polmin) polmin = r
        IF (r > polmax) polmax = r

      END IF

    END DO

    dmspolx = SQRT(dmspolx/icou)
    dmspoly = SQRT(dmspoly/icou)
    dmspolz = SQRT(dmspolz/icou)

! polAv = polAv/nc
    polav = polav/icou
! polStdDev = polStdDev/nc - polAv**2
    polstddev = polstddev/icou - polav**2

! OPEN(271,STATUS='unknown',FILE='ppolariz.'//outnam)

    WRITE (271, *) '#polarization of mob. substr. at iter1= ', iter1
    WRITE (271, *) '#'
    WRITE (271, *) '#mean polarization: ', polav
    WRITE (271, *) '#standard deviance: ', polstddev
    WRITE (271, *) '#minimum polarization: ', polmin
    WRITE (271, *) '#maximum polarization: ', polmax

    WRITE (6, *) '#polarization of substrate at iter1= ', iter1
    WRITE (6, *) '#'
    WRITE (6, *) '#mean polarization: ', polav
    WRITE (6, *) '##################################'
    WRITE (6, '(a,3e17.7)') 'dmspolx,y,z: ', dmspolx, dmspoly, dmspolz

    WRITE (6, *) '##################################'
    WRITE (6, *) '#standard deviance: ', polstddev
    WRITE (6, *) '#minimum polarization: ', polmin
    WRITE (6, *) '#maximum polarization: ', polmax

    DO ii = 1, nc
      rr = (xc(ii) - xe(ii))**2 + (yc(ii) - ye(ii))**2 + (zc(ii) - ze(ii))**2
      rr = SQRT(rr)
      rrr = xc(ii)**2 + yc(ii)**2
      rrr = SQRT(rrr)
      rrre = xe(ii)**2 + ye(ii)**2
      rrre = SQRT(rrre)

      polz = (zc(ii) - ze(ii))*zc(ii)/ABS(zc(ii))

! WRITE(271,'(4f17.7,2e17.7)') rrr,rrre,zc(ii),ze(ii),rr,polz
      WRITE (271, '(6e13.5)') xc(ii), yc(ii), zc(ii), &
        xc(ii) - xe(ii), yc(ii) - ye(ii), zc(ii) - ze(ii)

! rrr,rrre,zc(ii),ze(ii),rr,polz
    END DO

! CLOSE(271)

  ELSE

    DO icurlay = 1, nlay

      polav = 0D0
      polstddev = 0D0
      polmin = 1D6
      polmax = 0D0

      DO i = 1, nc

        IF (ilayerc(i) == icurlay) THEN
          rr = (xc(i) - xe(i))**2 + (yc(i) - ye(i))**2 + (zc(i) - ze(i))**2
          r = SQRT(rr)

          polav = polav + r
          polstddev = polstddev + rr

          IF (r < polmin) polmin = r
          IF (r > polmax) polmax = r

        END IF

      END DO

      polavl(icurlay) = polav/nc
      polstddevl(icurlay) = polstddev/nc - polav**2
      polminl(icurlay) = polmin
      polmaxl(icurlay) = polmax

    END DO

    OPEN (271, STATUS='unknown', FILE='ppolariz.'//outnam)

    WRITE (271, *) '#layer, mean pol.,std.dev.,min. pol, max. pol.'

    DO icurlay = 1, nlay

      WRITE (271, '(i6,4e17.7)') icurlay, polavl(icurlay), polstddevl(icurlay), &
        polminl(icurlay), polmaxl(icurlay)

    END DO

    CLOSE (271)

  END IF

  RETURN
END SUBROUTINE pripolariz
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE init_surftemp()

! initializes temperature of substrate constituents (MgO)

  USE params
  USE util, ONLY: givetemperature
  IMPLICIT NONE

  INTEGER :: i

  CALL givetemperature(pxc, pyc, pzc, nc, surftemp, mion*1836D0*ame, 1)
  DO i = 1, NE
    pxe(i) = me*1836D0*ame/mion*pxc(i)
    pye(i) = me*1836D0*ame/mion*pyc(i)
    pze(i) = me*1836D0*ame/mion*pzc(i)
  END DO
  CALL givetemperature(pxk, pyk, pzk, nk, surftemp, mkat*1836D0*ame, 3)

  RETURN
END SUBROUTINE init_surftemp

!------------------------------------------------------------

SUBROUTINE printsurfpot(iunit)
!------------------------------------------------------------
  USE params
  USE util, ONLY: printfield2
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iunit

  INTEGER :: i, ii, ind, indmin, indmax, ipoints
  REAL(DP) :: acc, du, pneg, pt, ptmin, ptmax, rmaxd, u, zhigh, zlow
  REAL(DP) :: field1(kdfull2), field2(kdfull2)
  REAL(DP) :: uhisto(501)
  REAL(DP), EXTERNAL :: getxval, getyval, getzval
  rmaxd = 15D0
! define interior:

  zlow = -12D0
  zhigh = -7.5D0

  DO i = 1, kdfull2
    field1(i) = 0D0
    field2(i) = 0D0
  END DO

  DO i = 1, 501
    uhisto(i) = 0D0
  END DO

  CALL addgsmpot(field1, 1)
  CALL addgsmpot(field1, 0)
  CALL addshortrepulsivepotonsubgrid(field2, 0)
  CALL addshortrepulsivepotonsubgrid(field2, 1)

  CALL printfield2(iunit, field1, field2)

! calculate some characteristic quantities

! mean value

  acc = 0D0
  ptmin = 1D9
  ptmax = -1D9
  indmax = 0
  indmin = 0
  ipoints = 0
  pneg = 0D0

  DO ii = 1, kdfull2

    IF (getxval(ii)**2 + getyval(ii)**2 <= rmaxd**2 .AND. &
        getzval(ii) > zlow .AND. getzval(ii) < zhigh) THEN

      ipoints = ipoints + 1

      pt = -field1(ii) - field2(ii)
      acc = acc + pt

      IF (pt < 0D0) pneg = pneg + 1

      IF (pt > ptmax) THEN
        ptmax = pt
        indmax = ii
      END IF
      IF (pt < ptmin) THEN
        ptmin = pt
        indmin = ii
      END IF

    END IF

  END DO
  acc = acc/ipoints

  WRITE (6, *) 'NUMBER of gridpoints inside = ', ipoints
  WRITE (iunit + 5, *) 'NUMBER of gridpoints inside = ', ipoints

  WRITE (6, *) 'rmaxd = ', rmaxd
  WRITE (iunit + 5, *) '#rmaxd = ', rmaxd

  WRITE (6, *) '% negative values: ', pneg/ipoints
  WRITE (iunit + 5, *) '#% negative values: ', pneg/ipoints

  WRITE (6, '(a,f15.5)') 'Average surface potential: ', acc
  WRITE (iunit + 5, '(a,f15.5)') '#Average surface potential: ', acc

  WRITE (6, '(a,e15.5,a,3f10.2)') &
    'maximum potential: ', ptmax, ' at ', getxval(indmax), &
    getyval(indmax), getzval(indmax)

  WRITE (iunit + 5, '(a,e15.5,a,3f10.2)') &
    '#maximum potential: ', ptmax, ' at ', getxval(indmax), &
    getyval(indmax), getzval(indmax)

  WRITE (6, '(a,e15.5,a,3f10.2)') &
    'minimum potential: ', ptmin, ' at ', getxval(indmin), &
    getyval(indmin), getzval(indmin)

  WRITE (iunit + 5, '(a,e15.5,a,3f10.2)') &
    '#minimum potential: ', ptmin, ' at ', getxval(indmin), &
    getyval(indmin), getzval(indmin)

  WRITE (iunit + 5, *) '# Histogram: '

! make histogram

  du = (ptmax - ptmin)/250

  DO ind = 1, kdfull2
    u = -field1(ind) - field2(ind)
    i = nint((u - ptmin)/du) + 1
    uhisto(i) = uhisto(i) + 1
  END DO

  DO i = 0, 250
    u = i*du + ptmin
    WRITE (iunit + 5, '(e17.7,e18.7)') u, uhisto(i + 1)
  END DO

  RETURN
END SUBROUTINE printsurfpot
!------------------------------------------------------------
#endif
