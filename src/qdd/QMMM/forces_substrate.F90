
!------------------------------------------------------------

#if(raregas)
SUBROUTINE getforceelgsm(iflag, rho)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iflag
  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)

  INTEGER :: ibegin, iend, ii, ii2, ind
  REAL(DP) :: prefc
  REAL(DP) :: rhotmp(2*kdfull2)

  INTEGER, EXTERNAL :: iconvlongtoshort, iptyp, isoutofbox
  REAL(DP), EXTERNAL :: v_soft, gauss, dgaussdr

! first Coulomb force from cluster-el density
! using pseudo-densities IF desired

  IF (nelect > 0) THEN

    IF (iflag /= 0) THEN ! ONLY force on valences
      ibegin = nc + 1
      iend = nc + NE
    ELSE ! forces on all particles
      ibegin = 1
      iend = nc + NE + nk
    END IF

    DO ii = ibegin, iend

      CALL getparas(ii)
      IF (.NOT. ipseudo) THEN
        IF (imobtmp /= 0) THEN
          IF (idielec == 1) THEN
            DO ind = 1, 2*kdfull2
              rhotmp(ind) = rho(ind)
            END DO
            CALL addimage(rho, 1)
          END IF
          CALL foldgradfunc(rho, v_soft, rvectmp(1), rvectmp(2), &
                            rvectmp(3), sq2*sigtmp)
! CALL dIntFieldFunc(rho,dVsdr,rVecTmp(1),rVecTmp(2),
! & rVecTmp(3),SQ2*sigTmp)
          IF (idielec == 1) THEN
            DO ind = 1, 2*kdfull2
              rho(ind) = rhotmp(ind)
            END DO
          END IF

          CALL addforce(ii, rvectmp(1)*e2*chgtmp, &
                        rvectmp(2)*e2*chgtmp, rvectmp(3)*e2*chgtmp)

! the plus sign is correct because rho is positive
! and we actually need electronic charge density
! which is negative...so its minus sign cancels
! with the minus sign from - nabla V

        END IF
      ELSE ! ipseudo = .TRUE.

        IF (imobtmp /= 0) THEN
          IF (isoutofbox(rvectmp(1), rvectmp(2), rvectmp(3)) == 2) CYCLE

          CALL foldgradfunconsubgrid(chpcoul, gauss, rvectmp(1), &
                                     rvectmp(2), rvectmp(3), sigtmp*sq2)
! CALL dIntFieldFuncOnSubGrid(chpcoul,dgaussdr,rVecTmp(1),
! & rVecTmp(2),rVecTmp(3),sigTmp*SQ2)

          prefc = chgtmp/(2*pi*sigtmp**2)**1.5D0
! no factor e2, because it is already IN the
! field chpcoul()
          CALL addforce(ii, prefc*rvectmp(1), prefc*rvectmp(2), prefc*rvectmp(3))
        END IF
      END IF

      ii2 = iconvlongtoshort(ii)

      CALL getshortforce(iptyp(ii), 5, ii2, 0, rho, 0)
! last PARAMETER=0 is a dummy PARAMETER because it is not needed
! here

    END DO

  END IF ! nelect

  RETURN
END SUBROUTINE getforceelgsm
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getforcenagsm(iflag)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iflag

  INTEGER :: ibegin, iend, ii, ii2, jj
  REAL(DP) :: dist, radfor, smix1, smix2, sigsig, xr, yr, zr
  REAL(DP) :: dummy(kdfull2*2)

  INTEGER, EXTERNAL :: iconvlongtoshort, iptyp, isoutofbox
  REAL(DP), EXTERNAL :: dv_softdr

  IF (iflag /= 0) THEN ! ONLY force on valences
    ibegin = nc + 1
    iend = nc + NE
  ELSE ! either forces on all particles or all but valences
    ibegin = 1
    iend = nc + NE + nk
  END IF

  DO ii = ibegin, iend

    CALL getparas(ii)

    IF (isoutofbox(rvectmp(1), rvectmp(2), rvectmp(3)) == 2) CYCLE

    DO jj = 1, nion

      xr = rvectmp(1) - cx(jj)
      yr = rvectmp(2) - cy(jj)
      zr = rvectmp(3) - cz(jj)

      dist = xr*xr + yr*yr + zr*zr
      dist = SQRT(dist)

      smix1 = sigsig(sigtmp, sgm1(np(jj)))
      smix2 = sigsig(sigtmp, sgm2(np(jj)))

      radfor = e2*chgtmp*(chg1(np(jj))*dv_softdr(dist, smix1*sq2) + &
                          chg2(np(jj))*dv_softdr(dist, smix2*sq2))/dist

! radfor = e2*chgTmp*(chg1(np(jj))*dVsdr(dist,smix1*SQ2)
! & + chg2(np(jj))*dVsdr(dist,smix2*SQ2))/dist

      forx = xr*radfor
      fory = yr*radfor
      forz = zr*radfor

      IF (idielec /= 0) THEN

        xr = 2D0*xdielec - rvectmp(1) - cx(jj)
        yr = rvectmp(2) - cy(jj)
        zr = rvectmp(3) - cz(jj)

        dist = xr*xr + yr*yr + zr*zr
        dist = SQRT(dist)

        smix1 = sigsig(sigtmp, sgm1(np(jj)))
        smix2 = sigsig(sigtmp, sgm2(np(jj)))

        radfor = -(epsdi - 1)/(epsdi + 1)*e2*chgtmp*(chg1(np(jj))* &
                                                     dv_softdr(dist, smix1*sq2) + chg2(np(jj))*dv_softdr(dist, smix2*sq2)) &
                 /dist

        forx = xr*radfor + forx
        fory = yr*radfor + fory
        forz = zr*radfor + forz

      END IF

      IF (iflag == 0) THEN
        fx(jj) = fx(jj) + forx
        fy(jj) = fy(jj) + fory
        fz(jj) = fz(jj) + forz
      END IF

      IF (idielec /= 0) THEN !because force GSMimage-ion is not
! -(force ionimage-GSM)
        IF (iflag == 0) THEN
          forx = forx - 2D0*xr*radfor
! fory = fory-2D0*yr*radfor
! forz = forz-2D0*zr*radfor
        END IF
      END IF

      CALL addforce(ii, -forx, -fory, -forz)

      ii2 = iconvlongtoshort(ii)

      CALL getshortforce(iptyp(ii), 4, ii2, jj, dummy, 0)

    END DO

  END DO

  RETURN
END SUBROUTINE getforcenagsm
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getforcegsmgsm(iflag)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iflag

  INTEGER :: ii, ii2, jj, jj2
  REAL(DP) :: cml, chi, dist, dist2, radfor, si, smix, sigsig, xi, xr, yi, yr, zi, zr
  REAL(DP) :: dummy(2*kdfull2)

  INTEGER, EXTERNAL :: iconvlongtoshort, ismobile, iptyp
  REAL(DP), EXTERNAL :: dv_softdr

  cml = 0.541382_DP

  IF (ibh == 0) THEN

    DO ii = 1, nc + NE + nk - 1

      CALL getparas(ii)

      xi = rvectmp(1)
      yi = rvectmp(2)
      zi = rvectmp(3)
      si = sigtmp
      chi = chgtmp

      DO jj = ii + 1, nc + NE + nk

        IF (iflag == 0 .OR. iptyp(ii) == 2 .OR. iptyp(jj) == 2) THEN

          IF (ismobile(ii) /= 0 .OR. ismobile(jj) /= 0) THEN

! IF ONLY forces on valences are TO be calculated, skip
! all other combinations of particle sorts

            CALL getparas(jj)

            xr = xi - rvectmp(1)
            yr = yi - rvectmp(2)
            zr = zi - rvectmp(3)

            dist2 = xr*xr + yr*yr + zr*zr
            dist = SQRT(dist2)

            smix = sigsig(si, sigtmp)

            IF (jj <= nc + NE .AND. ii + nc == jj) THEN
              IF (iflag /= 2) THEN
                radfor = cspr
              ELSE
                radfor = 0D0
              END IF

            ELSE

              IF (dist < smix*8*sq2) THEN
! cut-off for error FUNCTION: dist < smix * 8
! means: rel. error is smaller than 1e-18

                radfor = e2*chgtmp*chi*dv_softdr(dist, smix*sq2)/dist
! radfor = e2*chgTmp*chi*dVsdr(dist,smix*SQ2)/dist

              ELSE
                radfor = -e2*chgtmp*chi/(dist2*dist) ! minus is from derivative
              END IF

            END IF

            forx = xr*radfor
            fory = yr*radfor
            forz = zr*radfor

! IF (ii.le.nc .and. jj.gt.ne+nc) THEN
! WRITE(6,*) forx,fory,forz,iptyp(ii),iptyp(jj),ii,jj
! END IF

            IF (iflag == 0) THEN ! as usual
              CALL addforce(ii, -forx, -fory, -forz)
              CALL addforce(jj, forx, fory, forz)
            ELSE ! ONLY store valence forces
              IF (iptyp(ii) == 2) CALL addforce(ii, -forx, -fory, -forz)
              IF (iptyp(jj) == 2) CALL addforce(jj, forx, fory, forz)
            END IF

            ii2 = iconvlongtoshort(ii)
            jj2 = iconvlongtoshort(jj)

!test WRITE(*,*) '1: ii2,jj2 = ', ii2,jj2,ii,jj

            CALL getshortforce(iptyp(ii), iptyp(jj), ii2, jj2, dummy, 0)

          END IF ! ismobile

        END IF

      END DO

    END DO

! DO i=1,nc
! WRITE(6,*) xcinit(i),ycinit(i),zcinit(i)
! END DO

  ELSE ! ibh = 1

! first round with actual positions

    DO ii = 1, nc + NE + nk

      CALL getparas(ii)

      xi = rvectmp(1)
      yi = rvectmp(2)
      zi = rvectmp(3)
      si = sigtmp
      chi = chgtmp

      IF (imobtmp /= 0) THEN

        DO jj = 1, nc + NE + nk

          CALL getparas(jj)

          IF ((imobtmp /= 0 .OR. iptyp(ii) == 2) .AND. ii /= jj) THEN

            xr = xi - rvectmp(1)
            yr = yi - rvectmp(2)
            zr = zi - rvectmp(3)

            dist2 = xr*xr + yr*yr + zr*zr
            dist = SQRT(dist2)

            smix = sigsig(si, sigtmp)

            IF (MAX(ii, jj) <= nc + NE .AND. ABS(ii - jj) == nc) THEN
              IF (iflag /= 2) THEN
                radfor = cspr
              ELSE
                radfor = 0D0
              END IF

            ELSE

              IF (dist < smix*8*sq2) THEN
                radfor = e2*chgtmp*chi*dv_softdr(dist, smix*sq2)/dist
              ELSE
                radfor = -e2*chgtmp*chi/(dist2*dist) ! minus is from derivative
              END IF

            END IF

            forx = xr*radfor
            fory = yr*radfor
            forz = zr*radfor

            IF (iflag == 0) THEN ! as usual
              CALL addforce(ii, -forx, -fory, -forz)
            ELSE ! ONLY store valence forces
              IF (iptyp(ii) == 2) CALL addforce(ii, -forx, -fory, -forz)
            END IF

            ii2 = iconvlongtoshort(ii)
            jj2 = iconvlongtoshort(jj)

!test WRITE(*,*) '2: ii2,jj2 = ', ii2,jj2,ii,jj

            CALL getshortforce(iptyp(ii), iptyp(jj), ii2, jj2, dummy, 1)

          END IF

        END DO

      END IF

    END DO

    DO jj = 1, nc
      xcold(jj) = xc(jj)
      ycold(jj) = yc(jj)
      zcold(jj) = zc(jj)
    END DO
    DO jj = 1, NE
      xeold(jj) = xe(jj)
      yeold(jj) = ye(jj)
      zeold(jj) = ze(jj)
    END DO
    DO jj = 1, nk
      xkold(jj) = xk(jj)
      ykold(jj) = yk(jj)
      zkold(jj) = zk(jj)
    END DO

! goto 777

! second round with initial positions

    DO ii = 1, nc + NE + nk

      IF (iptyp(ii) /= 2) THEN

        CALL getparas(ii)

        xi = rvectmp(1)
        yi = rvectmp(2)
        zi = rvectmp(3)
        si = sigtmp
        chi = chgtmp

        DO jj = 1, nc
          IF (jj /= ii) THEN
            xc(jj) = xcinit(jj)
            yc(jj) = ycinit(jj)
            zc(jj) = zcinit(jj)
          END IF
        END DO
        DO jj = 1, NE
          IF (jj + nc /= ii) THEN
            xe(jj) = xeinit(jj)
            ye(jj) = yeinit(jj)
            ze(jj) = zeinit(jj)
          END IF
        END DO
        DO jj = 1, nk
          IF (jj + nc + NE /= ii) THEN
            xk(jj) = xkinit(jj)
            yk(jj) = ykinit(jj)
            zk(jj) = zkinit(jj)
          END IF
        END DO

        IF (imobtmp /= 0) THEN

          DO jj = 1, nc + NE + nk

            CALL getparas(jj)

            IF (imobtmp /= 0 .AND. ii /= jj) THEN

              xr = xi - rvectmp(1)
              yr = yi - rvectmp(2)
              zr = zi - rvectmp(3)

              dist2 = xr*xr + yr*yr + zr*zr
              dist = SQRT(dist2)

              smix = sigsig(si, sigtmp)

              IF (MAX(ii, jj) <= nc + NE .AND. ABS(ii - jj) == nc) THEN
                IF (iflag /= 2) THEN
                  radfor = cspr
                ELSE
                  radfor = 0D0
                END IF

              ELSE

                IF (dist < smix*8D0*sq2) THEN
                  radfor = e2*chgtmp*chi*dv_softdr(dist, smix*sq2)/dist
                ELSE
                  radfor = -e2*chgtmp*chi/(dist2*dist) ! minus is from derivative
                END IF

              END IF

              forx = xr*radfor
              fory = yr*radfor
              forz = zr*radfor

              IF (iflag == 0) THEN ! as usual
                CALL addforce(ii, forx, fory, forz)
              ELSE ! ONLY store valence forces
                IF (iptyp(ii) == 2) CALL addforce(ii, -forx, -fory, -forz)
              END IF

              ii2 = iconvlongtoshort(ii)
              jj2 = iconvlongtoshort(jj)

              CALL getshortforce(iptyp(ii), iptyp(jj), ii2, jj2, dummy, 2)
! CALL getShortForce(iptyp(jj),iptyp(ii),jj2,ii2,rho,1)

            END IF

          END DO ! jj-loop

        END IF ! imobTmp(i)=0

        DO jj = 1, nc
          IF (jj /= ii) THEN
            xc(jj) = xcold(jj)
            yc(jj) = ycold(jj)
            zc(jj) = zcold(jj)
          END IF
        END DO
        DO jj = 1, NE
          IF (jj + nc /= ii) THEN
            xe(jj) = xeold(jj)
            ye(jj) = yeold(jj)
            ze(jj) = zeold(jj)
          END IF
        END DO
        DO jj = 1, nk
          IF (jj + nc + NE /= ii) THEN
            xk(jj) = xkold(jj)
            yk(jj) = ykold(jj)
            zk(jj) = zkold(jj)
          END IF
        END DO

      END IF

    END DO

! third round: diagonal terms

    DO ii = 1, nc
      forx = -cml*(xc(ii) - xcinit(ii))
      fory = -cml*(yc(ii) - ycinit(ii))
      forz = -cml*(zc(ii) - zcinit(ii))
      CALL addforce(ii, forx, fory, forz)
    END DO
! DO ii=1,ne
! forx = -cml*(xe(ii)-xeinit(ii))
! fory = -cml*(ye(ii)-yeinit(ii))
! forz = -cml*(ze(ii)-zeinit(ii))
! CALL addforce(ii+nc,forx,fory,forz)
! END DO
    DO ii = 1, nk
      forx = -cml*(xk(ii) - xkinit(ii))
      fory = -cml*(yk(ii) - ykinit(ii))
      forz = -cml*(zk(ii) - zkinit(ii))
      CALL addforce(ii + NE + nc, forx, fory, forz)
    END DO

!$$$ WRITE(6,*) '-------------------------------------------------'
!$$$ dmf2c = 0.
!$$$ dmf2e = 0.
!$$$ dmf2k = 0.
!$$$ DO ii=1,nc
!$$$ IF (imobc(ii).ne.0) THEN
!$$$ WRITE(6,'(a,i,3e10.2)') 'c:',ii,fxc(ii),fyc(ii),fzc(ii)
!$$$ IF (fxc(ii)**2+fyc(ii)**2+fzc(ii)**2.gt.dmf2c) dmf2c =
!$$$ & fxc(ii)**2+fyc(ii)**2+fzc(ii)**2
!$$$ END IF
!$$$ END DO
!$$$ DO ii=1,ne
!$$$ IF (imobe(ii).ne.0) THEN
!$$$ WRITE(6,'(a,i,3e10.2)') 'e:',ii,fxe(ii),fye(ii),fze(ii)
!$$$ IF (fxe(ii)**2+fye(ii)**2+fze(ii)**2.gt.dmf2e) dmf2e =
!$$$ & fxe(ii)**2+fye(ii)**2+fze(ii)**2
!$$$ END IF
!$$$ END DO
!$$$ DO ii=1,nk
!$$$ IF (imobk(ii).ne.0) THEN
!$$$ WRITE(6,'(a,i,3e10.2)') 'k:',ii,fxk(ii),fyk(ii),fzk(ii)
!$$$ IF (fxk(ii)**2+fyk(ii)**2+fzk(ii)**2.gt.dmf2k) dmf2k =
!$$$ & fxk(ii)**2+fyk(ii)**2+fzk(ii)**2
!$$$
!$$$ END IF
!$$$
!$$$ END DO
!$$$c STOP
!$$$ WRITE(6,*) dmf2c**0.5,dmf2e**0.5,dmf2k**0.5
!$$$ WRITE(6,*) '-------------------------------------------------'

  END IF !ibh

!$$$ IF (ibh.eq.1) THEN
!$$$
!$$$ CALL cpArr(xc,xcold,nc)
!$$$ CALL cpArr(yc,ycold,nc)
!$$$ CALL cpArr(zc,zcold,nc)
!$$$
!$$$ CALL cpArr(xe,xeold,ne)
!$$$ CALL cpArr(ye,yeold,ne)
!$$$ CALL cpArr(ze,zeold,ne)
!$$$
!$$$ CALL cpArr(xk,xkold,nk)
!$$$ CALL cpArr(yk,ykold,nk)
!$$$ CALL cpArr(zk,zkold,nk)
!$$$
!$$$
!$$$ DO ii=1,nc+ne+nk
!$$$
!$$$c TO be commented OUT IN other cases:
!$$$ IF (ii.le.nc .or. (ii.gt.nc+ne)) THEN
!$$$
!$$$
!$$$ CALL getParas(ii)
!$$$
!$$$ xi = rVecTmp(1)
!$$$ yi = rVecTmp(2)
!$$$ zi = rVecTmp(3)
!$$$ si = sigTmp
!$$$ chi = chgTmp
!$$$
!$$$
!$$$ DO jj=1,nc+ne+nk
!$$$
!$$$ IF (ii.ne.jj) THEN
!$$$
!$$$ IF (jj.le.nc) THEN
!$$$ xc(jj)=xcinit(jj)
!$$$ yc(jj)=ycinit(jj)
!$$$ zc(jj)=zcinit(jj)
!$$$ elseif (jj.le.nc+ne) THEN
!$$$ xe(jj)=xeinit(jj)
!$$$ ye(jj)=yeinit(jj)
!$$$ ze(jj)=zeinit(jj)
!$$$ ELSE
!$$$ xk(jj)=xkinit(jj)
!$$$ yk(jj)=ykinit(jj)
!$$$ zk(jj)=zkinit(jj)
!$$$ END IF
!$$$
!$$$
!$$$ IF (iflag.eq.0 .or. iptyp(ii).eq.2 .or.iptyp(jj).eq.2) THEN
!$$$
!$$$ IF (isMobile(ii).ne.0 .or. isMobile(jj).ne.0) THEN
!$$$
!$$$
!$$$ ! IF ONLY forces on valences are TO be calculated, skip
!$$$ ! all other combinations of particle sorts
!$$$
!$$$
!$$$ CALL getParas(jj)
!$$$
!$$$
!$$$ xr = xi - rVecTmp(1)
!$$$ yr = yi - rVecTmp(2)
!$$$ zr = zi - rVecTmp(3)
!$$$
!$$$
!$$$ dist2 = xr*xr+yr*yr+zr*zr
!$$$ dist = sqrt(dist2)
!$$$
!$$$c$$$ IF (dist2.le.1e-8 .and.abs(ii-jj).ne.nc)
!$$$c$$$ & WRITE(6,*) dist2
!$$$
!$$$ smix = sigsig(si,sigTmp)
!$$$
!$$$ IF (ii.le.nc+ne .and.
!$$$ & jj.le.nc+ne .and. abs(ii-jj).eq.nc) THEN
!$$$
!$$$ radfor = cspr
!$$$
!$$$ ELSE
!$$$
!$$$ IF (dist.lt.smix*8*SQ2) THEN
!$$$ ! cut-off for error FUNCTION: dist < smix * 8
!$$$ ! means: rel. error is smaller than 1e-18
!$$$
!$$$ radfor = e2*chgTmp*chi*DV_softDr(dist,smix*SQ2)/dist
!$$$c radfor = e2*chgTmp*chi*dVsdr(dist,smix*SQ2)/dist
!$$$
!$$$ ELSE
!$$$ radfor = -e2*chgTmp*chi/(dist2*dist) ! minus is from derivative
!$$$ END IF
!$$$
!$$$ END IF
!$$$
!$$$ forx = xr*radfor
!$$$ fory = yr*radfor
!$$$ forz = zr*radfor
!$$$
!$$$c IF (ii.le.nc .and. jj.gt.ne+nc) THEN
!$$$c WRITE(6,*) forx,fory,forz,iptyp(ii),iptyp(jj),ii,jj
!$$$c END IF
!$$$
!$$$
!$$$
!$$$ IF (iflag.eq.0) THEN ! as usual
!$$$ CALL addForce(ii,forx,fory,forz)
!$$$c CALL addForce(jj,-forx,-fory,-forz)
!$$$ ELSE ! ONLY store valence forces
!$$$ IF (iptyp(ii).eq.2)CALL addForce(ii,forx,fory,forz)
!$$$c IF (iptyp(jj).eq.2)CALL addForce(jj,-forx,-fory,-forz)
!$$$ END IF
!$$$
!$$$ ii2 = iconvLongToShort(ii)
!$$$ jj2 = iconvLongToShort(jj)
!$$$
!$$$c WRITE(6,*) 'ii2,jj2 = ', ii2,jj2,ii,jj
!$$$
!$$$
!$$$ CALL getShortForce(iptyp(ii),iptyp(jj),ii2,jj2,rho,iflag,1)
!$$$ ! since rho is not needed, it is a dummy field
!$$$
!$$$
!$$$ END IF ! ismobile
!$$$
!$$$ END IF
!$$$
!$$$ IF (jj.le.nc) THEN
!$$$ xc(jj)=xcold(jj)
!$$$ yc(jj)=ycold(jj)
!$$$ zc(jj)=zcold(jj)
!$$$ elseif (jj.le.nc+ne) THEN
!$$$ xe(jj)=xeold(jj)
!$$$ ye(jj)=yeold(jj)
!$$$ ze(jj)=zeold(jj)
!$$$ ELSE
!$$$ xk(jj)=xkold(jj)
!$$$ yk(jj)=ykold(jj)
!$$$ zk(jj)=zkold(jj)
!$$$ END IF
!$$$
!$$$
!$$$ END IF
!$$$
!$$$ END DO
!$$$
!$$$ END IF ! skip valences
!$$$
!$$$ END DO
!$$$
!$$$c END IF
!$$$
!$$$ DO ii=1,nc
!$$$ fxc(ii)=fxc(ii)-cbh*(xc(ii)-xcinit(ii))
!$$$ fyc(ii)=fyc(ii)-cbh*(yc(ii)-ycinit(ii))
!$$$ fzc(ii)=fzc(ii)-cbh*(zc(ii)-zcinit(ii))
!$$$ END DO
!$$$
!$$$c$$$ DO ii=1,ne
!$$$c$$$ fxe(ii)=fxe(ii)-cbh*(xe(ii)-xeinit(ii))
!$$$c$$$ fye(ii)=fye(ii)-cbh*(ye(ii)-yeinit(ii))
!$$$c$$$ fze(ii)=fze(ii)-cbh*(ze(ii)-zeinit(ii))
!$$$c$$$ END DO
!$$$
!$$$ DO ii=1,nk
!$$$ fxk(ii)=fxk(ii)-cbh*(xk(ii)-xkinit(ii))
!$$$ fyk(ii)=fyk(ii)-cbh*(yk(ii)-ykinit(ii))
!$$$ fzk(ii)=fzk(ii)-cbh*(zk(ii)-zkinit(ii))
!$$$ END DO
!$$$
!$$$
!$$$
!$$$
!$$$ END IF

  IF (idielec /= 0) THEN

    DO ii = 1, nc + NE + nk

      CALL getparas(ii)

      xi = rvectmp(1)
      yi = rvectmp(2)
      zi = rvectmp(3)
      si = sigtmp
      chi = chgtmp

      DO jj = 1, nc + NE + nk

        CALL getparas(jj)

        xr = xi + rvectmp(1) - 2D0*xdielec
        yr = yi - rvectmp(2)
        zr = zi - rvectmp(3)

        dist2 = xr*xr + yr*yr + zr*zr
        dist = SQRT(dist2)

        IF (dist < sigtmp*8*sq2) THEN
          radfor = -(epsdi - 1)/(epsdi + 1)* &
                   e2*chgtmp*chi*dv_softdr(dist, sigtmp*sq2)/dist
        ELSE
          radfor = (epsdi - 1)/(epsdi + 1)*e2*chgtmp*chi/(dist2*dist)
        END IF

        forx = xr*radfor
        fory = yr*radfor
        forz = zr*radfor

        CALL addforce(ii, -forx, -fory, -forz)

      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE getforcegsmgsm
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getvdwforce(rho)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! THIS ONLY WORKS FOR THE ARGON CASE!!!

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)

  INTEGER :: i, ico, ii, ind, isign, ix, iy, iz
  REAL(DP) :: chpvdw, ds, r2, rr, x1, y1, z1
! DIMENSION frhop(ngpar,3),frhom(ngpar,3)
  REAL(DP) :: ri(3)
  REAL(DP) :: xcsave(ngpar), ycsave(ngpar), zcsave(ngpar)
  REAL(DP) :: f(ngpar, 3)

  REAL(DP), EXTERNAL :: v_vdw

  DO i = 1, nc
    xcsave(i) = xc(i)
    ycsave(i) = yc(i)
    zcsave(i) = zc(i)
  END DO

  DO i = 1, nc
    DO ico = 1, 3
      f(i, ico) = 0D0
    END DO
  END DO
  ds = 1D-5

  DO ico = 1, 3

    DO ISIGN = 1, -1, -2

      IF (ISIGN == 1) THEN
        CALL shiftpos(1, ico, ds)
      ELSE
        CALL shiftpos(1, ico, -ds - ds)
      END IF

      CALL calc_frho(rho)

      DO ii = 1, nc

        ind = 0
        DO iz = minz, maxz
          z1 = (iz - nzsh)*dz
          ri(3) = zc(ii) - z1

          DO iy = miny, maxy
            y1 = (iy - nysh)*dy
            ri(2) = yc(ii) - y1

            DO ix = minx, maxx
              x1 = (ix - nxsh)*dx
              ri(1) = xc(ii) - x1

              ind = ind + 1

              r2 = ri(1)*ri(1) + ri(2)*ri(2) + ri(3)*ri(3)
              rr = MAX(SQRT(r2), small)

              chpvdw = rho(ind)*v_vdw(ri, rr, ii, 0.5D0)/rr

! f(ii,ico) = f(ii,ico) - isign*chpvdw*ri(ico)
              f(ii, ico) = f(ii, ico) + ISIGN*chpvdw*ri(ico)

            END DO
          END DO
        END DO

      END DO ! nc

    END DO ! sign

    DO ii = 1, nc
      xc(ii) = xcsave(ii)
      yc(ii) = ycsave(ii)
      zc(ii) = zcsave(ii)
    END DO

  END DO ! coordinates

  DO ii = 1, nc
    fxc(ii) = fxc(ii) + f(ii, 1)*dvol/(ds + ds)
    fyc(ii) = fyc(ii) + f(ii, 2)*dvol/(ds + ds)
    fzc(ii) = fzc(ii) + f(ii, 3)*dvol/(ds + ds)
  END DO

  CALL calc_frho(rho)

  RETURN
END SUBROUTINE getvdwforce
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getforcegsmgsm1(iflag)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iflag

  INTEGER :: ii, ii2, imi, iw, jj, jj2, jw
  REAL(DP) :: dist, dist2, chi, radfor, si, sigsig, smix, xi, yi, zi, xr, yr, zr
  REAL(DP) :: dummy(2*kdfull2)

  INTEGER, EXTERNAL :: iconvlongtoshort, iptyp
  REAL(DP), EXTERNAL :: dv_softdr

  DO ii = 1, nc + NE + nk ! mobile and fixed

    iw = ipointer(ii)
    CALL getparas(iw)

    xi = rvectmp(1)
    yi = rvectmp(2)
    zi = rvectmp(3)
    si = sigtmp
    chi = chgtmp
    imi = imobtmp

    DO jj = 1, MIN(ii - 1, nmob) ! mobile

      jw = ipointer(jj)
      CALL getparas(jw)

      xr = xi - rvectmp(1)
      yr = yi - rvectmp(2)
      zr = zi - rvectmp(3)

      dist2 = xr*xr + yr*yr + zr*zr
      dist = SQRT(dist2)

      smix = sigsig(si, sigtmp)

!spring force or normal?

      IF (iptyp(jw) == 2 .AND. jw - nc == iw) THEN
        IF (iflag /= 2) THEN
          radfor = cspr
        ELSE
          radfor = 0D0
        END IF
      ELSE IF (iptyp(iw) == 2 .AND. iw - nc == jw) THEN
        IF (iflag /= 2) THEN
          radfor = cspr
        ELSE
          radfor = 0.0D0
        END IF

      ELSE

        IF (dist < smix*8*sq2) THEN
! TO SPEED UP TURN DOWN CUTOFF!
! cut-off for error FUNCTION: dist < smix * 8
! means: rel. error is smaller than 1e-18
          radfor = e2*chgtmp*chi*dv_softdr(dist, smix*sq2)/dist
        ELSE
          radfor = -e2*chgtmp*chi/(dist2*dist) ! minus is from derivative
        END IF

      END IF

      forx = xr*radfor
      fory = yr*radfor
      forz = zr*radfor

      IF (iflag == 0) THEN ! as usual
        CALL addforce(iw, -forx, -fory, -forz)
        CALL addforce(jw, forx, fory, forz)
      ELSE ! ONLY store valence forces
        IF (iptyp(iw) == 2) CALL addforce(iw, -forx, -fory, -forz)
        IF (iptyp(jw) == 2) CALL addforce(jw, forx, fory, forz)
      END IF

      ii2 = iconvlongtoshort(iw)
      jj2 = iconvlongtoshort(jw)

      CALL getshortforce(iptyp(iw), iptyp(jw), ii2, jj2, dummy, 0)

    END DO

  END DO

  RETURN
END SUBROUTINE getforcegsmgsm1
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE getforces_clust2cores(rho, iflag)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! GSM means
! a particle described by the Gaussian Shell Model, i.e.
! cores, valence shells and cations.

! - calculation of the forces between GSM particles and
! cluster cores, labelled Na
! - calculation of the force on GSM and cluster cores due
! TO the electronic density

! iflag ---> 0 calculate forces on all particles
! 1 forces on valence shells ONLY
! 2 forces on valence shells ONLY, but
! skipping the spring force
! (used IN adjustdip)

  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: iflag

  INTEGER :: ii, idielecsave
  LOGICAL :: ipseudosave

  IF (iflag == 0) THEN
    DO ii = 1, nion
      fx(ii) = 0.0D0
      fy(ii) = 0.0D0
      fz(ii) = 0.0D0
    END DO
    DO ii = 1, nc
      fxc(ii) = 0.0D0
      fyc(ii) = 0.0D0
      fzc(ii) = 0.0D0
    END DO
    DO ii = 1, nk
      fxk(ii) = 0.0D0
      fyk(ii) = 0.0D0
      fzk(ii) = 0.0D0
    END DO
  END IF

  DO ii = 1, NE
    fxe(ii) = 0.0D0
    fye(ii) = 0.0D0
    fze(ii) = 0.0D0
  END DO

! GSM relevant parts

  idielecsave = idielec
  ipseudosave = ipseudo
  idielec = 0
  ipseudo = .FALSE.

  CALL getforcenagsm(iflag)

  IF (nelect > 0) CALL getforceelgsm(iflag, rho)

  idielec = idielecsave
  ipseudo = ipseudosave

! van der Waals part
  IF (ivdw == 1) THEN
    CALL getvdwforce(rho)
  END IF ! IF ivdw.eq.2 vdw is done implicitly by vArElCore

  RETURN
END SUBROUTINE getforces_clust2cores
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE shiftpos(ityp, ico, ds)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ityp ! SELECT : cores / cation / shells /ions clusters
  ! 1 2 3 4
  INTEGER, INTENT(IN) :: ico ! SELECT coordinate : x / y / z
  ! 1 2 3
  REAL(DP), INTENT(IN) :: ds

  INTEGER :: i

  SELECT CASE (ityp)
!--------
  CASE (1)
! Cores
    SELECT CASE (ico)
    CASE (1) ! x
      DO i = 1, nc
        xc(i) = xc(i) + ds
      END DO
    CASE (2) ! y
      DO i = 1, nc
        yc(i) = yc(i) + ds
      END DO
    CASE (3) ! z
      DO i = 1, nc
        zc(i) = zc(i) + ds
      END DO
    CASE DEFAULT
      STOP 'Error IN shiftPos'
    END SELECT
!--------
  CASE (2)
! Cations
    SELECT CASE (ico)
    CASE (1) ! x
      DO i = 1, NE
        xe(i) = xe(i) + ds
      END DO
    CASE (2) ! y
      DO i = 1, NE
        ye(i) = ye(i) + ds
      END DO
    CASE (3) !z
      DO i = 1, NE
        ze(i) = ze(i) + ds
      END DO
    CASE DEFAULT
      STOP 'Error IN shiftPos'
    END SELECT
!--------
  CASE (3)
! Shells
    SELECT CASE (ico)
    CASE (1) ! x
      DO i = 1, nk
        xk(i) = xk(i) + ds
      END DO
    CASE (2) ! y
      DO i = 1, nk
        yk(i) = yk(i) + ds
      END DO
    CASE (3) !z
      DO i = 1, nk
        zk(i) = zk(i) + ds
      END DO
    CASE DEFAULT
      STOP 'Error IN shiftPos'
    END SELECT
!--------
  CASE (4)
! Ions clusters
    SELECT CASE (ico)
    CASE (1) ! x
      DO i = 1, nion
        cx(i) = cx(i) + ds
      END DO
    CASE (2) ! y
      DO i = 1, nion
        cy(i) = cy(i) + ds
      END DO
    CASE (3) ! z
      DO i = 1, nion
        cz(i) = cz(i) + ds
      END DO
    CASE DEFAULT
      STOP 'Error IN shiftPos'
    END SELECT
!--------
  CASE DEFAULT
    STOP 'Error IN shiftPos'
  END SELECT

  RETURN
END SUBROUTINE shiftpos
!------------------------------------------------------------

SUBROUTINE madelung
!------------------------------------------------------------
! calculates force on particle due TO electrostatic
! EXTERNAL Madelung potential (or pseudoparticles)
  USE params
  IMPLICIT NONE

  INTEGER :: i
  REAL(DP) :: prefac, xi, yi, zi

  REAL(DP), EXTERNAL :: v_soft, gauss

! WRITE(6,*) 'before: ', fx(1),fy(1),fz(1)

  IF (.NOT. ipseudo) THEN

    DO i = 1, nion
      CALL foldgradfunc(potstat, gauss, cx(i), cy(i), cz(i), sgm1(np(i))*sq2)

      fx(i) = fx(i) - e2*chg1(np(i))*rvectmp(1)
      fy(i) = fy(i) - e2*chg1(np(i))*rvectmp(2)
      fz(i) = fz(i) - e2*chg1(np(i))*rvectmp(3)

      CALL foldgradfunc(potstat, gauss, cx(i), cy(i), cz(i), sgm2(np(i))*sq2)

      fx(i) = fx(i) - e2*chg2(np(i))*rvectmp(1)
      fy(i) = fy(i) - e2*chg2(np(i))*rvectmp(2)
      fz(i) = fz(i) - e2*chg2(np(i))*rvectmp(3)

    END DO

    DO i = 1, nc + NE + nk
      CALL getparas(i)

      xi = rvectmp(1)
      yi = rvectmp(2)
      zi = rvectmp(3)

      CALL foldgradfunc(potstat, gauss, xi, yi, zi, sigtmp*sq2)

      CALL addforce(i, -e2*chgtmp*rvectmp(1), -e2*chgtmp*rvectmp(2), &
                    -e2*chgtmp*rvectmp(3))

    END DO

  ELSE ! ipseudo = .TRUE.

    DO i = 1, nion

      CALL foldgradfunconsubgrid(potstat, gauss, cx(i), cy(i), cz(i), &
                                 sgm1(np(1))*sq2)

      prefac = chg1(np(i))/(2*pi*sgm1(np(i))**2)**1.5D0
! check IF factor e2 is already IN phim!! yes it is!!

      fx(i) = fx(i) - prefac*rvectmp(1)
      fy(i) = fy(i) - prefac*rvectmp(2)
      fz(i) = fz(i) - prefac*rvectmp(3)

      CALL foldgradfunconsubgrid(potstat, gauss, cx(i), cy(i), cz(i), &
                                 sgm2(np(1))*sq2)

      prefac = chg2(np(i))/(2D0*pi*sgm2(np(i))**2)**1.5D0
! check IF factor e2 is already IN phim!! yes it is!!

      fx(i) = fx(i) - prefac*rvectmp(1)
      fy(i) = fy(i) - prefac*rvectmp(2)
      fz(i) = fz(i) - prefac*rvectmp(3)

    END DO

! WRITE(6,*) 'fx1',fzc(1)

    DO i = 1, nc + NE + nk
      CALL getparas(i)

      xi = rvectmp(1)
      yi = rvectmp(2)
      zi = rvectmp(3)

      CALL foldgradfunconsubgrid(potstat, gauss, xi, yi, zi, sigtmp*sq2)

      CALL addforce(i, -e2*chgtmp*rvectmp(1), -e2*chgtmp*rvectmp(2), &
                    -e2*chgtmp*rvectmp(3))

    END DO
! WRITE(6,*) 'fx2',fzc(1)

  END IF

! WRITE(6,*) 'after ', fx(1),fy(1),fz(1)
! STOP

END SUBROUTINE madelung
!------------------------------------------------------------

SUBROUTINE adjustdip(rho, it)
! computes the force on the GSM-dipoles and adjusts their dipole moments
! IN static iteration (with jdip=1 IN dynamic CASE)
!----------------------------------------------------------------------------
  USE params
  USE util, ONLY: prifld
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: it

  COMPLEX(DP) :: psidummy(1)
  REAL(DP), PARAMETER :: delmrmax = 2D-6
  INTEGER, PARAMETER :: maxadjust = 200
  LOGICAL, PARAMETER :: ttest = .false.

  INTEGER :: i, icoun, iit, irar, iter
  REAL(DP) :: delmr, delmrold, delt
  REAL(DP) :: dmspox, dmspoy, dmspoz, dmsfx, dmsfy, dmsfz, dmsfxold, dmsfyold, dmsfzold
  REAL(DP) :: posox, posoy, posoz, rmaxpol, rrrr
  REAL(DP) :: delmr_iter(maxadjust)
  REAL(DP), ALLOCATABLE :: xm(:) ! array for actual particle mass
  LOGICAL::converged = .false. ! will be true IF the loop exits before reaching maxadjust

  REAL(DP), EXTERNAL :: rgetmaxpol, rgetmeanzpol

  IF (NE == 0) RETURN

  IF (ttest) WRITE (6, *) 'Entering adjustdip. rhonorm=', SUM(rho(1:kdfull2))*dvol
  IF (ttest) CALL prifld(rho, 'density adj')

! the 'psidummy' trick works ONLY for local Gaussian PsP
  IF (ipsptyp == 1) STOP ' ADJUSTDIP must not be used with Goedecker PsP'

  delt = 0.35D0 !0.5D0

! SAVE positions for the CASE that the iteration does not
! converge

  DO i = 1, NE
    xeold(i) = xe(i)
    yeold(i) = ye(i)
    zeold(i) = ze(i)
  END DO

! goto 15

  delmr = 1.0D0
  delmrold = delmr

  DO iter = 1, maxadjust
    IF (SQRT(delmr) < delmrmax) THEN
      converged = .true.
      EXIT
    END IF
    delmr = 0D0
! compute forces of valence clouds
! CALL calcforce_dip(rho)

    CALL getforces(rho, psidummy, it, 2)
    icoun = 0
    rmaxpol = rgetmaxpol()
    rrrr = rgetmeanzpol()

    IF (ttest) WRITE (6, '(a,2e17.7)') ' maxpol, meanpol = ', rMaxPol, rrrr

! CALL rgetMaxForce(2)

! WRITE(6,'(a,3e17.7)') 'max. fx, fy, fz: ', rVecTmp(1),
! & rVecTmp(2), rVecTmp(3)

    DO irar = 1, NE
      posox = xe(irar)
      posoy = ye(irar)
      posoz = ze(irar)

      IF (ttest) WRITE (*, *) irar, xe(irar) - xc(irar), ye(irar) - yc(irar), &
        ze(irar) - zc(irar), fxe(irar), fye(irar), fze(irar)

      IF (imobe(irar) /= 0) THEN
        xe(irar) = xe(irar) - delt*(xe(irar) - xc(irar) - fxe(irar)/cspr) !(e2*C_dipmod)
        ye(irar) = ye(irar) - delt*(ye(irar) - yc(irar) - fye(irar)/cspr) !idem
        ze(irar) = ze(irar) - delt*(ze(irar) - zc(irar) - fze(irar)/cspr) !idem

!$$$ xe(irar)=
!$$$ & xc(irar)+fxe(irar)/cspr !(e2*C_dipmod)
!$$$ ye(irar)=
!$$$ & yc(irar)+fye(irar)/cspr !idem
!$$$ ze(irar)=
!$$$ & zc(irar)+fze(irar)/cspr !idem

        icoun = icoun + 1

        delmr = delmr + (posox - xe(irar))**2 + (posoy - ye(irar))**2 &
                + (posoz - ze(irar))**2

      END IF

      IF (ipseudo) THEN
        CALL updatesubgrids
      END IF

    END DO
! delmr=delmr/ne
    delmr = delmr/icoun
    delmr_iter(iter) = delmr
    IF (ttest) WRITE (6, *) ' it,delmr = ', iter, sqrt(delmr)

! IF (delmr.gt.delmrold) STOP 'adjustdip not converging'
! IF (delmr.gt.delmrold) goto 15

    delmrold = delmr

  END DO

  IF ((delmr .GT. delmrmax*10D0) .AND. .NOT. converged) THEN
    WRITE (*, '(a)') ' Iteration DELMR:'
    WRITE (*, '(i5,1pg13.5)') (iter, delmr_iter(iter), iter=1, maxadjust)
    STOP 'routine adjustdip did not converge'
  END IF

  IF (ifadiadip == 1) THEN
    DO irar = 1, NE
      pxe(irar) = 0D0
      pye(irar) = 0D0
      pze(irar) = 0D0
    END DO
  END IF

! CALL pripolariz

! WRITE(6,*) 'Leaving adjustdip'
  RETURN

! normal iteration did not converge; make cooling dynamics instead

15 WRITE (6, *) 'adjustdip did not converge. Making cooling dynamics.'

  DO i = 1, NE
    xe(i) = xeold(i)
    ye(i) = yeold(i)
    ze(i) = zeold(i)
  END DO

  DO iit = 1, 200

    IF (iit > 1) THEN
      dmsfxold = dmsfx
      dmsfyold = dmsfy
      dmsfzold = dmsfz
    ELSE
      dmsfxold = 1D6
      dmsfyold = 1D6
      dmsfzold = 1D6
    END IF

    dmsfx = 0D0
    dmsfy = 0D0
    dmsfz = 0D0
    dmspox = 0D0
    dmspoy = 0D0
    dmspoz = 0D0

    IF (ALLOCATED(xm)) DEALLOCATE (xm)
    ALLOCATE (xm(1:NE))

    xm(:) = me*ame*1836D0

    CALL leapfr(xe(1:ne), ye(1:ne), ze(1:ne), pxe(1:ne), pye(1:ne), pze(1:ne), dt1/2D0, xm, ne, 2)

    CALL getforces(rho, psidummy, it, 0)

    icoun = 0

    DO i = 1, NE
      IF (imobe(i) == 1) THEN
        icoun = icoun + 1

        fxe(i) = fxe(i) - 1D2*pxe(i)/xm(i)
        fye(i) = fye(i) - 1D2*pye(i)/xm(i)
        fze(i) = fze(i) - 1D2*pze(i)/xm(i)

        dmsfx = dmsfx + fxe(i)**2
        dmsfy = dmsfy + fye(i)**2
        dmsfz = dmsfz + fze(i)**2

        dmspox = (xc(i) - xe(i))**2 + dmspox
        dmspoy = (yc(i) - ye(i))**2 + dmspoy
        dmspoz = (zc(i) - ze(i))**2 + dmspoz

      END IF
    END DO

    xm = 1D0 ! fxe, fye, fze are forces

    CALL leapfr(pxe(1:ne), pye(1:ne), pze(1:ne), fxe(1:ne), fye(1:ne), fze(1:ne), dt1/2D0, xm, ne, 2)

    IF (dmsfxold < dmsfx .OR. dmsfyold < dmsfy .OR. dmsfzold < dmsfz) THEN
      DO i = 1, NE
        pxe(i) = 0D0
        pye(i) = 0D0
        pze(i) = 0D0
      END DO
    END IF

    WRITE (6, '(a,5e18.8,1e22.14)') 'dmsf: ', SQRT(dmsfx)/icoun, &
      SQRT(dmsfy)/icoun, SQRT(dmsfz)/icoun, SQRT(dmspox)/icoun, &
      SQRT(dmspoy)/icoun, SQRT(dmspoz)/icoun

  END DO

  WRITE (6, *) 'Leaving adjustdip'

  RETURN
END SUBROUTINE adjustdip

!--------------------------------------------------------------------------

SUBROUTINE adjustdipz(rho, it)
! computes the force on the GSM-dipoles and adjusts their dipole moments
! IN static iteration (with jdip=1 IN dynamic CASE)
!----------------------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: it

  INTEGER :: iter, icoun, irar
  REAL(DP) :: delmr, delmrold
  REAL(DP) ::rmaxpol, rrrr, posox, posoy, posoz
  COMPLEX(DP) :: psidummy(1)

  REAL(DP), EXTERNAL :: rgetmaxpol, rgetmeanzpol
  WRITE (6, *) 'Entering adjustdipz'

  delmr = 1.0D0
  delmrold = delmr

  DO iter = 1, 200
    IF (SQRT(delmr) < 2D-5) EXIT
    delmr = 0D0
! compute forces of valence clouds
! CALL calcforce_dip(rho)

    WRITE (6, *) 'Entering getforces'
    CALL getforces(rho, psidummy, it, 2)

    icoun = 0

    rmaxpol = rgetmaxpol()
    rrrr = rgetmeanzpol()

    CALL rgetmaxforce(2)

    WRITE (6, '(a,2e17.7)') ' maxpol, meanpol = ', rmaxpol, rrrr
    WRITE (6, '(a,3e17.7)') 'max. fx, fy, fz: ', rvectmp(1), rvectmp(2), rvectmp(3)

    DO irar = 1, NE
      posox = xe(irar)
      posoy = ye(irar)
      posoz = ze(irar)

      IF (imobe(irar) /= 0) THEN
! xe(irar)=xc(irar)+fxe(irar)/cspr !(e2*C_dipmod)
! ye(irar)=yc(irar)+fye(irar)/cspr !idem
        ze(irar) = zc(irar) + fze(irar)/cspr !idem

        icoun = icoun + 1

        delmr = delmr + (posox - xe(irar))**2 + (posoy - ye(irar))**2 &
                + (posoz - ze(irar))**2

      END IF

      IF (ipseudo) THEN
        CALL updatesubgrids
      END IF

    END DO
! delmr=delmr/ne
    delmr = delmr/icoun

    WRITE (6, *) ' delmr = ', SQRT(delmr)

! IF (delmr.gt.delmrold) STOP 'adjustdip not converging'
    delmrold = delmr

  END DO

  IF (ifadiadip == 1) THEN
    DO irar = 1, NE
      pxe(irar) = 0D0
      pye(irar) = 0D0
      pze(irar) = 0D0
    END DO
  END IF

! IF(mod(jter,jbasob).eq.0)THEN
! DO irar=1,nrare
! WRITE(57,'(i4,f7.2,3(a,E15.7))')iter,jter*dt,' ',
! & posrar(irar,1,2)-posrar(irar,1,1),' ',
! & posrar(irar,2,2)-posrar(irar,2,1),' ',
! & posrar(irar,3,2)-posrar(irar,3,1)
! END DO
! END IF

  WRITE (6, *) 'Leaving adjustdipz'

  RETURN
END SUBROUTINE adjustdipz
!------------------------------------------------------------
!

REAL(DP) FUNCTION dgaussdr(r, s)
!------------------------------------------------------------
! Derivative of a Gaussian of width 's' at coordinate 'r'.
  USE params, ONLY: DP
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: s
  REAL(DP)::ra

  ra = r/s
  dgaussdr = -2D0*ra*EXP(-(ra*ra))/s
  RETURN
END FUNCTION dgaussdr

#endif
