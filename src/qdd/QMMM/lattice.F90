
#if(raregas)
!------------------------------------------------------------

SUBROUTINE createouterpot
!------------------------------------------------------------
  USE params
!USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)

  REAL(DP), ALLOCATABLE :: potouter(:)
  ALLOCATE (potouter(kdfull2))

  ind = 0
  DO iz = minz, maxz
    z = (iz - nzsh)*dz
    DO iy = miny, maxy
      y = (iy - nysh)*dy
      DO ix = minx, maxx
        x = (ix - nxsh)*dx

        ind = ind + 1

        potouter(ind) = v_outer(x, y, z)

      END DO
    END DO
  END DO

  DEALLOCATE (potouter)

  RETURN
END SUBROUTINE createouterpot
!------------------------------------------------------------

!~ !------------------------------------------------------------

!~ SUBROUTINE debugmultexp(n,isor)
!~ !------------------------------------------------------------
!~ USE params
!~ !USE kinetic
!~ IMPLICIT REAL(DP) (A-H,O-Z)

!~ ! DIMENSION xatom(ngpar,3)

!~ INTEGER, INTENT(IN) :: n
!~ INTEGER, INTENT(IN OUT) :: isor(ngpar)

!~ REAL(DP) :: xr(3)

!~ DO iz=1,1000

!~ x=5.8D0
!~ y=5.8D0
!~ z=0D0+iz*0.1D0

!~ sumpotex = 0D0
!~ sumpotme = 0D0

!~ iup = 50
!~ idown = 0

!~ DO jx=-iup,iup

!~ DO jy=-iup,iup

!~ IF (ABS(jx) >= idown .OR. ABS(jy) >= idown) THEN

!~ rvx = jx*rlattvec(1)
!~ rvy = jy*rlattvec(2)
!~ rvz = 0

!~ rrm = 0D0
!~ potex = 0D0
!~ potme = 0D0

!~ DO i=1,n
!~ IF (ipsort(i) == 1) THEN
!~ cqq=chgc0
!~ ss=sigmac
!~ ELSE IF (ipsort(i) == 2) THEN
!~ cqq=chge0
!~ ss=sigmav
!~ ELSE IF (ipsort(i) == 3) THEN
!~ cqq=chgk0
!~ ss=sigmak
!~ END IF

!~ rr = (xatom(i,1)-x+rvx)**2+(xatom(i,2)-y+rvy)**2+ &
!~ (xatom(i,3)-z+rvz)**2
!~ rr=MAX(small,SQRT(rr))
!~ rrm = rrm + rr

!~ potex = potex + e2*cqq*v_soft(rr,ss)
!~ sumpotex=sumpotex + e2*cqq*v_soft(rr,ss)
!~ END DO

!~ rrm = rrm / n

!~ sumd = 0D0
!~ sumq = 0D0

!~ jz = 1

!~ xr(1) = x - jx*rlattvec(1)
!~ xr(2) = y - jy*rlattvec(2)
!~ xr(3) = z - (jz-1)*rlattvec(3)

!~ rr = xr(1)*xr(1) + xr(2)*xr(2) + xr(3)*xr(3)
!~ r = SQRT(rr)
!~ rrr = r*rr
!~ rrrrr = rrr*rr

!~ DO ico=1,3
!~ sumd = sumd + xr(i)*celldipole(i)/rrr
!~ END DO

!~ DO ico1=1,3
!~ DO ico2=1,3
!~ sumq = sumq + cellquad(ico1,ico2)*xr(ico1)*xr(ico2) /rrrrr
!~ END DO
!~ END DO

!~ ! WRITE(6,*) 'sum = ',0.5*sumq+sumd

!~ sumq = 0.5D0*sumq

!~ potme = (sumd + sumq)*e2
!~ sumpotme = sumpotme + (sumd + sumq)*e2

!~ ! WRITE(6,'(1f15.3,3e17.7)') rrm , potex, potME,(potME-potEx)/potex

!~ END IF
!~ END DO

!~ END DO
!~ WRITE(6,'(1f10.2,3e17.7)') z,sumpotex,sumpotme, (sumpotme-sumpotex)/sumpotex

!~ WRITE(407,'(1f10.2,3e17.7)') z,sumpotex,sumpotme, &
!~ (sumpotme-sumpotex)/sumpotex

!~ END DO
!~ STOP

!~ RETURN
!~ END SUBROUTINE debugmultexp
!~ !------------------------------------------------------------

!~ !------------------------------------------------------------

!~ REAL(DP) FUNCTION debugfme(n,isor)
!~ !------------------------------------------------------------
!~ USE params
!~ !USE kinetic
!~ IMPLICIT REAL(DP) (A-H,O-Z)

!~ INTEGER, INTENT(IN) :: n
!~ INTEGER, INTENT(IN OUT) :: isor(ngpar)

!~ REAL(DP) :: xr(3)

!~ x=1D0
!~ y=1D0
!~ z=1D0
!~ sumfxex=0D0
!~ sumfyex=0D0
!~ sumfzex=0D0
!~ sumfxme=0D0
!~ sumfyme=0D0
!~ sumfzme=0D0

!~ iup = 20
!~ idown = 10

!~ DO jx=-iup,iup

!~ DO jy=-iup,iup

!~ IF (ABS(jx) >= idown .OR. ABS(jy) >= idown) THEN
!~ rvx = jx*rlattvec(1)
!~ rvy = jy*rlattvec(2)
!~ rvz = 0D0

!~ rrm = 0D0
!~ potex = 0D0
!~ potme = 0D0
!~ fxex=0D0
!~ fyex=0D0
!~ fzex=0D0
!~ fxme=0D0
!~ fyme=0D0
!~ fzme=0D0
!~ rder=1.D-5

!~ DO i=1,n
!~ IF (ipsort(i) == 1) THEN
!~ cqq=chgc0
!~ ss=sigmac
!~ ELSE IF (ipsort(i) == 2) THEN
!~ cqq=chge0
!~ ss=sigmav
!~ ELSE IF (ipsort(i) == 3) THEN
!~ cqq=chgk0
!~ ss=sigmak
!~ END IF

!~ x1=x+rder
!~ y1=y
!~ z1=z
!~ x2=x-rder
!~ y2=y
!~ z2=z

!~ rr1 = (xatom(i,1)-x1+rvx)**2+(xatom(i,2)-y1+rvy)**2+ &
!~ (xatom(i,3)-z1+rvz)**2
!~ rr1=MAX(small,SQRT(rr1))

!~ rr2 = (xatom(i,1)-x2+rvx)**2+(xatom(i,2)-y2+rvy)**2+ &
!~ (xatom(i,3)-z2+rvz)**2
!~ rr2=MAX(small,SQRT(rr2))

!~ radforx = (v_soft(rr1,ss)-v_soft(rr2,ss))/2D0/rder

!~ fxex = fxex - cqq*e2*radforx
!~ sumfxex=sumfxex- cqq*e2*radforx

!~ x1=x
!~ y1=y+rder
!~ z1=z
!~ x2=x
!~ y2=y-rder
!~ z2=z

!~ rr1 = (xatom(i,1)-x1+rvx)**2+(xatom(i,2)-y1+rvy)**2+ &
!~ (xatom(i,3)-z1+rvz)**2
!~ rr1=MAX(small,SQRT(rr1))

!~ rr2 = (xatom(i,1)-x2+rvx)**2+(xatom(i,2)-y2+rvy)**2+ &
!~ (xatom(i,3)-z2+rvz)**2
!~ rr2=MAX(small,SQRT(rr2))

!~ radfory = (v_soft(rr1,ss)-v_soft(rr2,ss))/2D0/rder

!~ fyex = fyex - cqq*e2*radfory
!~ sumfyex=sumfyex- cqq*e2*radfory

!~ x1=x
!~ y1=y
!~ z1=z+rder
!~ x2=x
!~ y2=y
!~ z2=z-rder

!~ rr1 = (xatom(i,1)-x1+rvx)**2+(xatom(i,2)-y1+rvy)**2+ &
!~ (xatom(i,3)-z1+rvz)**2
!~ rr1=MAX(small,SQRT(rr1))

!~ rr2 = (xatom(i,1)-x2+rvx)**2+(xatom(i,2)-y2+rvy)**2+ &
!~ (xatom(i,3)-z2+rvz)**2
!~ rr2=MAX(small,SQRT(rr2))

!~ radforz = (v_soft(rr1,ss)-v_soft(rr2,ss))/2D0/rder

!~ fzex = fzex - cqq*e2*radforz
!~ sumfzex=sumfzex - cqq*e2*radforz

!~ END DO

!~ x1=x+rder
!~ y1=y
!~ z1=z

!~ sumd = 0D0
!~ sumq = 0D0

!~ jz = 1

!~ xr(1) = x1 - jx*rlattvec(1)
!~ xr(2) = y1 - jy*rlattvec(2)
!~ xr(3) = z1 - (jz-1)*rlattvec(3)

!~ rr = xr(1)*xr(1) + xr(2)*xr(2) + xr(3)*xr(3)
!~ r = SQRT(rr)
!~ rrr = r*rr
!~ rrrrr = rrr*rr

!~ DO ico=1,3
!~ sumd = sumd + xr(i)*celldipole(i)/rrr
!~ END DO

!~ DO ico1=1,3
!~ DO ico2=1,3
!~ sumq = sumq + cellquad(ico1,ico2)*xr(ico1)*xr(ico2) /rrrrr
!~ END DO
!~ END DO

!~ ! WRITE(6,*) 'sum = ',0.5*sumq+sumd

!~ sumq = 0.5D0*sumq

!~ potme = (sumd + sumq)*e2

!~ x1=x-rder
!~ y1=y
!~ z1=z

!~ sumd = 0D0
!~ sumq = 0D0

!~ jz = 1

!~ xr(1) = x1 - jx*rlattvec(1)
!~ xr(2) = y1 - jy*rlattvec(2)
!~ xr(3) = z1 - (jz-1)*rlattvec(3)

!~ rr = xr(1)*xr(1) + xr(2)*xr(2) + xr(3)*xr(3)
!~ r = SQRT(rr)
!~ rrr = r*rr
!~ rrrrr = rrr*rr

!~ DO ico=1,3
!~ sumd = sumd + xr(i)*celldipole(i)/rrr
!~ END DO

!~ DO ico1=1,3
!~ DO ico2=1,3
!~ sumq = sumq + cellquad(ico1,ico2)*xr(ico1)*xr(ico2) /rrrrr
!~ END DO
!~ END DO

!~ ! WRITE(6,*) 'sum = ',0.5*sumq+sumd

!~ sumq = 0.5D0*sumq

!~ potme = potme - (sumd + sumq)*e2

!~ fxme = fxme - potme/2D0/rder
!~ sumfxme=sumfxme - potme/2D0/rder

!~ potme = 0D0

!~ x1=x
!~ y1=y+rder
!~ z1=z

!~ sumd = 0D0
!~ sumq = 0D0

!~ jz = 1

!~ xr(1) = x1 - jx*rlattvec(1)
!~ xr(2) = y1 - jy*rlattvec(2)
!~ xr(3) = z1 - (jz-1)*rlattvec(3)

!~ rr = xr(1)*xr(1) + xr(2)*xr(2) + xr(3)*xr(3)
!~ r = SQRT(rr)
!~ rrr = r*rr
!~ rrrrr = rrr*rr

!~ DO ico=1,3
!~ sumd = sumd + xr(i)*celldipole(i)/rrr
!~ END DO

!~ DO ico1=1,3
!~ DO ico2=1,3
!~ sumq = sumq + cellquad(ico1,ico2)*xr(ico1)*xr(ico2) /rrrrr
!~ END DO
!~ END DO

!~ ! WRITE(6,*) 'sum = ',0.5*sumq+sumd

!~ sumq = 0.5D0*sumq

!~ potme = (sumd + sumq)*e2

!~ x1=x
!~ y1=y-rder
!~ z1=z

!~ sumd = 0D0
!~ sumq = 0D0

!~ jz = 1

!~ xr(1) = x1 - jx*rlattvec(1)
!~ xr(2) = y1 - jy*rlattvec(2)
!~ xr(3) = z1 - (jz-1)*rlattvec(3)

!~ rr = xr(1)*xr(1) + xr(2)*xr(2) + xr(3)*xr(3)
!~ r = SQRT(rr)
!~ rrr = r*rr
!~ rrrrr = rrr*rr

!~ DO ico=1,3
!~ sumd = sumd + xr(i)*celldipole(i)/rrr
!~ END DO

!~ DO ico1=1,3
!~ DO ico2=1,3
!~ sumq = sumq + cellquad(ico1,ico2)*xr(ico1)*xr(ico2) /rrrrr
!~ END DO
!~ END DO

!~ ! WRITE(6,*) 'sum = ',0.5*sumq+sumd

!~ sumq = 0.5D0*sumq

!~ potme = potme - (sumd + sumq)*e2

!~ fyme = fyme - potme/2D0/rder
!~ sumfyme=sumfyme - potme/2D0/rder

!~ potme = 0D0

!~ x1=x
!~ y1=y
!~ z1=z+rder

!~ sumd = 0D0
!~ sumq = 0D0

!~ jz = 1

!~ xr(1) = x1 - jx*rlattvec(1)
!~ xr(2) = y1 - jy*rlattvec(2)
!~ xr(3) = z1 - (jz-1)*rlattvec(3)

!~ rr = xr(1)*xr(1) + xr(2)*xr(2) + xr(3)*xr(3)
!~ r = SQRT(rr)
!~ rrr = r*rr
!~ rrrrr = rrr*rr

!~ DO ico=1,3
!~ sumd = sumd + xr(i)*celldipole(i)/rrr
!~ END DO

!~ DO ico1=1,3
!~ DO ico2=1,3
!~ sumq = sumq + cellquad(ico1,ico2)*xr(ico1)*xr(ico2) /rrrrr
!~ END DO
!~ END DO

!~ ! WRITE(6,*) 'sum = ',0.5*sumq+sumd

!~ sumq = 0.5D0*sumq

!~ potme = (sumd + sumq)*e2

!~ x1=x
!~ y1=y
!~ z1=z-rder

!~ sumd = 0D0
!~ sumq = 0D0

!~ jz = 1

!~ xr(1) = x1 - jx*rlattvec(1)
!~ xr(2) = y1 - jy*rlattvec(2)
!~ xr(3) = z1 - (jz-1)*rlattvec(3)

!~ rr = xr(1)*xr(1) + xr(2)*xr(2) + xr(3)*xr(3)
!~ r = SQRT(rr)
!~ rrr = r*rr
!~ rrrrr = rrr*rr

!~ DO ico=1,3
!~ sumd = sumd + xr(i)*celldipole(i)/rrr
!~ END DO

!~ DO ico1=1,3
!~ DO ico2=1,3
!~ sumq = sumq + cellquad(ico1,ico2)*xr(ico1)*xr(ico2) /rrrrr
!~ END DO
!~ END DO

!~ ! WRITE(6,*) 'sum = ',0.5*sumq+sumd

!~ sumq = 0.5D0*sumq

!~ potme = potme - (sumd + sumq)*e2

!~ fzme = fzme - potme/2D0/rder
!~ sumfzme=sumfzme - potme/2D0/rder

!~ WRITE(6,'(1f12.3,6e16.5)') rrm, fxex,fxme,fyex,fyme,fzex,fzme

!~ ! WRITE(6,'(1f15.3,3e17.7)') rrm , potex, potME,(potME-potEx)/potME

!~ END IF
!~ END DO

!~ END DO

!~ WRITE(6,'(3e17.7)') sumfxex,sumfxme,(sumfxme-sumfxex)/sumfxex
!~ WRITE(6,'(3e17.7)') sumfyex,sumfyme,(sumfyme-sumfyex)/sumfyex
!~ WRITE(6,'(3e17.7)') sumfyex,sumfyme,(sumfzme-sumfzex)/sumfzex
!~ STOP

!~ RETURN
!~ END FUNCTION debugfme
!~ !------------------------------------------------------------

!------------------------------------------------------------

FUNCTION v_outer(x, y, z)
!------------------------------------------------------------
  USE params
!USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)

  REAL, INTENT(IN) :: x
  REAL, INTENT(IN) :: y
  REAL, INTENT(IN) :: z
  REAL(DP) :: xr(3)

  COMPLEX(DP) :: sumo

  sumo = CMPLX(0D0, 0D0, DP)

  sumd = 0D0
  sumq = 0D0

! WRITE(6,*) 'HERE', ncells1x,ncells2x

  DO jx = -(ncells1x + ncells2x), (ncells1x + ncells2x)

    DO jy = -(ncells1y + ncells2y), (ncells1y + ncells2y)

      DO jz = ncells1z + 1, ncells2z + ncells1z + 1

        IF (ABS(jx) > ncells1y .OR. ABS(jy) > ncells1z) THEN

          xr(1) = x - jx*rlattvec(1)
          xr(2) = y - jy*rlattvec(2)
          xr(3) = z - (jz - 1)*rlattvec(3)

          rr = xr(1)*xr(1) + xr(2)*xr(2) + xr(3)*xr(3)
          r = SQRT(rr)
          rrr = r*rr
          rrrrr = rrr*rr

          DO ico = 1, 3
            sumd = sumd + xr(i)*celldipole(i)/rrr
          END DO

          DO ico1 = 1, 3
            DO ico2 = 1, 3
              sumq = sumq + cellquad(ico1, ico2)*xr(ico1)*xr(ico2)/rrrrr
            END DO
          END DO

!$$$ th = acos(xr(3)/r)
!$$$
!$$$ IF (abs(th).gt.1e-15) THEN
!$$$ ph = acos(xr(1)/r/sin(th))
!$$$ ELSE
!$$$ ph = 0.
!$$$ END IF
!$$$
!$$$
!$$$ DO ilm=-2,2
!$$$
!$$$ sumo = sumo + 0.8*PI*cellMult(2,ilm)*spherHarm(2,ilm,
!$$$ & th,ph)/rrr
!$$$
!$$$ END DO
!$$$
!$$$ WRITE(6,*) 'sum = ',0.5*sumq+sumd
!$$$ WRITE(6,*) 'sumo = ',sumo

        END IF
      END DO
    END DO
  END DO

  sumq = 0.5D0*sumq

  v_outer = sumd + sumq

! WRITE(6,*) 'STOP PROGRAM'
! STOP

  RETURN
END FUNCTION v_outer
!------------------------------------------------------------

#if(raregas)
!------------------------------------------------------------

SUBROUTINE buildlatt(iside, ilayers)
!------------------------------------------------------------
! constructs MgO lattice according TO parameters IN FILE box.par
  USE params
!USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)

  WRITE (6, *) 'iside,ilayers', iside, ilayers
  ischicht = 0
  icntc = 0
  icntk = 0

  IF (MOD(iside, 2) == 1) THEN
    STOP 'nside must be even for neutral box!'
  END IF

  m = iside

  DO k = 1, ilayers
    ischicht = ischicht + 1

    IF (MOD(ischicht, 2) == 0) THEN
      itypion = 0
    ELSE
      itypion = 1
    END IF

    z = zsurf - (k - 1)*a

    DO i = 1, m

      itypion = itypion + 1

      IF (MOD(itypion, 2) == 0) THEN
        itypion2 = 0
      ELSE
        itypion2 = 1
      END IF

      y = 0.5D0*a*(m - 1) - (i - 1)*a

      DO j = 1, m
        itypion2 = itypion2 + 1

        x = 0.5D0*a*(m - 1) - (j - 1)*a

        IF (MOD(itypion2, 2) == 0) THEN ! anion

          icntc = icntc + 1
          xc(icntc) = x
          yc(icntc) = y
          zc(icntc) = z
          chgc(icntc) = chgcore

          xe(icntc) = x + delx
          ye(icntc) = y + dely
          ze(icntc) = z + delz
          chge(icntc) = chgval

        ELSE ! cation

          icntk = icntk + 1

          chgk(icntk) = chgkat
          xk(icntk) = x
          yk(icntk) = y
          zk(icntk) = z

        END IF

      END DO

    END DO

  END DO

  nc = icntc
  NE = nc
  nk = icntk

  dlx = iside*a
  dly = dlx
  dlz = ilayers*a

! WRITE (6,*) 'Grid construction done.... writing TO FILE...'

! OPEN (43,STATUS='unknown',FILE='gitter.res')

! itmp = icntc+icntk

! WRITE (43,'(1i10.1,4f10.1)') itmp,dLx,dLy,dLz,zsurf

! chganion = chgcore + chgval

! DO i=1,icntc
! WRITE (43,'(3f13.3,1f10.1)') xc(i),yc(i),zc(i),chganion
!c WRITE (6,'(3f13.3,1f10.1)') x(i),y(i),z(i),chganion
! END DO

! DO i=1,icntk
! WRITE (43,'(3f13.3,1f10.1)') xk(i),yk(i),zk(i),chgkat
!c WRITE (6,'(3f13.3,1f10.1)') xk(i),yk(i),zk(i),chgkat
! END DO

! CLOSE(43)

! WRITE(6,*) ianions,ikations,chgtot

END SUBROUTINE buildlatt
!------------------------------------------------------------
#endif

#endif
