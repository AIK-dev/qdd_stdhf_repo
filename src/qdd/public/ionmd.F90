
!------------------------------------------------------------

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


  enii = sumion


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



! propagation of positions


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

! compute forces on molecule ions with new positions


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

  IF (ityp == 4) THEN
    DO i = 1, n
      x(i) = x(i) + xprop(i)*(ddt/xm(i))
      y(i) = y(i) + yprop(i)*(ddt/xm(i))
      z(i) = z(i) + zprop(i)*(ddt/xm(i))
    END DO
  END IF

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

  IF (ttestprint) WRITE (*, *) 'enter ITSTEPV: time[fs],cpx,cpy,cpz=', &
    it*dt1*0.0484D0, cpx(nion), cpy(nion), cpz(nion)


! PROPAGATION OF POSITIONS


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


! second half of propagation of momenta


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

  IF (ityp == 4) THEN
    DO i = 1, n
      px(i) = px(i) + ddth*fox(i)
      py(i) = py(i) + ddth*foy(i)
      pz(i) = pz(i) + ddth*foz(i)
    END DO
  END IF

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


  CALL getforces(rho, psi, -1, 0)


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

