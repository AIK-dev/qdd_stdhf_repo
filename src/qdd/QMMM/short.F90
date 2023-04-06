
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! The Erlangen - Toulouse

! SHORT RANGE INTERACTION PACKAGE

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!------------------------------------------------------------

SUBROUTINE getshortforce(itypi, itypj, ii, jj, rho, iflag2)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! distributing short force routine
! ind1 and ind2 are of the short TYPE!

  INTEGER, INTENT(IN) :: itypi
  INTEGER, INTENT(IN) :: itypj
  INTEGER, INTENT(IN) :: ii
  INTEGER, INTENT(IN) :: jj
  REAL(DP), INTENT(IN) :: rho(kdfull2*2)
  INTEGER, INTENT(IN) :: iflag2

  INTEGER :: ityp1, ityp2
  INTEGER :: ind1, ind2
  REAL(DP) :: rder
#if(raregas)
  INTEGER :: ind
  INTEGER :: itmp, itmp1, itmp2
  INTEGER :: ix, iy, iz
  INTEGER :: nsgsize
  REAL(DP) :: chpddr, dist, dist2
  REAL(DP) :: radfor, signum, sumx, sumy, sumz
  REAL(DP) :: r2, rr, xr, yr, zr, x1, y1, z1
#endif

  INTEGER, EXTERNAL :: conv3to1, getnearestgridpoint, iconvshorttolong
  REAL(DP), EXTERNAL :: fbornmayer, fbornmayermod, getdistance2, vararlj, &
                        varnashort, varelcore, v_ar_ar

  ityp1 = itypi !Because the RESULT of a FUNCTION or a scalar expression cannot be INTENT(OUT)/INTENT(IN OUT) arguments.
  ityp2 = itypj
  ind1 = ii ! Because ii and jj can be DO-variables, that should not be modified inside a DO-loop.
  ind2 = jj

  rder = 1.0D-5

#if(raregas)

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IF (isrtyp(ityp1, ityp2) == 1) THEN ! MgO-CASE

    signum = 1D0

! determine parameters

    IF (ityp1 <= 3 .AND. ityp2 <= 3) THEN ! GSM-GSM is of Born-Mayer-TYPE

      itmp1 = iconvshorttolong(ityp1, ind1)
      itmp2 = iconvshorttolong(ityp2, ind2)

      rr = getdistance2(itmp1, itmp2)
      rr = SQRT(MAX(rr, 1D-14))

      CALL setbornparas(ityp1, ityp2) ! paras are stored IN rVecTmp

      radfor = fbornmayer(rr, rvectmp(1), rvectmp(2), rvectmp(3))

      CALL getrelvec(itmp1, itmp2, rr) ! direction stored IN rVecTmp

      forx = radfor*rvectmp(1)*signum
      fory = radfor*rvectmp(2)*signum
      forz = radfor*rvectmp(3)*signum

      IF (iflag2 == 0) THEN ! as usual
        CALL addforce(itmp1, -forx, -fory, -forz)
        CALL addforce(itmp2, forx, fory, forz)
      ELSE IF (iflag2 == 1) THEN ! ibh second part
        CALL addforce(itmp1, -forx, -fory, -forz)
      ELSE
        CALL addforce(itmp1, forx, fory, forz)
      END IF

    ELSE IF ((ityp1 <= 3 .AND. ityp2 == 4) .OR. &
             (ityp1 == 4 .AND. ityp2 <= 3)) THEN !GSM-Na or Na-Na force

      IF (ityp1 == 1) THEN !

        xr = xc(ind1) - cx(ind2)
        yr = yc(ind1) - cy(ind2)
        zr = zc(ind1) - cz(ind2)

        rr = MAX(SQRT(xr*xr + yr*yr + zr*zr), small)

        radfor = fbornmayermod(rr, bcn, ecn, fcn, sigcn, ccncud, ccnd, dcn)

        forx = radfor*xr/rr*signum
        fory = radfor*yr/rr*signum
        forz = radfor*zr/rr*signum

        fxc(ind1) = fxc(ind1) - forx
        fyc(ind1) = fyc(ind1) - fory
        fzc(ind1) = fzc(ind1) - forz

        fx(ind2) = fx(ind2) + forx
        fy(ind2) = fy(ind2) + fory
        fz(ind2) = fz(ind2) + forz
      ELSE IF (ityp2 == 1) THEN
        xr = cx(ind1) - xc(ind2)
        yr = cy(ind1) - yc(ind2)
        zr = cz(ind1) - zc(ind2)

        rr = MAX(SQRT(xr*xr + yr*yr + zr*zr), small)

        radfor = fbornmayermod(rr, bcn, ecn, fcn, sigcn, ccncud, ccnd, dcn)

        forx = radfor*xr/rr*signum
        fory = radfor*yr/rr*signum
        forz = radfor*zr/rr*signum

        fx(ind1) = fx(ind1) - forx
        fy(ind1) = fy(ind1) - fory
        fz(ind1) = fz(ind1) - forz

        fxc(ind2) = fxc(ind2) + forx
        fyc(ind2) = fyc(ind2) + fory
        fzc(ind2) = fzc(ind2) + forz

      ELSE IF (ityp1 == 3) THEN

        xr = xk(ind1) - cx(ind2)
        yr = yk(ind1) - cy(ind2)
        zr = zk(ind1) - cz(ind2)

        rr = MAX(SQRT(xr*xr + yr*yr + zr*zr), small)

        radfor = fbornmayermod(rr, bkn, ekn, fkn, sigkn, ckncud, cknd, dkn)

        forx = radfor*xr/rr*signum
        fory = radfor*yr/rr*signum
        forz = radfor*zr/rr*signum

        fxk(ind1) = fxk(ind1) - forx
        fyk(ind1) = fyk(ind1) - fory
        fzk(ind1) = fzk(ind1) - forz

        fx(ind2) = fx(ind2) + forx
        fy(ind2) = fy(ind2) + fory
        fz(ind2) = fz(ind2) + forz

      ELSE IF (ityp2 == 3) THEN

        xr = cx(ind1) - xk(ind2)
        yr = cy(ind1) - yk(ind2)
        zr = cz(ind1) - zk(ind2)

        rr = MAX(SQRT(xr*xr + yr*yr + zr*zr), small)

        radfor = fbornmayermod(rr, bkn, ekn, fkn, sigkn, ckncud, cknd, dkn)

        forx = radfor*xr/rr*signum
        fory = radfor*yr/rr*signum
        forz = radfor*zr/rr*signum

        fx(ind1) = fx(ind1) - forx
        fy(ind1) = fy(ind1) - fory
        fz(ind1) = fz(ind1) - forz

        fxk(ind2) = fxk(ind2) + forx
        fyk(ind2) = fyk(ind2) + fory
        fzk(ind2) = fzk(ind2) + forz

      ELSE ! Na-Na
! no short range interaction as IN MgO CASE
      END IF

    ELSE IF (ityp2 == 5) THEN ! GSM-DFT force

      IF (ityp1 <= 3) THEN

! TO be implemented by fitting the Muenchen DATA

      ELSE
! nothing TO DO
      END IF

    ELSE
      STOP 'Error IN MgO CASE of getShortForce'
    END IF

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(ityp1, ityp2) == 2) THEN ! Ar-CASE

    signum = 1D0

! sort by particle TYPE
    IF (ityp1 > ityp2) THEN
      itmp = ityp1
      ityp1 = ityp2
      ityp2 = itmp
      itmp = ind1
      ind1 = ind2
      ind2 = itmp
      signum = -1D0
      STOP 'warning: signum is -1'
    END IF

    IF (ityp1 == 1) THEN ! 1=c
      IF (ityp2 == 1) THEN ! 2=c
! Ar CASE: Core-Core Interaction

        xr = xc(ind1) - xc(ind2)
        yr = yc(ind1) - yc(ind2)
        zr = zc(ind1) - zc(ind2)

        dist2 = xr*xr + yr*yr + zr*zr
        dist = SQRT(dist2)

! radfor = VArArShort(dist+rder)-VArArShort(dist-rder)
        IF (iararlj == 0) THEN
          radfor = v_ar_ar(dist + rder) - v_ar_ar(dist - rder)
        ELSE
          radfor = vararlj(dist + rder) - vararlj(dist - rder)
        END IF

        radfor = radfor/2D0/rder*signum

        forx = radfor*xr/dist
        fory = radfor*yr/dist
        forz = radfor*zr/dist

        fxc(ind1) = fxc(ind1) - forx
        fyc(ind1) = fyc(ind1) - fory
        fzc(ind1) = fzc(ind1) - forz

        fxc(ind2) = fxc(ind2) + forx
        fyc(ind2) = fyc(ind2) + fory
        fzc(ind2) = fzc(ind2) + forz

      ELSE IF (ityp2 == 4) THEN ! 2=Na
! Ar CASE: Core-Na interaction

        xr = xc(ind1) - cx(ind2)
        yr = yc(ind1) - cy(ind2)
        zr = zc(ind1) - cz(ind2)

! WRITE(6,*) 'xryrzr:',xr,yr,zr

        dist2 = xr*xr + yr*yr + zr*zr
        dist = SQRT(dist2)

        radfor = varnashort(dist + rder) - varnashort(dist - rder)

        radfor = radfor/(rder + rder)*signum

        forx = radfor*xr/dist
        fory = radfor*yr/dist
        forz = radfor*zr/dist

        fxc(ind1) = fxc(ind1) - forx
        fyc(ind1) = fyc(ind1) - fory
        fzc(ind1) = fzc(ind1) - forz

        fx(ind2) = fx(ind2) + forx
        fy(ind2) = fy(ind2) + fory
        fz(ind2) = fz(ind2) + forz

      ELSE IF (ityp2 == 5) THEN !2=DFT
! Ar CASE: Core-electron interaction

        ind = 0
        sumx = 0D0
        sumy = 0D0
        sumz = 0D0

        IF (.NOT. ipseudo) THEN

          DO iz = minz, maxz
            z1 = (iz - nzsh)*dz
            zr = z1 - zc(ind1)

            DO iy = miny, maxy
              y1 = (iy - nysh)*dy
              yr = y1 - yc(ind1)

              DO ix = minx, maxx
                x1 = (ix - nxsh)*dx
                xr = x1 - xc(ind1)

                ind = ind + 1

                r2 = xr*xr + yr*yr + zr*zr
                rr = MAX(SQRT(r2), small)

                chpddr = rho(ind)*(varelcore(rr + rder) &
                                   - varelcore(rr - rder))/(rr*2.*rder)

! attention again: rho is positive!!!

                sumx = sumx - chpddr*xr
                sumy = sumy - chpddr*yr
                sumz = sumz - chpddr*zr

              END DO
            END DO
          END DO

        ELSE ! ipseudo = .TRUE.

          ind = getnearestgridpoint(xc(ind1), yc(ind1), zc(ind1))

          CALL conv1to3(ind)

          nsgsize = 10 ! just for the argon CASE!

          DO iz = iindtmp(3) - nsgsize, iindtmp(3) + nsgsize
            IF (iz >= minz .AND. iz <= maxz) THEN
              z1 = (iz - nzsh)*dz
              zr = zc(ind1) - z1

              DO iy = iindtmp(2) - nsgsize, iindtmp(2) + nsgsize
                IF (iy >= miny .AND. iy <= maxy) THEN

                  y1 = (iy - nysh)*dy
                  yr = yc(ind1) - y1

                  DO ix = iindtmp(1) - nsgsize, iindtmp(1) + nsgsize
                    IF (ix >= minx .AND. ix <= maxx) THEN

                      x1 = (ix - nxsh)*dx
                      xr = xc(ind1) - x1

                      ind = conv3to1(ix, iy, iz)

                      r2 = xr*xr + yr*yr + zr*zr
                      rr = MAX(SQRT(r2), small)

                      chpddr = rho(ind)*(varelcore(rr + rder) &
                                         - varelcore(rr - rder))/(rr*2D0*rder)

! chpddr = rho(ind)*(V_Ar_el_core(rr+rder)
! & - V_Ar_el_core(rr-rder))/(rr*2.*rder)
! attention again: rho is positive!!!

                      sumx = sumx + chpddr*xr
                      sumy = sumy + chpddr*yr
                      sumz = sumz + chpddr*zr

                    END IF
                  END DO
                END IF
              END DO
            END IF
          END DO

        END IF

        fxc(ind1) = fxc(ind1) + sumx*dvol
        fyc(ind1) = fyc(ind1) + sumy*dvol
        fzc(ind1) = fzc(ind1) + sumz*dvol

      END IF
    END IF

  ELSE !other cases

! OTHER CASES TO BE IMPLEMENTED HERE

  END IF
#endif

  RETURN
END SUBROUTINE getshortforce
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE check_isrtyp
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! Checks consistency of the short range interaction TYPE matrix
! isrtyp(i,j)

! i,j each range from 1..5, meaning:

! i ---> 1 ---> GSM core
! 2 ---> GSM valence shell
! 3 ---> GSM cation
! 4 ---> Na core
! 5 ---> DFT electron

! The value of the matrix elements can be:

! 0 ---> no short range interaction
! 1 ---> Born-Mayer TYPE
! 2 ---> Argon CASE

! The matrix must be symmetric, of course.

#if(raregas)
  INTEGER :: i, j

  IF (isrtypall /= 0) THEN
    WRITE (6, *) 'Setting short range parameters...'
    DO i = 1, 5
      DO j = 1, 5
        isrtyp(i, j) = isrtypall
      END DO
    END DO

  END IF

  DO i = 1, 5
    DO j = 1, 5
      IF (isrtyp(i, j) /= isrtyp(j, i)) THEN
        STOP 'Short range interaction TYPE matrix not symmetric'
      END IF
    END DO
  END DO

! check IF it is necessary TO handle FUNCTION parameters by input FILE

  DO i = 1, 5
    DO j = i, 5
      IF (isrtyp(i, j) == 3) THEN
        IF (fermia == 0D0) STOP 'no parameters for fermi FUNCTION given'
      END IF
    END DO
  END DO
#endif

  RETURN
END SUBROUTINE check_isrtyp
!------------------------------------------------------------

#if(raregas)
!------------------------------------------------------------

SUBROUTINE getrelvec(ind1, ind2, r)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE
! Calculates the direction vector between part ind1 and ind2, divided by length r.
! output direction stored IN rVecTmp
  INTEGER, INTENT(IN) :: ind1
  INTEGER, INTENT(IN) :: ind2
  REAL(DP), INTENT(IN) :: r

  REAL(DP) :: x1, y1, z1

  CALL getparas(ind1)
  x1 = rvectmp(1)
  y1 = rvectmp(2)
  z1 = rvectmp(3)

  CALL getparas(ind2)

  rvectmp(1) = (x1 - rvectmp(1))/r
  rvectmp(2) = (y1 - rvectmp(2))/r
  rvectmp(3) = (z1 - rvectmp(3))/r

  RETURN
END SUBROUTINE getrelvec
!------------------------------------------------------------

!-----VArArShort------------------------------------------------
! ONLY short range part of Ar-Ar interaction

REAL(DP) FUNCTION vararshort(r)

  USE params
  IMPLICIT NONE

! Ar-Ar potential

  REAL(DP), INTENT(IN) :: r

  REAL(DP) :: core, rabs
  REAL(DP) :: alpha_ar, beta_ar, c6_ar, c8_ar, a_ar
  DATA alpha_ar, beta_ar, c6_ar, c8_ar, a_ar &
    /1.7301D0, 1.7966D0, 55.465D0, 3672.9D0, 794.21D0/

  rabs = MAX(ABS(r), small)
  core = a_ar*EXP(-alpha_ar*rabs)/rabs - 2D0/(1D0 + EXP((beta_ar/rabs)**2)) &
         *(c6_ar/rabs**6 + c8_ar/rabs**8)

  vararshort = e2*core

  RETURN
END FUNCTION vararshort
!------------------------------------------------------------

!-----VArArLJ------------------------------------------------
! ONLY short range part of Ar-Ar interaction

REAL(DP) FUNCTION vararlj(r)

  USE params
  IMPLICIT NONE

! Ar-Ar potential

  REAL(DP), INTENT(IN) :: r
  REAL(DP) :: rabs
  REAL(DP) :: epslj, sigmalj
  DATA epslj, sigmalj/0.0007647D0, 6.42503D0/

  rabs = MAX(ABS(r), small)

  vararlj = 4.0D0*epslj*((sigmalj/rabs)**12 - (sigmalj/rabs)**6)

  RETURN
END FUNCTION vararlj
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION vbornmayer(r, a, s, c)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: a
  REAL(DP), INTENT(IN) :: s
  REAL(DP), INTENT(IN) :: c

  vbornmayer = a*EXP(-r/s) - c/r**6

  RETURN
END FUNCTION vbornmayer
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION fbornmayer(r, a, s, c)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: a
  REAL(DP), INTENT(IN) :: s
  REAL(DP), INTENT(IN) :: c

  REAL(DP):: rder

  rder = 1D-5

! ftmp = VBornMayer(r+rder,A,s,C) - VBornMayer(r-rder,A,s,C)
! ftmp = ftmp/(2.0*rder)

  fbornmayer = -a/s*EXP(-r/s) + 6*c*(erf(r)/r)**7

! FBornMayer = ftmp

  RETURN
END FUNCTION fbornmayer
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION fbornmayer46(r, a, s, c4, c6)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: a
  REAL(DP), INTENT(IN) :: s
  REAL(DP), INTENT(IN) :: c4
  REAL(DP), INTENT(IN) :: c6

  fbornmayer46 = -a/s*EXP(-r/s) + 4*c4/r**5 + 6*c6/r**7

! FBornMayer = ftmp

  RETURN
END FUNCTION fbornmayer46
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION fbornmayer2345678d(r, a, e, f, s, s2, c2, c3, c4, c5, c6, c7, c8, c9, c10, cd, d)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: a
  REAL(DP), INTENT(IN) :: e
  REAL(DP), INTENT(IN) :: f
  REAL(DP), INTENT(IN) :: s
  REAL(DP), INTENT(IN) :: s2
  REAL(DP), INTENT(IN) :: c2
  REAL(DP), INTENT(IN) :: c3
  REAL(DP), INTENT(IN) :: c4
  REAL(DP), INTENT(IN) :: c5
  REAL(DP), INTENT(IN) :: c6
  REAL(DP), INTENT(IN) :: c7
  REAL(DP), INTENT(IN) :: c8
  REAL(DP), INTENT(IN) :: c9
  REAL(DP), INTENT(IN) :: c10
  REAL(DP), INTENT(IN) :: cd
  REAL(DP), INTENT(IN) :: d

  REAL(DP) :: acc, rder, rrp, rrm, sum2
  rder = 1D-5

  acc = 0D0

  sum2 = -a/s*EXP((r - f)/s)/(e + EXP((r - f)/s))**2

  rrp = r + rder
  rrm = r - rder

  IF (c2 /= 0D0) acc = acc + c2*((erf(rrp/s2)/rrp)**2 - (erf(rrm/s2)/rrm)**2)
  IF (c3 /= 0D0) acc = acc + c3*((erf(rrp/s2)/rrp)**3 - (erf(rrm/s2)/rrm)**3)
  IF (c4 /= 0D0) acc = acc + c4*((erf(rrp/s2)/rrp)**4 - (erf(rrm/s2)/rrm)**4)
  IF (c5 /= 0D0) acc = acc + c5*((erf(rrp/s2)/rrp)**5 - (erf(rrm/s2)/rrm)**5)
  IF (c6 /= 0D0) acc = acc + c6*((erf(rrp/s2)/rrp)**6 - (erf(rrm/s2)/rrm)**6)
  IF (c7 /= 0D0) acc = acc + c7*((erf(rrp/s2)/rrp)**7 - (erf(rrm/s2)/rrm)**7)
  IF (c8 /= 0D0) acc = acc + c8*((erf(rrp/s2)/rrp)**8 - (erf(rrm/s2)/rrm)**8)
  IF (c9 /= 0D0) acc = acc + c9*((erf(rrp/s2)/rrp)**9 - (erf(rrm/s2)/rrm)**9)
  IF (c10 /= 0D0) acc = acc + c10*((erf(rrp/s2)/rrp)**10 - (erf(rrm/s2)/rrm)**10)
  IF (cd /= 0D0) acc = acc + cd*((erf(rrp/s2)/rrp)**d - (erf(rrm/s2)/rrm)**d)

  acc = acc/2D0/rder

!$$$ IF (C2.ne.0.0) acc=acc+2*C2/r**3
!$$$ IF (C3.ne.0.0) acc=acc+3*C3/r**4
!$$$ IF (C4.ne.0.0) acc=acc+4*C4/r**5
!$$$ IF (C5.ne.0.0) acc=acc+5*C5/r**6
!$$$ IF (C6.ne.0.0) acc=acc+6*C6/r**7
!$$$ IF (C7.ne.0.0) acc=acc+7*C7/r**8
!$$$ IF (C8.ne.0.0) acc=acc+8*C8/r**9
!$$$ IF (C9.ne.0.0) acc=acc+9*C9/r**10
!$$$ IF (C10.ne.0.0) acc=acc+10*C10/r**11
!$$$ IF (CD.ne.0.0) acc=acc+D*CD/r**(D+1)

  fbornmayer2345678d = acc + sum2

  RETURN
END FUNCTION fbornmayer2345678d
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION fbornmayermod(r, a, e, f, s, s2, cd, d)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: a
  REAL(DP), INTENT(IN) :: e
  REAL(DP), INTENT(IN) :: f
  REAL(DP), INTENT(IN) :: s
  REAL(DP), INTENT(IN) :: s2
  REAL(DP), INTENT(IN) :: cd
  REAL(DP), INTENT(IN) :: d

  REAL(DP) :: acc, rder, rrm, rrp, sum2
  rder = 1D-5

  acc = 0D0

  sum2 = -a/s*EXP((r - f)/s)/(e + EXP((r - f)/s))**2

  rrp = r + rder
  rrm = r - rder

  IF (cd /= 0D0) THEN
    acc = acc - cd*((erf(rrp/s2)/rrp)**d - (erf(rrm/s2)/rrm)**d)
  END IF

  acc = acc/2D0/rder

  fbornmayermod = acc + sum2

  RETURN
END FUNCTION fbornmayermod
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION varnashort(r)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

! Ar-Na potential

! parameters for repulsive core (energies IN Ha, rescale TO Ry by e2 !)

! DATA vrep/214.14/ ! strength of repulsive core
! DATA rcor/1.6906/ ! inverse radius of repulsive core
! DATA vdip2/308.98/ ! 2*strength of higher dipole part
! DATA rdip2/2.20/ ! radius of cutoff IN higher dipole part

!----------------------------------------------------------------------

! rabs = abs(r)

! V_Ar_Na =
! & e2*(vrep*exp(-rcor*rabs)/rabs
! & - vdip2/((1.0d0+exp(rdip2/rabs)**2)*rabs**6)
! & + ch(18)*ch(11)*V_soft(r,sgm1(18))
! & )

  REAL(DP), INTENT(IN) :: r

  REAL(DP) :: a1, a2, a3, b1, b2

  DATA a1/334.85D0/
  DATA a2/52.5D0/
  DATA a3/1383.0D0/
  DATA b1/-1.7624D0/
  DATA b2/1.815D0/

  varnashort = a1*EXP(b1*r)/r - 2.0D0/(1.d0 + EXP(b2/r))*(a2/r**6 + a3/r**8)
  varnashort = e2*varnashort

  RETURN
END FUNCTION varnashort

!------------------------------------------------------------

SUBROUTINE getshortenergy(ityp1, ityp2, ind1, ind2)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: ityp1
  INTEGER, INTENT(INOUT) :: ityp2
  INTEGER, INTENT(INOUT) :: ind1
  INTEGER, INTENT(INOUT) :: ind2

  INTEGER :: itmp, itmp1, itmp2
  REAL(DP) :: dist, dist2, rr, signum, sumion

  INTEGER, EXTERNAL :: iconvshorttolong
  REAL(DP), EXTERNAL :: getdistance2, vbornmayer, vararshort, varnashort, vararlj
! sort by particle TYPE
  IF (ityp1 > ityp2) THEN
    itmp = ityp1
    ityp1 = ityp2
    ityp2 = itmp
    itmp = ind1
    ind1 = ind2
    ind2 = itmp
    signum = -1D0
! STOP 'warning: signum is -1'
  END IF

  sumion = 0D0

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IF (isrtyp(ityp1, ityp2) == 1) THEN ! MgO-CASE

    IF (ityp2 <= 3) THEN ! GSM-GSM CASE
      itmp1 = iconvshorttolong(ityp1, ind1)
      itmp2 = iconvshorttolong(ityp2, ind2)
      rr = getdistance2(itmp1, itmp2)
      rr = SQRT(MAX(rr, 1D-14))

      CALL setbornparas(ityp1, ityp2) ! paras are stored IN rVecTmp

      sumion = sumion + vbornmayer(rr, rvectmp(1), rvectmp(2), rvectmp(3))

    ELSE IF (ityp2 == 4) THEN ! GSM-Na or Na-Na

      IF (ityp1 <= 3) THEN

        IF (ityp1 == 1) THEN
          rr = (xc(ind1) - cx(ind2))**2 + (yc(ind1) - cy(ind2))**2 + &
               (zc(ind1) - cz(ind2))**2
          rr = MAX(small, SQRT(rr))
          sumion = sumion + bcn/(ecn + EXP((rr - fcn)/sigcn))
          IF (ccn2 /= 0D0) sumion = sumion - &
                                    ccn2*(erf(rr/ccncu2)/rr)**2*EXP(-rr/ccncul2)
          IF (ccn3 /= 0D0) sumion = sumion - &
                                    ccn3*(erf(rr/ccncu3)/rr)**3*EXP(-rr/ccncul3)
          IF (ccn4 /= 0D0) sumion = sumion - &
                                    ccn4*(erf(rr/ccncu4)/rr)**4*EXP(-rr/ccncul4)
          IF (ccn5 /= 0D0) sumion = sumion - &
                                    ccn5*(erf(rr/ccncu5)/rr)**5*EXP(-rr/ccncul5)
          IF (ccn6 /= 0D0) sumion = sumion - &
                                    ccn6*(erf(rr/ccncu6)/rr)**6*EXP(-rr/ccncul6)
          IF (ccn7 /= 0D0) sumion = sumion - &
                                    ccn7*(erf(rr/ccncu7)/rr)**7*EXP(-rr/ccncul7)
          IF (ccn8 /= 0D0) sumion = sumion - &
                                    ccn8*(erf(rr/ccncu8)/rr)**8*EXP(-rr/ccncul8)
          IF (ccn9 /= 0D0) sumion = sumion - &
                                    ccn9*(erf(rr/ccncu9)/rr)**9*EXP(-rr/ccncul9)
          IF (ccn10 /= 0D0) sumion = sumion - &
                                     ccn10*(erf(rr/ccncu10)/rr)**10*EXP(-rr/ccncul10)
          IF (ccnd /= 0D0) sumion = sumion - &
                                    ccnd*(erf(rr/ccncud)/rr)**dcn*EXP(-rr/ccnculd)

        ELSE IF (ityp1 == 2) THEN
          rr = (xe(ind1) - cx(ind2))**2 + (ye(ind1) - cy(ind2))**2 + &
               (ze(ind1) - cz(ind2))**2
        ELSE
          rr = (xk(ind1) - cx(ind2))**2 + (yk(ind1) - cy(ind2))**2 + &
               (zk(ind1) - cz(ind2))**2
          rr = MAX(small, SQRT(rr))
! WRITE(6,*) 'DEBUG: ',rr,ind1
          sumion = sumion + bkn/(ekn + EXP((rr - fkn)/sigkn))
          IF (bkn2 /= 0D0) sumion = sumion + bkn2/(ekn2 + EXP((rr - fkn2)/sigkn2))
          IF (ckn2 /= 0D0) sumion = sumion - &
                                    ckn2*(erf(rr/ckncu2)/rr)**2*EXP(-rr/ckncul2)
          IF (ckn3 /= 0D0) sumion = sumion - &
                                    ckn3*(erf(rr/ckncu3)/rr)**3*EXP(-rr/ckncul3)
          IF (ckn4 /= 0D0) sumion = sumion - &
                                    ckn4*(erf(rr/ckncu4)/rr)**4*EXP(-rr/ckncul4)
          IF (ckn5 /= 0D0) sumion = sumion - &
                                    ckn5*(erf(rr/ckncu5)/rr)**5*EXP(-rr/ckncul5)
          IF (ckn6 /= 0D0) sumion = sumion - &
                                    ckn6*(erf(rr/ckncu6)/rr)**6*EXP(-rr/ckncul6)
          IF (ckn7 /= 0D0) sumion = sumion - &
                                    ckn7*(erf(rr/ckncu7)/rr)**7*EXP(-rr/ckncul7)
          IF (ckn8 /= 0D0) sumion = sumion - &
                                    ckn8*(erf(rr/ckncu8)/rr)**8*EXP(-rr/ckncul8)
          IF (ckn9 /= 0D0) sumion = sumion - &
                                    ckn9*(erf(rr/ckncu9)/rr)**9*EXP(-rr/ckncul9)
          IF (ckn10 /= 0D0) sumion = sumion - &
                                     ckn10*(erf(rr/ckncu10)/rr)**10*EXP(-rr/ckncul10)
          IF (cknd /= 0D0) sumion = sumion - &
                                    cknd*(erf(rr/ckncud)/rr)**dkn*EXP(-rr/cknculd)
        END IF ! IF ityp1.eq.1

      ELSE ! Na-Na
! no short range interaction IN this CASE
      END IF

    ELSE IF (ityp2 == 5) THEN ! DFT

      IF (ityp1 <= 3) THEN ! DFT-GSM

        IF (ityp1 == 1) THEN ! core - DFT

        END IF

        IF (ityp1 == 3) THEN ! cation - DFT

        END IF

      END IF

    END IF

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(ityp1, ityp2) == 2) THEN ! Ar CASE

    IF (ityp1 == 1 .AND. ityp2 == 4) THEN

      dist2 = (cx(ind2) - xc(ind1))**2 + (cy(ind2) - yc(ind1))**2 + &
              (cz(ind2) - zc(ind1))**2
      dist = SQRT(dist2)
      dist = MAX(dist, small)

      sumion = varnashort(dist)

    ELSE IF (ityp1 == 1 .AND. ityp2 == 1) THEN

      dist2 = (xc(ind1) - xc(ind2))**2 + (yc(ind1) - yc(ind2))**2 + &
              (zc(ind1) - zc(ind2))**2
      dist = SQRT(dist2)
      dist = MAX(dist, small)

      IF (iararlj == 0) THEN
        sumion = vararshort(dist)
      ELSE
        sumion = vararlj(dist)
      END IF

    END IF

  ELSE ! other cases

! TO be implemented

  END IF

  rscaltmp = sumion

  RETURN
END SUBROUTINE getshortenergy
!------------------------------------------------------------

!************************************************************

SUBROUTINE addshortrepulsivepot(field, iswitch)
!************************************************************
  USE params
  IMPLICIT NONE
! gives each GSM particle its short range repulsive (Pauli)
! potential, i.e. short range interaction GSM-el
! right now ONLY the argon CASE is implemented

  REAL(DP), INTENT(IN OUT) :: field(kdfull2)
  INTEGER, INTENT(IN) :: iswitch

  INTEGER :: is
  REAL(DP) :: delinv

  REAL(DP), EXTERNAL :: varelcore, vfermi, funkfermi

  delinv = 1733D0/dx

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IF (isrtyp(1, 5) == 2) THEN ! The Argon CASE; c-DFT

    DO is = 1, nc

      IF (imobc(is) == iswitch) THEN

        CALL addtabtofield(field, varelcore0, xc(is), yc(is), zc(is), 1D0)
!old CALL addFuncToField(field,vArElCore,xc(is),yc(is),zc(is),1D0)

      END IF

    END DO !nc-loop

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(1, 5) == 1) THEN ! the MgO-CASE; c-DFT

    DO is = 1, nc

      IF (imobc(is) == iswitch) THEN

        CALL addfunctofield3(field, funkfermi, xc(is), yc(is), zc(is), 1D0, &
                             fermiac, fermibc, fermicc)

        CALL addfunctofield3(field, funkfermi, xc(is), yc(is), zc(is), 1D0, &
                             fermia2c, fermib2c, fermic2c)

      END IF

    END DO ! nc-loop

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(1, 5) == 3) THEN ! the general Fermi FUNCTION; c-DFT
    DO is = 1, nc

      IF (imobc(is) == iswitch) THEN

        CALL addfunctofield(field, vfermi, xc(is), yc(is), zc(is), 1D0)

      END IF

    END DO
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE ! c-DFT

! TO be implemented for other cases

  END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IF (isrtyp(2, 5) == 1) THEN ! the MgO-CASE, e-DFT

! nothing TO DO IN this model

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(2, 5) == 2) THEN ! the Ar-CASE, e-DFT

! nothing TO DO IN this model

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE ! e-DFT

! TO be implemented for other cases

  END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IF (isrtyp(3, 5) == 1) THEN ! the MgO-CASE, k-DFT

    DO is = 1, nk

      IF (imobk(is) == iswitch) THEN

        CALL addfunctofield3(field, funkfermi, xk(is), yk(is), zk(is), 1D0, &
                             fermiak, fermibk, fermick)

        CALL addfunctofield3(field, funkfermi, xk(is), yk(is), zk(is), 1D0, &
                             fermia2k, fermib2k, fermic2k)

      END IF

    END DO ! nk-loop

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(3, 5) == 2) THEN ! the Ar-CASE, k-DFT

! there are no cations IN this model

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(3, 5) == 3) THEN ! k-DFT by Fermi FUNCTION
    DO is = 1, nk

      IF (imobk(is) == iswitch) THEN
        CALL addfunctofield(field, vfermi, xk(is), yk(is), zk(is), 1D0)
      END IF
    END DO
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE ! k-DFT

! TO be implemented for other cases

  END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  RETURN
END SUBROUTINE addshortrepulsivepot
!------------------------------------------------------------

!************************************************************

SUBROUTINE addshortrepulsivepotonsubgrid(field, iswitch)
!************************************************************
  USE params
  IMPLICIT NONE
! gives each GSM particle its short range repulsive (Pauli)
! potential, i.e. short range interaction GSM-el
! right now ONLY the argon CASE is implemented

! nsgsize is the length of the subgrid IN one direction
! the total length of the subgrid box is 2*nsgsize+1

  REAL(DP), INTENT(IN OUT) :: field(kdfull2)
  INTEGER, INTENT(IN) :: iswitch

  INTEGER :: is
  REAL(DP) :: delinv

  REAL(DP), EXTERNAL :: varelcore, vfermi, funkfermi, funkpower

  delinv = 1733D0/dx

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IF (isrtyp(1, 5) == 1) THEN ! MgO-CASE ! c-DFT

    DO is = 1, nc

      IF (imobc(is) == iswitch) THEN

        CALL addfunctofieldonsubgrid3(field, funkfermi, &
                                      xc(is), yc(is), zc(is), 1D0, fermiac, fermibc, fermicc, 15)

        CALL addfunctofieldonsubgrid3(field, funkfermi, &
                                      xc(is), yc(is), zc(is), 1D0, fermia2c, fermib2c, fermic2c, 15)

        IF (ccel6 /= 0) THEN ! Hazardous TO compare a REAL TO an INTEGER...

          CALL addfunctofieldonsubgrid1(field, funkpower, &
                                        xc(is), yc(is), zc(is), ccel6, 6D0, 15)

        ELSE IF (ccel8 /= 0) THEN

          CALL addfunctofieldonsubgrid1(field, funkpower, &
                                        xc(is), yc(is), zc(is), ccel8, 8D0, 15)

        ELSE IF (ccel10 /= 0) THEN

          CALL addfunctofieldonsubgrid1(field, funkpower, &
                                        xc(is), yc(is), zc(is), ccel10, 10D0, 15)

        END IF

      END IF

    END DO ! nc-loop

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(1, 5) == 2) THEN ! The Argon CASE, c-DFT

    DO is = 1, nc

      IF (imobc(is) == iswitch) THEN ! choose ONLY cores which are
! mobile or immobile

        CALL addtabtofieldonsubgrid(field, varelcore0, xc(is), &
                                    yc(is), zc(is), 1D0, 10)

      END IF

    END DO ! nc-loop

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(1, 5) == 3) THEN ! c-DFT by Fermi FUNCTION

    CALL addfunctofieldonsubgrid(field, vfermi, xc(is), yc(is), zc(is), 1D0, 10)

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE ! c-DFT
! TO be implemented for other CASE

  END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IF (isrtyp(2, 5) == 1) THEN ! The MgO CASE, e-DFT
! nothing TO DO IN this model
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(2, 5) == 2) THEN ! The Ar CASE, e-DFT
! nothing TO DO IN this model
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE ! e-DFT
! TO be implemented for other cases
  END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  IF (isrtyp(3, 5) == 1) THEN ! The MgO CASE, k-DFT

    DO is = 1, nk

      IF (imobk(is) == iswitch) THEN
        CALL addfunctofieldonsubgrid3(field, funkfermi, xk(is) &
                                      , yk(is), zk(is), 1D0, fermiak, fermibk, fermick, 15)

        CALL addfunctofieldonsubgrid3(field, funkfermi, xk(is) &
                                      , yk(is), zk(is), 1D0, fermia2k, fermib2k, fermic2k, 15)

        IF (ckel6 /= 0) THEN

          CALL addfunctofieldonsubgrid1(field, funkpower, &
                                        xk(is), yk(is), zk(is), ckel6, 6D0, 15)

        ELSE IF (ckel8 /= 0) THEN

          CALL addfunctofieldonsubgrid1(field, funkpower, &
                                        xk(is), yk(is), zk(is), ckel8, 8D0, 15)

        ELSE IF (ckel10 /= 0) THEN

          CALL addfunctofieldonsubgrid1(field, funkpower, &
                                        xk(is), yk(is), zk(is), ckel10, 10D0, 15)

        END IF

      END IF

! CALL addFuncToFieldOnSubgrid(potion,vArElCore,xk(is),
! & yk(is),zk(is),1D0,10)
    END DO

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(3, 5) == 2) THEN ! The Ar CASE, k-DFT
! nothing TO DO IN this model, no cations present
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE IF (isrtyp(3, 5) == 3) THEN ! k-DFT, by Fermi FUNCTION

    CALL addfunctofieldonsubgrid(field, vfermi, xk(is), yk(is), zk(is), 1D0, 10)

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ELSE ! k-DFT
! TO be implemented for other cases
  END IF
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! STOP

  RETURN
END SUBROUTINE addshortrepulsivepotonsubgrid
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE setbornparas(ityp1, ityp2)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ityp1
  INTEGER, INTENT(IN) :: ityp2

  IF (ityp1 == 1) THEN ! 1:c
    IF (ityp2 == 1) THEN ! 2:c
      rvectmp(1) = bcc
      rvectmp(2) = sigcc
      rvectmp(3) = ccc
    ELSE IF (ityp2 == 2) THEN ! 2:e
      rvectmp(1) = bcv
      rvectmp(2) = sigcv
      rvectmp(3) = ccv
    ELSE IF (ityp2 == 3) THEN ! 2:k
      rvectmp(1) = bck
      rvectmp(2) = sigck
      rvectmp(3) = cck
    END IF
  ELSE IF (ityp1 == 2) THEN ! 1:e
    IF (ityp2 == 2) THEN ! 2:e
      rvectmp(1) = bvv
      rvectmp(2) = sigvv
      rvectmp(3) = cvv
    ELSE IF (ityp2 == 3) THEN ! 2:k
      rvectmp(1) = bkv
      rvectmp(2) = sigkv
      rvectmp(3) = ckv
    END IF
  ELSE IF (ityp1 == 3) THEN ! 1:k
    rvectmp(1) = bkk
    rvectmp(2) = sigkk
    rvectmp(3) = ckk
  END IF

  RETURN
END SUBROUTINE setbornparas
!------------------------------------------------------------
#endif

!------------------------------------------------------------

REAL(DP) FUNCTION funkpower(r, p)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) ::r
  REAL(DP), INTENT(IN) ::p

  funkpower = erf(r/2D0)/(r**p)

  RETURN
END FUNCTION funkpower
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION funkfermi(r, a, b, c)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: a
  REAL(DP), INTENT(IN) :: b
  REAL(DP), INTENT(IN) :: c

  funkfermi = a/(1D0 + EXP(b*(r - c)))

  RETURN
END FUNCTION funkfermi
!------------------------------------------------------------

!------------------------------------------------------------

REAL(DP) FUNCTION funkderfermi(r, a, b, c)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: r
  REAL(DP), INTENT(IN) :: a
  REAL(DP), INTENT(IN) :: b
  REAL(DP), INTENT(IN) :: c

  funkderfermi = -a*b*EXP(b*(r - c))/(1D0 + EXP(b*(r - c)))**2D0

  RETURN
END FUNCTION funkderfermi
!------------------------------------------------------------
!------------------------------------------------------------

SUBROUTINE checkstability(rho)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: rho(kdfull2*2)

  INTEGER :: i
  REAL(DP) :: dmaxfxc, dmaxfyc, dmaxfzc, dmaxfxe, dmaxfye, dmaxfze, dmaxfxk, dmaxfyk, dmaxfzk, dmaxfxn, dmaxfyn, dmaxfzn

  COMPLEX(DP) :: psidummy(1) !! bugfix for force

! the 'psidummy' trick works ONLY for local Gaussian PsP
  IF (ipsptyp == 1) STOP ' CHECKSTABILITY must not be used with Goedecker PsP'

  CALL getforces(rho, psidummy, -1, 0)

  dmaxfxc = 0D0
  dmaxfyc = 0D0
  dmaxfzc = 0D0
  dmaxfxe = 0D0
  dmaxfye = 0D0
  dmaxfze = 0D0
  dmaxfxk = 0D0
  dmaxfyk = 0D0
  dmaxfzk = 0D0
  dmaxfxn = 0D0
  dmaxfyn = 0D0
  dmaxfzn = 0D0

#if(raregas)
  DO i = 1, nc
    IF (imobc(i) == 1) THEN
      IF (ABS(fxc(i)) > ABS(dmaxfxc)) dmaxfxc = fxc(i)
      IF (ABS(fyc(i)) > ABS(dmaxfyc)) dmaxfyc = fyc(i)
      IF (ABS(fzc(i)) > ABS(dmaxfzc)) dmaxfzc = fzc(i)
    END IF
  END DO

  DO i = 1, NE
    IF (imobe(i) == 1) THEN
      IF (ABS(fxe(i)) > ABS(dmaxfxe)) dmaxfxe = fxe(i)
      IF (ABS(fye(i)) > ABS(dmaxfye)) dmaxfye = fye(i)
      IF (ABS(fze(i)) > ABS(dmaxfze)) dmaxfze = fze(i)
    END IF
  END DO

  DO i = 1, nk
    IF (imobk(i) == 1) THEN
      IF (ABS(fxk(i)) > ABS(dmaxfxk)) dmaxfxk = fxk(i)
      IF (ABS(fyk(i)) > ABS(dmaxfyk)) dmaxfyk = fyk(i)
      IF (ABS(fzk(i)) > ABS(dmaxfzk)) dmaxfzk = fzk(i)
    END IF
  END DO
#endif

  DO i = 1, nion
    IF (ABS(fx(i)) > ABS(dmaxfxn)) dmaxfxn = fx(i)
    IF (ABS(fy(i)) > ABS(dmaxfyn)) dmaxfyn = fy(i)
    IF (ABS(fz(i)) > ABS(dmaxfzn)) dmaxfzn = fz(i)
  END DO

  WRITE (6, *) '##############################################'
#if(raregas)
  WRITE (6, '(a,3e17.7)') 'max fc: ', dmaxfxc, dmaxfyc, dmaxfzc
  WRITE (6, '(a,3e17.7)') 'max fe: ', dmaxfxe, dmaxfye, dmaxfze
  WRITE (6, '(a,3e17.7)') 'max fk: ', dmaxfxk, dmaxfyk, dmaxfzk
#endif
  WRITE (6, '(a,3e17.7)') 'max fN: ', dmaxfxn, dmaxfyn, dmaxfzn
  WRITE (6, *) '##############################################'

  RETURN
END SUBROUTINE checkstability
!------------------------------------------------------------

!------------------------------------------------------------

SUBROUTINE checkespot(x, y, iunit)
!------------------------------------------------------------
  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::iunit
  REAL(DP), INTENT(IN) ::x
  REAL(DP), INTENT(IN) ::y

#if(raregas)

  INTEGER :: i, iz
  REAL(DP) :: potex, z
  REAL(DP) :: rr, xr, yr, zr
  REAL(DP), EXTERNAL :: v_soft

  DO iz = 1, 1000
    z = -10D0*0.1D0

    potex = 0D0

    DO i = 1, nc + NE + nk
      CALL getparas(i)

      xr = x - rvectmp(1)
      yr = y - rvectmp(2)
      zr = z - rvectmp(3)

! WRITE(6,*) sigTmp,chgTmp

      rr = SQRT(xr*xr + yr*yr + zr*zr)
      rr = MAX(rr, small)

      potex = e2*chgtmp*v_soft(rr, sigtmp) + potex

    END DO

    WRITE (iunit, '(1f15.4,1e17.7)') z, potex

  END DO
#endif

  RETURN
END SUBROUTINE checkespot
!------------------------------------------------------------
