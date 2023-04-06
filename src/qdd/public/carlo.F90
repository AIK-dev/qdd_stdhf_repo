
SUBROUTINE init_simann()

! Initializes variables for simulated annealing.

  USE params
  IMPLICIT NONE

  INTEGER :: ifall
  LOGICAL, PARAMETER :: tMCparamsin = .FALSE.

! This is the master routine for the Metropolis algorithm,
! see M.P. Allen, D.J. Tildesley, Computer Simulation of Liquids,
! Oxford University Press, New York 1987.

! The Monte-Carlo (MC) scheme is run IN two layers. IN both layers, ionic
! positions are shifted randomly and are tracked by randomized path
! twoard lowering energy.
!
! The outer loop applies the MC minimization scheme using the total energies
! from the full KS solution and calls another search IN the inner loop.
!
! The inner loop applies the MC minimization scheme computing energies
! with fixed KS wavefunctions using the Hellmann-Feynman theorem.
!
! IN both loops, the temperature (used IN MC decision) and stepsize
! for ionic shifts is reduced gradually during the process. Ideally,
! one ends up with zero shifst and temperature zero.
!
! There are two ways TO enter the MC driver PARAMETER:
! 1. One sets IN he above statement 'tMCparamsin=.TRUE.' and reads
! the parameters from the FILE 'ionen-IN.<NAME>' ('NAME' is the
! FILE qualifie from the 'for005.<NAME>). A typical FILE looks
! like that
! 9,0,1.5 ionsin,iknow,facann
! 10,50,100 nrun,nloop1,nloop2
! 0.02,0.05,6.E-5,2.0e-3 cptemp,delpos,ERR,errks0
! 121466 ifall
! 0.01,0.01,0.001,10 trfac2,prfac2,errsim,ncsim
! 0.002,8 errtot,ncon
! 0.0,0.01,0.01 erfac1,trfac1,prfac1
! (of course without the initial "!").
! 2. One sets 'tMCparamsin=.FALSE.' and defines the parameters
! IN the initialization below.

! The parameters driving the MC scheme are initialized here:
!
  ionsin = 3 ! obsolete
  ifall = 121466 ! obsolete
  facann = 1.5D0 ! obsolete
  iknow = 0 ! indicator for a priori knowledge of ionic config.
  nrun = 100 ! NUMBER of full runs of outer and inner loop
  nloop1 = 100 ! NUMBER of outer runs (with KS optimization)
  nloop2 = 200 ! NUMBER of inner runs (fixed KS wavefunctions)
  cptemp = 0.008D0 ! initial temperature T
  delpos = 0.2D0 ! inital SIZE of trial shift of ions
  errtot = 0.002D0 ! terminator outer loop: change IN E must be 'ncon' times below 'errtot'
  ncon = 8 ! required NUMBER of steps below 'errtot' IN outer loop
  trfac1 = 0.08D0 ! reduction 'cptemp' IN outer step (Metropolis)
  prfac1 = 0.05D0 ! reduction 'delpos' IN outer step (Metropolis)
  trfac2 = 0.03D0 ! reduction 'cptemp' IN inner step (Metropolis)
  prfac2 = 0.05D0 ! reduction 'delpos' IN inner step
  errsim = 0.001D0 ! terminator inner loop: change IN E must be 'ncsim' times below 'errsim'
  ncsim = 10 ! required NUMBER of steps below 'errsim' IN inner loop
! The following parameters are acually not used:
  ERR = 6.0D-5 ! lower limit of accuracy IN outer step
  errks0 = 5.0D-3 ! initial error on KS
  erfac1 = 0.1D0 ! reduction factor for 'errKS'

! Initialize random seed.
  CALL RANDOM_SEED()
  CALL RANDOM_SEED(SIZE=isize_seed)
  ALLOCATE (rand_seed(isize_seed))

  IF (myn == 0) THEN

    IF (tMCparamsin) &
      OPEN (UNIT=44, STATUS='old', FORM='FORMATTED', FILE='ionen-IN.'//outnam)
    WRITE (6, '(/a)') 'PARAMETERS FOR METROPOLIS:'
    IF (tMCparamsin) READ (44, *) ionsin, iknow, facann
    WRITE (6, '(a,2i2,f7.3)') ' ionsin,iknow,facann: ', ionsin, iknow, facann
    IF (tMCparamsin) READ (44, *) nrun, nloop1, nloop2
    WRITE (6, '(a,3i4)') ' nrun,nloop1,nloop2: ', nrun, nloop1, nloop2
    IF (tMCparamsin) READ (44, *) cptemp, delpos, ERR, errks0
    WRITE (6, '(a,4f12.6)') ' cptemp,delpos,ERR,errks0: ', cptemp &
      , delpos, ERR, errks0
    IF (tMCparamsin) READ (44, *) ifall
    WRITE (6, '(a,i10)') ' ifall (argument no longer used):', ifall
    IF (tMCparamsin) READ (44, *) trfac2, prfac2, errsim, ncsim
    WRITE (6, '(a,3f8.3,i3)') ' trfac2,prfac2,errsim,ncsim: ' &
      , trfac2, prfac2, errsim, ncsim
    IF (tMCparamsin) READ (44, *) errtot, ncon
    WRITE (6, '(a,f8.3,i3)') ' errtot,ncon: ', errtot, ncon
    IF (tMCparamsin) READ (44, *) erfac1, trfac1, prfac1
    WRITE (6, '(a,3f8.3)') ' erfac1,trfac1,prfac1: ', erfac1, trfac1, prfac1
    WRITE (6, *) ' '
    WRITE (6, *) ' '
    IF (tMCparamsin) CLOSE (44)

! Take random seed from first node.
    CALL RANDOM_SEED(GET=rand_seed)

  END IF


! Start all nodes with the same random seed.
  CALL RANDOM_SEED(PUT=rand_seed)

END SUBROUTINE init_simann

! *************************************

SUBROUTINE simann(psir, rho, aloc)

! Optimization of ionic positions by simulated annealing.
!
! Input/Output:
! psir = set of s.p. wavefunctions
! rho = local densities (spin-up and spin-down)
! aloc = local KS potentials

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: psir(kdfull2, kstate)
  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)

  INTEGER :: ion, jrun
  REAL(DP) :: cptem0, delbin, delps0, ebold, hph2
  REAL(DP), EXTERNAL :: energ_ions

! SAVE starting values of sim. annealing for further use:
  delps0 = delpos
  cptem0 = cptemp

! loop over various runs of the metropolis part:

  IF (tstat) THEN
    CALL resume(psir, outnam)
  END IF
  DO jrun = 1, nrun

! reset annealing schedule for further sim. annealing:

    cptemp = cptem0
    delpos = delps0
! nloop2 = nloop2/(2.0*facann)
    errks = errks0
    nyes = 0 !counter for convergence IN main loop
    ebold = 1D20
    hph2 = 1D0

    metropolis: DO loop1 = 1, nloop1

! reach convergence in Kohn-Sham loop

      CALL calcpseudo()

      ecorr = energ_ions()

! calculate the projectors around the new ionic positions

      IF (ipsptyp == 1) THEN
        DO ion = 1, nion
          CALL calc_proj(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
          WRITE (6, '(a,i2,a,3f9.4)') 'ion', ion, '=', cx(ion), cy(ion), cz(ion)
        END DO
        CALL checkproj(ion)
      END IF
      ismax = 200
      CALL statit(psir, rho, aloc)


      delbin = ABS(ebold - binerg)
      IF (loop1 >= 2) THEN
        WRITE (6, '(/a,f8.5,a,f8.5)') 'new binding energy =', binerg, &
          ' energy difference with last loop =', binerg - ebold
      END IF

      IF (delbin < errtot) THEN
        nyes = nyes + 1
        IF (nyes == ncon) EXIT metropolis
      ELSE
        nyes = 0
      END IF
      WRITE (6, '(a,f9.4,i3)') 'deltaE, nyes =', binerg - ebold, nyes
      ebold = binerg

! if convergence is not achieved, reduce annealing schedule
! and Kohn-Sham accuracy 'errks':

      errks = errks*(1D0 - erfac1)
      IF (errks < ERR) errks = ERR !maximum accuracy
      cptemp = cptemp*(1D0 - trfac1)
      delpos = delpos*(1D0 - prfac1)

      WRITE (6, '(a)') '-->NOW STARTING TO OPTIMIZE THE IONIC POSITIONS keep waiting ...'

      CALL minpos(rho, psir)

      CALL cenmass()

      OPEN (27, POSITION='append', FILE='pminpos.'//outnam)
      DO ion = 1, nion
        WRITE (27, '(i3,4f13.5)') loop1, cx(ion), cy(ion), cz(ion), binerg
      END DO
      FLUSH (27)

      WRITE (6, '(A/)') '-->ROLLING THE DICE COMPLETED'
    END DO metropolis

! if you come here, energetic accuracy is finally achieved

    WRITE (6, '(A)') ' NORMAL TERMINATION IN METROPOLIS'
  END DO
  CALL rsave(psir, 0, outnam)

  RETURN
END SUBROUTINE simann
! *************************************

SUBROUTINE minpos(rho, psimc)

! Change ionic positions according TO Monte-Carlo estimator

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(kdfull2)
  REAL(DP), INTENT(IN) :: psimc(kdfull2, kstate)

  REAL(DP), ALLOCATABLE :: q1(:)
  REAL(DP), ALLOCATABLE :: rhoion(:)

  INTEGER :: i, ii, ion, ionvar, iknowl, j, loop2, loop3
  INTEGER :: nb, nbut, nmin, nup, ndown
  REAL(DP) :: delps3, dii, dist, diffen, difii, diftot, deltax, deltay, deltaz, delta2
  REAL(DP) :: eold, eionic, facmo, oh, radone, sumnl, teml3
  REAL(DP) :: gxr, gyr, gzr, gxl, gyl, gzl
  REAL(DP) :: xred2, yred2, zred2, xled2, yled2, zled2, xion, yion, zion, xvar, yvar, zvar
  REAL(DP) :: rand0, rand3(3)

  LOGICAL :: ans

  ALLOCATE (q1(kdfull2), rhoion(kdfull2))

  eloctot = 0D0
  enoloctot = 0D0
  eiontot = 0D0
  DO ion = 1, nion

    xion = cx(ion)
    yion = cy(ion)
    zion = cz(ion)

    CALL rhopsg(xion, yion, zion, rhoion, ion)
    eloc(ion) = 0D0
    DO i = 1, nxyz
      eloc(ion) = eloc(ion) + rho(i)*rhoion(i)
    END DO
    eloc(ion) = eloc(ion)*dvol
    eloctot = eloctot + eloc(ion)

    IF (ipsptyp == 1) THEN
      CALL calc_proj(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
      enoloc(ion) = 0D0
      DO nb = 1, nstate
        CALL nonlocalr(psimc(1, nb), q1)
        sumnl = 0D0
        DO i = 1, ifin(ion)
          ii = icount(i, ion)
          sumnl = sumnl + psimc(ii, nb)*q1(ii)
        END DO
        enoloc(ion) = enoloc(ion) + sumnl*occup(nb)*dvol
      END DO
      enoloctot = enoloctot + enoloc(ion)
    END IF

    DO j = ion + 1, nion
      dist = SQRT((cx(j) - xion)**2 + (cy(j) - yion)**2 + (cz(j) - zion)**2)
      eion(ion, j) = e2*ch(np(ion))**2/dist
      eiontot = eiontot + eion(ion, j)
    END DO
  END DO

  etot = -eloctot + enoloctot + eiontot
  eold = etot

  delps3 = delpos
  teml3 = cptemp
  nmin = 0 ! counts loops with energy change < errsim

! NOW THE FUN BEGINS
! the outer dice loop:

  diftot = 0D0 ! total energy gain through moving ions
  difii = 0D0 ! ionic energy gain through moving ions

  DO loop2 = 1, nloop2

    ndown = 0 ! steps downwards
    nbut = 0 ! steps upwards although energy gets worse
    nup = 0 ! steps without any success
    oh = 0.5D0

! the inner loop:

    DO loop3 = 1, nion

! roll the dice TO decide which ion is TO be moved:

      CALL RANDOM_NUMBER(rand0)
      ionvar = INT(nion*rand0) + 1

! vary the coordinates of this ion:

      infinite: DO ! Loop until the ion as been displaced IN a valid location.
        CALL RANDOM_NUMBER(rand3)
        deltax = 2D0*delps3*(rand3(1) - oh)
        deltay = 2D0*delps3*(rand3(2) - oh)
        deltaz = 2D0*delps3*(rand3(3) - oh)
        delta2 = deltax*deltax + deltay*deltay + deltaz*deltaz
        IF (delta2 > delps3*delps3) CYCLE infinite

        xvar = cx(ionvar) + deltax
        yvar = cy(ionvar) + deltay
        zvar = cz(ionvar) + deltaz

        ! Ion within sensible grid area?

        radone = 3D0
        facmo = 2D0

        xred2 = xvar + facmo*radone
        yred2 = yvar + facmo*radone
        zred2 = zvar + facmo*radone
        xled2 = xvar - facmo*radone
        yled2 = yvar - facmo*radone
        zled2 = zvar - facmo*radone

        gxr = REAL(nx, DP)*dx
        gyr = REAL(ny, DP)*dy
        gzr = REAL(nz, DP)*dz
        gxl = -gxr
        gyl = -gyr
        gzl = -gzr
        IF (((gxr > xred2) .AND. (gxl < xled2)) .AND. &
        &((gyr > yred2) .AND. (gyl < yled2)) .AND. &
        &((gzr > zred2) .AND. (gzl < zled2))) EXIT infinite
      END DO infinite
! calculate new one-ion-electron energy of ion no. 'ionvar'

      CALL rhopsg(xvar, yvar, zvar, rhoion, ionvar)
      eloc(0) = 0D0
      DO i = 1, nxyz
        eloc(0) = eloc(0) + rho(i)*rhoion(i)
      END DO
      eloc(0) = eloc(0)*dvol

      IF (ipsptyp == 1) THEN
        CALL calc_proj(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
        enoloc(0) = 0D0
        DO nb = 1, nstate
          CALL nonlocalr(psimc(1, nb), q1)
          sumnl = 0D0
          DO i = 1, ifin(0)
            ii = icount(i, 0)
            sumnl = sumnl + psimc(ii, nb)*q1(ii)
          END DO
          enoloc(0) = enoloc(0) + sumnl*occup(nb)*dvol
        END DO
      END IF

! compute energy-difference diffen=etvar-etold caused by the move:

      CALL move(ionvar, xvar, yvar, zvar, diffen, dii)

! The mysterious decision: go down or CALL the ORACLE

      ans = (diffen < 0D0)
      IF (ans) THEN
        ndown = ndown + 1 !good CALL
      ELSE
        iknowl = iknow
        CALL metrop(diffen, teml3, iknowl, ans)
        IF (ans) nbut = nbut + 1
      END IF
      IF (ans) THEN
        diftot = diftot + diffen
        difii = difii + dii
      END IF

! IF 'yes', swap coordinates and energies

      IF (ans) THEN
        cx(ionvar) = xvar
        cy(ionvar) = yvar
        cz(ionvar) = zvar
        eloc(ionvar) = eloc(0)
        enoloc(ionvar) = enoloc(0)
        DO ion = 1, ionvar - 1
          eion(ion, ionvar) = eiinew(ion)
        END DO
        DO ion = ionvar + 1, nion
          eion(ionvar, ion) = eiinew(ion)
        END DO
      ELSE
        nup = nup + 1 !bad CALL
      END IF
    END DO ! END of dice loop at given T

! NEW ONE-ION ENERGY AFTER ORACLE:

    eloctot = 0D0 !reset energies for new calculation
    eiontot = 0D0
    enoloctot = 0D0

    DO i = 1, nion
      eloctot = eloctot + eloc(i)
      enoloctot = enoloctot + enoloc(i)
      DO j = i + 1, nion
        eiontot = eiontot + eion(i, j)
      END DO
    END DO
    eionic = -eloctot + enoloctot + eiontot

! lower annealing parameters for new dice loop:

    delps3 = delps3*(1D0 - prfac2)
    teml3 = teml3*(1D0 - trfac2)

! Criteria for convergence: END of loop IF ionic energy has been
! changed for 'ncsim' consecutive times by less than 'errsim':

    IF (ABS(eold - eionic) < errsim) THEN
      nmin = nmin + 1
    ELSE
      nmin = 0
    END IF
    IF (nmin == ncsim) EXIT

    eold = eionic
  END DO !END of annealing

! history of dice loops:

  WRITE (6, '(/a,i3)') 'no. of dice loops: ', loop2
  WRITE (6, '(a,f12.4)') 'total energy gain= ', diftot
  WRITE (6, '(a,f12.4)') 'ionic energy gain= ', difii
  WRITE (6, '(3x,a,/,5x,a,2f8.5,/,5x,a,3g15.6)') 'Leaving minpos with:', &
    'cptemp, delpos = ', teml3, delps3, 'eloctot,enoloctot,eiontot= ', &
    eloctot, enoloctot, eiontot
  WRITE (6, '(5x,a,f12.4)') 'eionic= ', eionic

  DEALLOCATE (q1, rhoion)

  RETURN
END SUBROUTINE minpos

! *************************************************

SUBROUTINE move(ionvar, xvar, yvar, zvar, diffen, dii)

! Computes energy-difference 'diffen' resulting from moving 'ion'
! to (xvar,yvar,zvar).
! Prepare energies on fieldindex '0' for a possible swapping
! of coordinates after oracle.
! Calculate new ion-ion-energy of ion no. 'ionvar' with
! all the other ions.

  USE params
  IMPLICIT NONE

  INTEGER :: ion, ionvar
  REAL(DP) :: dieloc, dienl, diffen, dii, dist, eiivar, eionold, etold, etvar, xvar, yvar, zvar

  eiivar = 0D0
  DO ion = 1, nion
    IF (ion /= ionvar) THEN
      dist = SQRT((cx(ion) - xvar)**2 + (cy(ion) - yvar)**2 + (cz(ion) - zvar)**2)
      eiinew(ion) = e2*ch(np(ion))**2/dist
      eiivar = eiivar + eiinew(ion)
    END IF
  END DO

! calculate old ion-ion-energy of ion no. 'ionvar' with
! all the other ions

  eionold = 0D0
  DO ion = 1, ionvar - 1
    eionold = eionold + eion(ion, ionvar)
  END DO
  DO ion = ionvar + 1, nion
    eionold = eionold + eion(ionvar, ion)
  END DO

! compare the energies for the oracle

  etold = eionold - eloc(ionvar) + enoloc(ionvar)
  etvar = eiivar - eloc(0) + enoloc(0)
  diffen = etvar - etold
  dii = eiivar - eionold
  dieloc = eloc(0) - eloc(ionvar)
  dienl = enoloc(0) - enoloc(ionvar)

  RETURN
END SUBROUTINE move

! *******************

REAL(DP) FUNCTION ran0(idum)
  USE params, ONLY: DP
  IMPLICIT NONE

! *******************

! Numerical recipes function ran0 (p. 270):
! Creates a random number from input 'idum'

  INTEGER, INTENT(OUT) :: idum

!REAL(8):: ran0 !danger, DO not declare ran0 as REAL*8
  INTEGER, PARAMETER :: ia = 16807
  INTEGER, PARAMETER :: im = 2147483647
  REAL(DP), PARAMETER :: am = 1D0/im
  INTEGER, PARAMETER :: iq = 127773
  INTEGER, PARAMETER :: ir = 2836
  INTEGER, PARAMETER :: mask = 123459876
  INTEGER :: k

  idum = IEOR(idum, mask)
  k = idum/iq
  idum = ia*(idum - k*iq) - ir*k
  IF (idum < 0) idum = idum + im
  ran0 = am*idum
  idum = IEOR(idum, mask)

  RETURN
END FUNCTION ran0

! **************************************

SUBROUTINE metrop(diffen, t, iknowi, ans)

! The Metropolis algorithm

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: diffen
  REAL(DP), INTENT(IN) :: t
  INTEGER, INTENT(IN) :: iknowi
  LOGICAL, INTENT(OUT) :: ans

  REAL(DP) :: bolfac, expon, rand0

  IF (t > 0D0) THEN
    expon = diffen/t
    IF (expon < 80D0) THEN
      bolfac = EXP(-expon)
    ELSE
      bolfac = 0D0
    END IF
  ELSE
    bolfac = 0D0
  END IF

  IF (iknowi == 0) THEN
    CALL RANDOM_NUMBER(rand0)
    ans = (rand0 < bolfac)

! proper choice whenever the user has not the
! faintest idea about the shape of the cluster!
! -> provides the possibility to escape a local
! minimum in order to reach the global one.

  ELSE IF (iknowi == 1) THEN
    IF (loop1 == 0) THEN
      ans = (0.95D0 < bolfac) ! provides good starting point
    ELSE IF (loop1 > 0) THEN
      ans = ((1.8D0 - 0.955D0**(LOG(loop1 + 0.00001D0))) < bolfac)

! this choice cares for a fast iteration ->
! DANGER: USE ONLY IF starting point is chosen well,
! meaning the starting configuration is not too
! far away from the desired final result.
! If this is not the case, being trapped in a local
! minimum is unavoidable!!

    END IF
  END IF

  RETURN
END SUBROUTINE metrop

! ********************

SUBROUTINE cenmass()

! Move ions globally to shift center-of-mass to origin.

  USE params
  IMPLICIT NONE

  INTEGER :: icoo, ion

  REAL(DP) :: gamu
  REAL(DP) :: cenmas(3)
  REAL(DP) :: pos(ng, 3)

  DO ion = 1, nion
    pos(ion, 1) = cx(ion)
    pos(ion, 2) = cy(ion)
    pos(ion, 3) = cz(ion)
  END DO

! compute center of mass

  cenmas(1) = zero
  cenmas(2) = zero
  cenmas(3) = zero
  gamu = zero

  DO ion = 1, nion
    DO icoo = 1, 3
      cenmas(icoo) = cenmas(icoo) + amu(np(ion))*pos(ion, icoo)
    END DO
    gamu = gamu + amu(np(ion))
  END DO

  DO icoo = 1, 3
    cenmas(icoo) = cenmas(icoo)/gamu
  END DO

! transform TO center of mass coordinates

  DO ion = 1, nion
    DO icoo = 1, 3
      pos(ion, icoo) = pos(ion, icoo) - cenmas(icoo)
    END DO
  END DO

  DO ion = 1, nion
    cx(ion) = pos(ion, 1)
    cy(ion) = pos(ion, 2)
    cz(ion) = pos(ion, 3)
  END DO

  RETURN
END SUBROUTINE cenmass
