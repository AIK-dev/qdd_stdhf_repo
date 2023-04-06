
! ---scan_coll-----------------------------------------------------------

SUBROUTINE scan_coll(psir)
  USE params
!USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)

! Scans the total excitation probabilities after collision with very fast
! charged ion.
! The static wavefunctions is handed through IN 'psdir' via list.
! The switching parameters are handled through COMMON, they are:
! bcol1 = first collision PARAMETER
! bcol2 = last collision PARAMETER
! dbcol = step IN collision PARAMETER
! ntheta = NUMBER of steps IN angle theta
! 0 --> no computation
! 1 --> ONLY angle 0
! n>1 --> n bins between 0 and PI
! nphi = NUMBER of steps IN angle phi
! 0 --> no computation
! 1 --> ONLY angle 0
! n>1 --> n bins between 0 and 2*PI
! betacol = v/c for colliding ion
! chgcol = charge for colliding ion
! Results are printed IN FILE 'collision.res'

  REAL, INTENT(IN OUT) :: psir(kdfull2, kstate)

!-------------------------------------------------------------------------

! nothing happens no angle is asked for

  IF (ntheta == 0 .OR. nphi == 0) RETURN

  OPEN (8, STATUS='unknown', FORM='FORMATTED', FILE='collision.'//outnam)
  WRITE (8, '(a,2(1pg13.5)/a)') &
    '# total excitation probabilies for colliding ion, v/c,charge=', &
    betacol, chgcol, '# b theta phi probability '
  nbcol = (bcol2 - bcol1)/dbcol + 1
  dbcol = (bcol2 - bcol1)/MAX(nbcol - 1, 1)
  dtheta = pi/MAX(ntheta - 1, 1)
  dphi = pi/MAX(nphi - 1, 1)
  DO itheta = 1, ntheta
    thetacol = (itheta - 1)*dtheta
    DO iphi = 1, nphi
      phicol = (iphi - 1)*dphi
      DO ibcol = 1, nbcol
        bcol = bcol1 + (ibcol - 1)*dbcol
        pr = coll_total(psir, bcol, betacol, thetacol, phicol, chgcol)
        WRITE (8, '(1x,3f8.5,1pg13.5)') bcol, thetacol, phicol, pr
      END DO
    END DO
  END DO
  CLOSE (8)

  RETURN
END SUBROUTINE scan_coll
! ---coll_total-----------------------------------------------------------

FUNCTION coll_total(psir, bcol, betcol, thetacol, phicol, chargecol)

! Computes total excitation probability for collision with a very fast
! charged ion. The formula for small probability is used which employs
! just the variance of the phase OPERATOR.
! Input through list is:
! psir = REAL array containing all occupied s.p. wavefunctions
! bcol = collision PARAMETER (distance of closest approach)
! betcol = velocity/c of colliding ion
! thetacol = angle of collision trajectory with respect TO x-axis
! phicol = angle of collision trajectory IN y-z plane
! chargecol = charge of colliding ion
! Input through COMMON is:
! occup = array of occupation numbers
! ispin = array of spin orientations
! The PARAMETER 'damp' for the Coulomb cutoff is set at compile time
! within the header of this routine.
! The routine uses REAL wavefunctions and should be called after
! static iteration.

! !! The collision angles are not yet activated.
! !! The routine assumes presently a collision trajectory parallel
! !! TO the x-axis.

  USE params
!USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)

  REAL, INTENT(IN) :: psir(kdfull2, kstate)
  REAL, INTENT(IN) :: bcol
  REAL, INTENT(IN) :: betcol
  REAL, INTENT(IN OUT) :: thetacol
  REAL, INTENT(IN OUT) :: phicol
  REAL, INTENT(IN) :: chargecol
  REAL, PARAMETER :: damp = 1.0

  LOGICAL, PARAMETER :: ttest = .false.

!-------------------------------------------------------------------------

! check workspace
! w1 is used for intermediate s.p. state
! w2 is used for phase field

  IF (usew1) STOP ' IN COLL_TOTAL: workspace W1 already active '
  IF (usew2) STOP ' IN COLL_TOTAL: workspace W2 already active '
  usew1 = .true.
  usew2 = .true.

! prepare phase field

  damp2 = damp**2
  denom = 1.0/(damp2 + bcol**2)
  colfac = chargecol*e2/betcol/274.0
  ind = 0
  DO jz = 1, nz2
    z = (jz - nz)*dz
    zb2 = (z - bcol)**2
    DO jy = 1, ny2
      y = (jy - ny)*dy
      y2 = y*y
      DO jx = 1, nx2
! x = (jx-NX)*DX
        ind = ind + 1
        w2(ind) = LOG(denom*(zb2 + y2 + damp2))*colfac
      END DO
    END DO
  END DO

! loop over 'alpha' states

  sumdir = 0.0
  sumex = 0.0
  DO na = 1, nstate
    ispa = ispin(nrel2abs(na))
    occa = occup(na)

! multiply phase profile and accumulate variance

    acc = 0.0
    DO ind = 1, kdfull2
      w1(ind) = w2(ind)*psir(ind, na)
      acc = w1(ind)*w1(ind) + acc
    END DO
    sumdir = occa*acc*dvol + sumdir
    IF (ttest) WRITE (8, '(a,i5,2(1pg13.5))') &
      ' DIRECT: na,occ,acc=', na, occa, acc*dvol

! subtract contributions from occupied intermediate states

    DO nb = 1, nstate
      ispb = ispin(nrel2abs(nb))
      occb = occup(nb)
      acc = 0.0
      DO ind = 1, kdfull2
        acc = w1(ind)*psir(ind, nb) + acc
      END DO
      sumex = occa*occb*(acc*dvol)**2 + sumex
      IF (ttest) WRITE (8, '(a,2i5,3(1pg13.5))') &
        ' exchange: na,nb,occa,occb,acc=', na, nb, occa, occb, (acc*dvol)**2
    END DO

  END DO

! finally compose RESULT

  IF (numspin == 2) THEN
    coll_total = sumdir - sumex
  ELSE
    coll_total = sumdir - sumex/2.0
  END IF

! release workspace

  usew1 = .false.
  usew2 = .false.

  RETURN
END FUNCTION coll_total
