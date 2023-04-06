
! presently ONLY static version of 'sstep_lsic'
#ifdef REALSWITCH
#if(twostsic)
MODULE localize_rad

  USE params
  USE kinetic
  USE orthmat
  USE util, ONLY: wfovlp, project

  INTEGER, SAVE, PRIVATE :: kdim

! matrices of radial moments

  REAL(DP), PRIVATE, ALLOCATABLE :: rrmatr(:, :, :) ! matrix of r**2
  REAL(DP), PRIVATE, ALLOCATABLE :: xxmatr(:, :, :) ! matrix of x**2
  REAL(DP), PRIVATE, ALLOCATABLE :: yymatr(:, :, :) ! matrix of y**2
  REAL(DP), PRIVATE, ALLOCATABLE :: zzmatr(:, :, :) ! matrix of z**2
  REAL(DP), PRIVATE, ALLOCATABLE :: xmatr(:, :, :) ! matrix of x
  REAL(DP), PRIVATE, ALLOCATABLE :: ymatr(:, :, :) ! matrix of y
  REAL(DP), PRIVATE, ALLOCATABLE :: zmatr(:, :, :) ! matrix of z
  REAL(DP), PRIVATE, ALLOCATABLE :: vecsr(:, :, :) ! searched eigenvectors
  INTEGER, SAVE, PRIVATE :: ndim(2)

  INTERFACE superpose_state
! automated choice between REAL, COMPLEX and complexsic versions of superpose_state
    MODULE PROCEDURE superpose_state_r, superpose_state_c, superpose_state_rc!, superpose_state_rc_dummy
  END INTERFACE superpose_state

  INTERFACE avermatrix
    MODULE PROCEDURE avermatrix_r, avermatrix_c
  END INTERFACE avermatrix

CONTAINS
! ******************************

  SUBROUTINE init_radmatrix()
    IMPLICIT NONE
    kdim = kstate
    ALLOCATE (rrmatr(kdim, kdim, 2))
    ALLOCATE (xxmatr(kdim, kdim, 2))
    ALLOCATE (yymatr(kdim, kdim, 2))
    ALLOCATE (zzmatr(kdim, kdim, 2))
    ALLOCATE (xmatr(kdim, kdim, 2))
    ALLOCATE (ymatr(kdim, kdim, 2))
    ALLOCATE (zmatr(kdim, kdim, 2))
    ALLOCATE (vecsr(kdim, kdim, 2))

  END SUBROUTINE init_radmatrix

  SUBROUTINE sstep_lsic(q0, akv, aloc, iter, qnew)

! ******************************

! One static step with localized SIC.
! Requires previous CALL of 'calc_locsic' which
! produces 'qnew' as set of localized states on
! which the localized SIC potentials has acted
! and which are superposed TO the actual state
! using the unitary transformation IN 'vecs'
! IN COMMON /radmatrixr/.

! - ONLY serial version -

    IMPLICIT NONE

    REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(IN) :: akv(kdfull2)
    REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
    INTEGER, INTENT(IN) :: iter
    REAL(DP), INTENT(IN OUT) :: qnew(kdfull2, kstate)

    INTEGER :: i, is, ishift, istdiag, na, nae, nb, nbe, noff
    REAL(DP) :: adamP, cf, FFTnorm, sum0, sume, sumk, sum2, vol
    REAL(DP) :: time_init, time_fin, time_cpu
! workspaces

#if(gridfft)
    COMPLEX(DP), ALLOCATABLE :: q2(:), psipr(:)
#else
    REAL(DP), ALLOCATABLE :: q2(:)
#endif
    REAL(DP) :: coeff(kstate), eigval(kstate, 2)
    REAL(DP) :: h_matrix(kstate, kstate, 2)
    REAL(DP) :: eigvec(kstate, kstate, 2)
    REAL(DP), ALLOCATABLE :: q1(:)
!EQUIVALENCE (q1(1),w1(1))
!EQUIVALENCE (q2(1),w2(1)) ! occupies also w3

    LOGICAL, PARAMETER :: tcpu = .true., ttest = .false.
!DATA tcpu,ttest/.true.,.false./
    DATA fftnorm/0D0/

    LOGICAL, PARAMETER :: tdiag = .true.

    INTEGER, PARAMETER :: iterdiag = 100

!------------------------------------------------------------------------

! WRITE(*,*) ' SSTEP with IFSICP=',ifsicp
    IF (tcpu) CALL cpu_time(time_init)
    IF (ifsicp /= 7) STOP 'SSTEP_LSIC ONLY for localized SIC'
    istdiag = 1 ! 4*istinf

#if(gridfft)
    ALLOCATE (q2(kdfull2), psipr(kdfull2))
#else
    ALLOCATE (q2(kdfull2))
#endif
    ALLOCATE (q1(kdfull2))

    dvol = dx*dy*dz
    DO nbe = 1, nstate
      is = ispin(nbe)
      ishift = (ispin(nbe) - 1)*nxyz ! store spin=2 IN upper BLOCK

! ACTION of mean-field Hamiltonian on 'hwfr'

      IF (ipsptyp == 1) THEN
        CALL nonlocalr(q0(1, nbe), q1)
        enonlo(nbe) = wfovlp(q0(:, nbe), q1)
        DO i = 1, nxyz
          q1(i) = q1(i) + q0(i, nbe)*(aloc(i + ishift) - amoy(nbe))
        END DO
      ELSE

        DO i = 1, nxyz
          q1(i) = q0(i, nbe)*(aloc(i + ishift) - amoy(nbe))
        END DO
      END IF

! subtract SIC potential for state NBE

      DO na = 1, ndim(is)
        nb = nbe - (is - 1)*ndim(1)
        nae = na + (is - 1)*ndim(1)
        cf = vecsr(nb, na, is)
        DO i = 1, nxyz
          q1(i) = q1(i) - qnew(i, nae)*cf
        END DO
      END DO

! optionally compute expectation value of potential energy

      IF (MOD(iter, istinf) == 0) epotsp(nbe) = wfovlp(q0(:, nbe), q1) + amoy(nbe)

#if(gridfft)

! ACTION of the kinetic energy IN momentum space

      CALL rftf(q0(1, nbe), psipr)
      CALL rftf(q1, q2)
      DO i = 1, nxyz
        q2(i) = psipr(i)*akv(i) + q2(i)
      END DO

! optionally compute expectation value of kinetic energy
! and variance of mean field Hamiltonian

      IF (MOD(iter, istinf) >= 0) THEN
        sum0 = 0D0
        sumk = 0D0
        sume = 0D0
        sum2 = 0D0
        DO i = 1, nxyz
          vol = REAL(psipr(i), DP)*REAL(psipr(i), DP) + AIMAG(psipr(i))*AIMAG(psipr(i))
          sum0 = vol + sum0
          sumk = vol*akv(i) + sumk
          sume = REAL(q2(i), DP)*REAL(psipr(i), DP) + AIMAG(q2(i))*AIMAG(psipr(i)) + sume
          sum2 = REAL(q2(i), DP)*REAL(q2(i), DP) + AIMAG(q2(i))*AIMAG(q2(i)) + sum2
        END DO
        ekinsp(nbe) = sumk/sum0
        sume = sume/sum0
        sum2 = sum2/sum0
        fftnorm = sum0
! WRITE(6,'(a,i3,2(1pg12.4))')
! & ' nbe,norm,spe=',nbe,sum0,sume
! amoy(nbe) = sume
        evarsp(nbe) = SQRT(MAX(sum2 - sume**2, small))
! WRITE(6,'(a,i3,2(1pg12.4))')
! & ' nbe,spe,var=',nbe,sume,evarsp(nbe)
        CALL rfftback(q2, q1)
        CALL project(q1, q1, ispin(nbe), q0)
        evarsp2(nbe) = SQRT(wfovlp(q1, q1))
      ELSE IF (tdiag .AND. fftnorm > 0D0) THEN
        sume = 0D0
        DO i = 1, nxyz
          sume = REAL(q2(i), DP)*REAL(psipr(i), DP) + AIMAG(q2(i))*AIMAG(psipr(i)) + sume
        END DO
        sume = sume/fftnorm
        WRITE (6, '(a,i3,2(1pg12.4))') ' nbe,spe,amoy=', nbe, sume, amoy(nbe)
      END IF
      IF (tdiag .AND. fftnorm > 0D0 .AND. iter > iterdiag &
          .AND. MOD(iter, istdiag) == 0) THEN
        DO i = 1, nxyz
          q2(i) = q2(i) - sume*q0(i, nbe)
        END DO
      END IF

! collect Hamiltonian matrix

      IF (tdiag .AND. iter > iterdiag .AND. ndim(1) > 0 &
          .AND. MOD(iter, istdiag) == 0) THEN
        CALL rfftback(q2, q1)
        nb = nbe - (is - 1)*ndim(1)
        DO na = 1, ndim(is)
          nae = na + (is - 1)*ndim(1)
          h_matrix(na, nb, is) = wfovlp(q0(:, nae), q1)
          IF (na == nb) h_matrix(na, nb, is) = &
            h_matrix(na, nb, is) + amoy(nbe)
        END DO
        IF (ttest) WRITE (6, '(a,2i3,10(1pg12.4))') &
          'diag: ', is, nbe, (h_matrix(na, nb, is), na=1, ndim(is))
      END IF

! perform the damped gradient step

      IF (idyniter /= 0 .AND. iter > 100) e0dmp = MAX(ABS(amoy(nbe)), 0.5D0)
      IF (e0dmp > small) THEN
        DO i = 1, nxyz
          adamp = epswf/(akv(i) + e0dmp)
          psipr(i) = psipr(i) - adamp*q2(i)
        END DO
      ELSE
        DO i = 1, nxyz
          psipr(i) = psipr(i) - epswf*q2(i)
        END DO
      END IF

      CALL rfftback(psipr, q0(1, nbe))

#else
      STOP ' full SIC requires fast Fourier code '
#endif

    END DO

    DEALLOCATE (q1)
#if(gridfft)
    DEALLOCATE (q2, psipr)
#else
    DEALLOCATE (q2)
#endif

! Schmidt ortho-normalisation

    CALL schmidt(q0)

! optionally diagonalization amongst occupied states

    IF (tdiag .AND. iter > iterdiag .AND. ndim(1) > 0 &
        .AND. MOD(iter, istdiag) == 0) THEN
      DO is = 1, 2
        DO na = 1, ndim(is)
          DO nb = 1, na - 1
            h_matrix(na, nb, is) = 0.5D0*(h_matrix(na, nb, is) + h_matrix(nb, na, is))
            h_matrix(nb, na, is) = h_matrix(na, nb, is)
          END DO
        END DO
        CALL vecgradstep(h_matrix(1, 1, is), eigvec(:, :, is), ndim(is), kstate, eigval(1, is))
        IF (ttest) THEN
          WRITE (6, '(a,i3)') 'after diag: is=', is
          DO nb = 1, ndim(is)
            WRITE (6, '(10(1pg12.4))') (eigvec(na, nb, is), na=1, ndim(is))
          END DO
          WRITE (6, '(10(1pg12.4))') (eigval(na, is), na=1, ndim(is))
        END IF
      END DO

      DO nbe = 1, nstate
        is = ispin(nbe)
        nb = nbe - (is - 1)*ndim(1)
        noff = (is - 1)*ndim(1) + 1
        CALL superpose_state(qnew(:, nbe), eigvec(:, nb, is), q0, is)
        amoy(nbe) = eigval(nb, is)
        epotsp(nbe) = amoy(nbe) - ekinsp(nbe)
      END DO
      DO nbe = 1, nstate
        DO i = 1, nxyz
          q0(i, nbe) = qnew(i, nbe)
        END DO
      END DO
    END IF

    IF (tcpu) THEN
      CALL cpu_time(time_fin)
      time_cpu = time_fin - time_init
! WRITE(6,'(a,1pg13.5)') ' CPU time IN SSTEP',time_cpu
! WRITE(7,'(a,1pg13.5)') ' CPU time IN SSTEP',time_cpu
    END IF
    WRITE (6, '(a,i5,6(f10.4))') 'iter,up/down,CPU=', iter, se(4), se(5), time_cpu
    WRITE (7, '(a,i5,6(f10.4))') 'iter,up/down,CPU=', iter, se(4), se(5), time_cpu

    RETURN
  END SUBROUTINE sstep_lsic

!-----calc_locsic--------------------------------------------------

  SUBROUTINE calc_locsic(q0, qsic)

! localized SIC:
! input is set of wavefunctions on 'q0'
! output are SIC s.p. wavefunctions on 'qsic'

    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(OUT) :: qsic(kdfull2, kstate)

    INTEGER :: idx, ind, is, ishift, nb, nbeff
    REAL(DP) :: enrear1, enrear2, enrearsave, enerpwsave, enpw1, enpw2, encadd

! workspaces

    REAL(DP), ALLOCATABLE :: usicsp(:), rhosp(:)
!COMMON /sicwork/ rhospu(2*kdfull2),rhospd(2*kdfull2), &
! chpdftspu(2*kdfull2),chpdftspd(2*kdfull2)

    LOGICAL, PARAMETER :: ttest = .false.

!------------------------------------------------------------------

! Compute matrix of radial moments and determine
! optimally localizing transformation.
! Results is unitary matrix 'vecs' communicated via
! COMMON /radmatrixr/.

    CALL spmomsmatrixo(q0)
    CALL locgradstep(1, 0)
    CALL locgradstep(2, 0)

! loop over localized basis:
! computation of localized-state SIC potentials
! and
! accumulation of energy contributions

    ALLOCATE (usicsp(2*kdfull2), rhosp(2*kdfull2))

    enrearsave = enrear
    enerpwsave = enerpw
    enrear1 = 0D0
    enrear2 = 0D0
    enpw1 = 0D0
    enpw2 = 0D0
    encadd = 0D0

    DO nb = 1, nstate
      IF (occup(nb) /= 1.0) THEN
        STOP 'localized SIC requires all OCCUP=1'
      ELSE
        is = ispin(nrel2abs(nb))
        ishift = (is - 1)*nxyz ! store spin=2 IN upper BLOCK
        nbeff = nb - (is - 1)*ndim(1)
        CALL superpose_state(qsic(:, nb), vecsr(:, nbeff, is), q0, is)
#ifdef REALSWITCH
        CALL calc_sicspr(rhosp, usicsp, qsic(1, nb), nb)
#else
        STOP 'localized SIC not yet for dynamics CASE'
#endif

        IF (ispin(nrel2abs(nb)) == 1) THEN
          enrear1 = enrear1 + enrear*occup(nb)
          enpw1 = enpw1 + enerpw*occup(nb)
        ELSE
          enrear2 = enrear2 + enrear*occup(nb)
          IF (directenergy) THEN
            enpw2 = enpw2 + enerpw*occup(nb)
          END IF
        END IF
        IF (directenergy) THEN
          encadd = encadd + encoulsp*occup(nb)
        END IF

        IF (ispin(nrel2abs(nb)) == 1) THEN
          DO ind = 1, nxyz
            qsic(ind, nb) = usicsp(ind)*qsic(ind, nb)
          END DO
        ELSE
          DO ind = 1, nxyz
            idx = ind + nxyz
            qsic(ind, nb) = usicsp(idx)*qsic(ind, nb)
          END DO
        END IF
      END IF
    END DO
    encadd = encadd/2D0
    enrear = enrearsave - enrear1 - enrear2
    IF (directenergy) THEN
      enerpw = enerpwsave - enpw1 - enpw2 - encadd
    END IF

    DEALLOCATE (usicsp, rhosp)

    RETURN
  END SUBROUTINE calc_locsic

!-----analyze_mom--------------------------------------------------

  SUBROUTINE analyze_mom(rhoin, uin, urmoms, tprint)

! Analyzes center of gravity for density 'rhoin' times
! potential 'uin'. The switch 'tprint' regulates the
! level of output.
! The accumulated moments are returned on 'urmoms'.

!USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: rhoin(kdfull2)
    REAL(DP), INTENT(IN) :: uin(kdfull2)
    REAL(DP), INTENT(OUT) :: urmoms(0:3)
    LOGICAL, INTENT(IN) :: tprint

    INTEGER :: i, ind, ix, iy, iz
    REAL(DP) :: x1, y1, z1
    REAL(DP), ALLOCATABLE :: urho(:)

!---------------------------------------------------------------------

    ALLOCATE (urho(kdfull2))
    DO i = 0, 3
      urmoms(i) = 0D0
    END DO
    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy
        DO ix = minx, maxx
          ind = ind + 1
          urho(ind) = rhoin(ind)*uin(ind)
          IF ((ix /= nx2) .AND. (iy /= ny2) .AND. (iz /= nz2)) THEN
            x1 = (ix - nxsh)*dx
            urmoms(0) = urmoms(0) + urho(ind)
            urmoms(1) = urmoms(1) + urho(ind)*x1
            urmoms(2) = urmoms(2) + urho(ind)*y1
            urmoms(3) = urmoms(3) + urho(ind)*z1
            IF (tprint .AND. iz == nzsh .AND. iy == nysh) THEN
              WRITE (6, '(i3,2(1pg13.5))') ix, rhoin(ind), urho(ind)
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i = 0, 3
      urmoms(i) = urmoms(i)*dvol
    END DO
    DEALLOCATE (urho)
    RETURN

  END SUBROUTINE analyze_mom

!-----superpose_state--------------------------------------------------
!
! Superposition TO new state:
! wfsup new single-particle state
! coeff superposition coefficients
! q0 set of s.p.states TO be combined
! is spin of states
!----------------------------------------------------------------------
! REAL version
!----------------------------------------------------------------------
  SUBROUTINE superpose_state_r(wfsup, coeff, q0, is)
    USE params
    USE kinetic
    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: wfsup(kdfull2)
    REAL(DP), INTENT(IN) :: coeff(kstate)
    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
    INTEGER, INTENT(IN) ::is

    INTEGER :: i, na, nas
!---------------------------------------------------------------------
    wfsup(1:nxyz) = 0.0_DP

    DO nas = 1, ndims(is)
      na = nas + (is - 1)*ndims(1)
      DO i = 1, nxyz
        wfsup(i) = wfsup(i) + coeff(nas)*q0(i, na)
      END DO
    END DO

    RETURN
  END SUBROUTINE superpose_state_r

!----------------------------------------------------------------------
! COMPLEX version
!----------------------------------------------------------------------
  SUBROUTINE superpose_state_c(wfsup, coeff, q0, is)

    USE params
    USE kinetic
    IMPLICIT NONE

    COMPLEX(DP), INTENT(OUT) :: wfsup(kdfull2)
    COMPLEX(DP), INTENT(IN) :: coeff(kstate)
    COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
    INTEGER, INTENT(IN) ::is

    INTEGER :: i, na, nas
!---------------------------------------------------------------------

    wfsup(1:nxyz) = (0.0_DP, 0.0_DP)

    DO nas = 1, ndims(is)
      na = nas + (is - 1)*ndims(1)
      DO i = 1, nxyz
        wfsup(i) = wfsup(i) + coeff(nas)*q0(i, na)
      END DO
    END DO

    RETURN
  END SUBROUTINE superpose_state_c

!----------------------------------------------------------------------
! cmplxsic version
!----------------------------------------------------------------------
  SUBROUTINE superpose_state_rc(wfsup, coeff, q0, is)

    USE params
    USE kinetic
    IMPLICIT NONE

    COMPLEX(DP), INTENT(OUT) :: wfsup(kdfull2)
    COMPLEX(DP), INTENT(IN) :: coeff(kstate)
    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
    INTEGER, INTENT(IN) ::is

    INTEGER :: i, na, nas
!---------------------------------------------------------------------

    wfsup(1:nxyz) = (0.0_DP, 0.0_DP)

    DO nas = 1, ndims(is)
      na = nas + (is - 1)*ndims(1)
      DO i = 1, nxyz
        wfsup(i) = wfsup(i) + coeff(nas)*q0(i, na)
      END DO
    END DO

    RETURN
  END SUBROUTINE superpose_state_rc

!-----spmomsmatrix----------------------------------------------spmoms

  SUBROUTINE spmomsmatrixo(wfr)

! Matrix of spatial moments between single-particle states
! from REAL wf's:
! wfr = set of REAL single particle wavefunctions
! The RESULT is stored IN COMMON/radmatrixr/ for further
! USE IN localization transformation.

!USE params
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: wfr(kdfull2, kstate)

    INTEGER, PARAMETER :: iunit = 0

    LOGICAL, SAVE :: tfirst = .true.

    INTEGER :: i, ind, is, ix, iy, iz, ma, mb, na, nb, noff
    REAL(DP) :: s, wfmom, xmom, ymom, zmom, xxmom, yymom, zzmom
    REAL(DP) :: x1, x2, y1, y2, z1, z2
!----------------------------------------------------------------------

    OPEN (iunit, POSITION='append', FILE='pstatmomsmatrix.'//outnam)

! check spin of states

    ndim(1) = 0
    ndim(2) = 0
    DO na = 1, nstate
      IF (ispin(nrel2abs(na)) == 1) THEN
        IF (ndim(2) /= 0) STOP ' spins OUT of order'
        ndim(1) = na
      ELSE
        IF (ndim(1) == 0) STOP ' spins OUT of order'
        ndim(2) = na - ndim(1)
      END IF
    END DO
    noff = ndim(1)
    IF (iunit > 0) WRITE (iunit, '(a,2i5)') 'NDIM=', ndim

! compute matrix elements and store

    IF (iunit > 0) THEN
      WRITE (iunit, '(a)') 'matrix of s.p. moments:', &
        ' na nb x y z xx yy zz'
      tfirst = .false.
    END IF
    DO na = 1, nstate
      DO nb = 1, na
        IF (ispin(nrel2abs(nb)) == ispin(nrel2abs(na))) THEN
          wfmom = 0D0
          xmom = 0D0
          ymom = 0D0
          zmom = 0D0
          xxmom = 0D0
          yymom = 0D0
          zzmom = 0D0

          ind = 0
          DO iz = minz, maxz
            z1 = (iz - nzsh)*dz
            z2 = z1*z1
            DO iy = miny, maxy
              y1 = (iy - nysh)*dy
              y2 = y1*y1
              DO ix = minx, maxx
                ind = ind + 1
                IF ((ix /= nx2) .AND. (iy /= ny2) .AND. (iz /= nz2)) THEN
                  x1 = (ix - nxsh)*dx
                  x2 = x1*x1
                  s = wfr(ind, nb)*wfr(ind, na)
! monopole
                  wfmom = wfmom + s
! dipole
                  xmom = xmom + s*x1
                  ymom = ymom + s*y1
                  zmom = zmom + s*z1
! quadrupole
                  xxmom = xxmom + s*x2
                  yymom = yymom + s*y2
                  zzmom = zzmom + s*z2
                END IF
              END DO
            END DO
          END DO

          wfmom = dvol*wfmom
          xmom = dvol*xmom
          ymom = dvol*ymom
          zmom = dvol*zmom
          xxmom = dvol*xxmom
          yymom = dvol*yymom
          zzmom = dvol*zzmom

          IF (iunit > 0) WRITE (iunit, '(3i4,6(1pg13.5))') &
            na, nb, ispin(nrel2abs(nb)), xmom, ymom, zmom, xxmom, yymom, zzmom

          is = ispin(nrel2abs(nb))
          IF (is == 1) THEN
            ma = na
            mb = nb
          ELSE
            ma = na - noff
            mb = nb - noff
          END IF
          IF (iunit > 0) WRITE (iunit, '(a,3i4)') 'ma,mb=', ma, mb
          xmatr(ma, mb, is) = xmom
          xmatr(mb, ma, is) = xmom
          ymatr(ma, mb, is) = ymom
          ymatr(mb, ma, is) = ymom
          zmatr(ma, mb, is) = zmom
          zmatr(mb, ma, is) = zmom
          xxmatr(ma, mb, is) = xxmom
          xxmatr(mb, ma, is) = xxmom
          yymatr(ma, mb, is) = yymom
          yymatr(mb, ma, is) = yymom
          zzmatr(ma, mb, is) = zzmom
          zzmatr(mb, ma, is) = zzmom
        END IF
      END DO
    END DO
    IF (iunit > 0) CLOSE (iunit)

    RETURN
  END SUBROUTINE spmomsmatrixo

!-----locgradstep-----------------------------------------

  SUBROUTINE locgradstep(is, iprint)
    USE params
    USE orthmat, ONLY: ovlpmatrix
    IMPLICIT NONE
! Nonlinear gradient iteration TO optimally localized states:
! 'vecs' system of eigen-vectors TO be determined
! 'is' isospin
! 'iprint' PRINT level: <0 --> no PRINT at all
! 0 --> ONLY final RESULT
! >0 --> PRINT modulus
! The matrices and DIMENSION are communicated via
! COMMON /radmatrixr/.

    INTEGER, INTENT(IN) :: is
    INTEGER, INTENT(IN) :: iprint

    INTEGER, PARAMETER :: itermax = 2000
    REAL(DP), PARAMETER :: step = 0.5D0
    REAL(DP), PARAMETER :: precis = 1.0D-14

    LOGICAL, PARAMETER :: ttest = .false.

    REAL(DP) :: xaver(kdim), yaver(kdim), zaver(kdim)
    REAL(DP) :: raver(kdim), rvary(kdim)
    REAL(DP) :: xlambda(kdim, kdim)
    REAL(DP) :: variance, variance2 ! variance of step
    REAL(DP) :: radvary ! variance of radii
    REAL(DP) :: radmax ! max. squared radius

    REAL(DP) :: w(kdim) ! workspace
    REAL(DP) :: varstate(kdim), averstate(kdim)
    REAL(DP) :: acc, actstep
    REAL(DP) :: xbar, ybar, zbar
    INTEGER :: i, j ! eigenstates
    INTEGER :: na, nb ! matrix indices
    INTEGER :: iter

!-------------------------------------------------------

    IF (ttest) WRITE (6, '(10(/3i4,6(1pg13.5)))') ((na, nb, is, &
                                                    xmatr(na, nb, is), ymatr(na, nb, is), zmatr(na, nb, is), &
                                                    xxmatr(na, nb, is), yymatr(na, nb, is), zzmatr(na, nb, is), &
                                                    na=1, ndim(is)), nb=1, ndim(is))

! initialize eigen-vectors

    DO na = 1, ndim(is)
      DO nb = 1, ndim(is)
        rrmatr(na, nb, is) = xxmatr(na, nb, is) + yymatr(na, nb, is) &
                             + zzmatr(na, nb, is)
        IF (na == nb) THEN
          vecsr(na, nb, is) = 1.0D0
        ELSE
          vecsr(na, nb, is) = 0.0D0
        END IF
      END DO
    END DO

    IF (iprint > 0) WRITE (6, '(a)') ' iter variance '
    DO iter = 1, itermax

! averages of x,y,z

      radvary = 0.0D0
      radmax = 0.0D0
      DO i = 1, ndim(is)
        xaver(i) = avermatrix(xmatr(:, :, is), vecsr(:, i, is), ndim(is), kdim)
        yaver(i) = avermatrix(ymatr(:, :, is), vecsr(:, i, is), ndim(is), kdim)
        zaver(i) = avermatrix(zmatr(:, :, is), vecsr(:, i, is), ndim(is), kdim)
        raver(i) = avermatrix(rrmatr(:, :, is), vecsr(:, i, is), ndim(is), kdim)
        rvary(i) = raver(i)**2 - xaver(i)**2 - yaver(i)**2 - zaver(i)**2
        radvary = radvary + rvary(i)
        radmax = MAX(radmax, raver(i))
      END DO
      IF (iprint > 0 .AND. iter == 1) THEN
        WRITE (6, '(a,i3)') ' initial for spin IS=', is, &
          ' state x y z variance '
        WRITE (6, '(4(1pg12.4))') (xaver(i), yaver(i), zaver(i), SQRT(rvary(i)), &
                                   i=1, ndim(is))
        WRITE (6, '(36x,1pg12.4)') SQRT(radvary/ndim(is))
      END IF
      IF (ttest) THEN
        WRITE (6, '(a)') 'aver x y z rr rvary:'
        WRITE (6, '(5(1pg12.4))') (xaver(i), yaver(i), zaver(i), raver(i), rvary(i), &
                                   i=1, ndim(is))
      END IF

! matrix of Lagrangian multipliers

      DO i = 1, ndim(is)
        DO j = 1, ndim(is)
          xlambda(i, j) = ovlpmatrix(rrmatr(:, :, is), vecsr(:, i, is), &
                                     vecsr(:, j, is), ndim(is), kdim) - 0.5D0*(xaver(i) + xaver(j))* &
                          ovlpmatrix(xmatr(:, :, is), vecsr(:, i, is), vecsr(:, j, is), ndim(is), kdim) &
                          - 0.5D0*(yaver(i) + yaver(j))*ovlpmatrix(ymatr(:, :, is), vecsr(:, i, is), &
                                                                   vecsr(:, j, is), ndim(is), kdim) - 0.5D0*(zaver(i) + zaver(j))* &
                          ovlpmatrix(zmatr(:, :, is), vecsr(:, i, is), vecsr(:, j, is), ndim(is), kdim)
        END DO
      END DO
      IF (ttest) THEN
        WRITE (6, '(a)') 'xlambda:'
        WRITE (6, '(20(1pg12.4))') ((xlambda(i, j), i=1, ndim(is)), j=1, ndim(is))
      END IF

! application of radii-matrix

      variance = 0.0D0
      variance2 = 0.0D0
      DO i = 1, ndim(is)
        DO na = 1, ndim(is)
! apply matrix
          acc = 0.0D0
          DO nb = 1, ndim(is)
            acc = (rrmatr(na, nb, is) - xaver(i)*xmatr(na, nb, is) &
                   - yaver(i)*ymatr(na, nb, is) - zaver(i)*zmatr(na, nb, is)) &
                  *vecsr(nb, i, is) + acc
          END DO
! subtract constraints
          DO j = 1, ndim(is)
            acc = acc - xlambda(i, j)*vecsr(na, j, is)
          END DO

          w(na) = acc
        END DO
        averstate(i) = vecovlp(vecsr(:, i, is), w, ndim(1))
        varstate(i) = vecnorm(w, ndim(1))
        variance = variance + varstate(i) - averstate(i)**2
        variance2 = variance2 + varstate(i)
        IF (ttest) THEN
          WRITE (6, '(a,(2(1pg12.4)))') i, averstate(i), varstate(i)
        END IF
! step
        actstep = step/radmax
        DO na = 1, ndim(is)
          vecsr(na, i, is) = vecsr(na, i, is) - actstep*w(na)
        END DO
      END DO
      IF (iprint > 0 .AND. MOD(iter, iprint) == 0) WRITE (6, '(i5,4(1pg12.4))') &
        iter, SQRT(variance), SQRT(variance2), SQRT(radvary), radmax
      IF (variance < precis) EXIT

! ortho-normalization

      CALL orthnorm(vecsr(:, :, is), ndim(is), kdim)

    END DO

    IF (iprint >= 0) THEN
      WRITE (6, '(a,i4)') ' at iteration nr.', iter, &
        ' results for spin IS=', is, ' state x y z variance '
      xbar = 0D0
      ybar = 0D0
      zbar = 0D0
      DO i = 1, ndim(is)
        WRITE (6, '(i3,4(1pg15.7))') i, xaver(i), yaver(i), zaver(i), SQRT(rvary(i))
        xbar = xbar + xaver(i)
        ybar = ybar + yaver(i)
        zbar = zbar + zaver(i)
      END DO
      WRITE (6, '(3x,4(1pg15.7))') xbar, ybar, zbar, SQRT(radvary/ndim(is))
    END IF

    RETURN
  END SUBROUTINE locgradstep

!-----vecgradstep-----------------------------------------

  SUBROUTINE vecgradstep(a, vecs, ndim, kdimin, averstate)

    IMPLICIT NONE

    REAL(DP), INTENT(IN OUT) :: a(kdimin, kdimin)
    REAL(DP), INTENT(IN OUT) :: vecs(kdimin, kdimin)
    INTEGER, INTENT(IN) :: ndim
    INTEGER, INTENT(IN) :: kdimin
    REAL(DP), INTENT(OUT) :: averstate(kdimin)

! One gradient iteration step:
! 'a' matrix TO be diagionalized
! 'vecs' system of vectors on which the step acts
! 'ndim' DIMENSION of 'a' and 'vecs'

    INTEGER, PARAMETER :: itermax = 2000

    INTEGER, PARAMETER :: mprint = -1
    REAL(DP), PARAMETER :: step = 0.5D0
    REAL(DP), PARAMETER :: precis = 1.0D-14

    REAL(DP) :: vecsr(kdimin, kdimin)

    REAL(DP) :: variance

    REAL(DP) :: w(kdimin) ! workspace
    REAL(DP) :: varstate(kdimin)
    REAL(DP) :: aver, aver2, actstep, avmax, avmin
!REAL(DP) :: vecovlp,vecnorm ! FUNCTION names ????
    INTEGER :: i, j, n, iter

!-------------------------------------------------------

    IF (mprint > 0) WRITE (6, '(a)') ' iter variance '
    DO i = 1, ndim
! WRITE(6,'(10(1pg12.4))') (a(j,i),j=1,ndim)
      DO j = 1, ndim
        IF (i == j) THEN
          vecsr(i, j) = 1.0D0
        ELSE
          vecsr(i, j) = 0.0D0
        END IF
      END DO
    END DO
    DO iter = 1, itermax
      IF (iter == 1) THEN
        actstep = step/10.0D0
      ELSE
        actstep = step/(avmax - avmin)
      END IF
      variance = 0.0D0
      avmax = -1.0D30
      avmin = 1.0D30
      DO i = 1, ndim
        CALL operate(a, vecsr(1, i), w, ndim, kdimin)
        aver = vecovlp(vecsr(:, i), w, ndim)
        averstate(i) = aver
        avmax = MAX(aver, avmax)
        avmin = MIN(aver, avmin)
        aver2 = vecnorm(w, ndim)
        varstate(i) = (aver2 - aver**2)
        variance = variance + varstate(i)
        DO n = 1, ndim
          vecsr(n, i) = vecsr(n, i) - step*(w(n) - aver*vecsr(n, i))
        END DO
      END DO
! variance = sqrt(variance)
      IF (mprint > 0 .AND. MOD(iter, mprint) == 0) WRITE (6, '(i5,2(1pg12.4))') &
        iter, variance, SQRT(variance)
! WRITE(6,*) averstate

! ortho-normalization

      CALL orthnorm(vecs, ndim, kdimin)
      IF (variance < precis) EXIT

    END DO
    IF (mprint >= 0) WRITE (6, '(a,i5,1pg12.4)') &
      'DIAG: iter,variance=', iter, SQRT(variance)

    RETURN
  END SUBROUTINE vecgradstep

!-----operate-----------------------------------------

  SUBROUTINE operate(a, vecin, vecout, ndim, kdimin)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: a(kdimin, kdimin)
    REAL(DP), INTENT(IN) :: vecin(kdimin)
    REAL(DP), INTENT(OUT) :: vecout(kdimin)
    INTEGER, INTENT(IN) :: ndim
    INTEGER, INTENT(IN) :: kdimin

! apply OPERATOR 'a' on vector 'vecin' and RETURN
! RESULT on 'vecout'.
! DIMENSION 'ndim'.

    REAL(DP) :: acc
    INTEGER :: n, m

!-------------------------------------------------------

    DO n = 1, ndim
      vecout(n) = 0.0D0
      DO m = 1, ndim
        vecout(n) = vecout(n) + a(n, m)*vecin(m)
      END DO
    END DO

    RETURN
  END SUBROUTINE operate

!-----avermatrix-----------------------------------------
!
! Average of matrix 'a' taken with state 'vec' of
! actual length 'ndim', dimensioned with 'kdimin'.
!-------------------------------------------------------
! REAL version
!-------------------------------------------------------
  REAL(DP) FUNCTION avermatrix_r(a, vec, ndim, kdimin)
    USE params, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: a(kdimin, kdimin)
    REAL(DP), INTENT(IN) :: vec(kdimin)
    INTEGER, INTENT(IN) :: ndim
    INTEGER, INTENT(IN) :: kdimin
    REAL(DP) :: aver
    INTEGER :: i, j

    aver = 0.0D0
    DO i = 1, ndim
      DO j = 1, ndim
        aver = aver + vec(j)*a(j, i)*vec(i)
      END DO
    END DO
    avermatrix_r = aver

    RETURN
  END FUNCTION avermatrix_r
!-------------------------------------------------------
! COMPLEX version
!-------------------------------------------------------
  REAL(DP) FUNCTION avermatrix_c(a, vec, ndim, kdimin)
    USE params, ONLY: DP
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: a(kdimin, kdimin)
    COMPLEX(DP), INTENT(IN) :: vec(kdimin)
    INTEGER, INTENT(IN) :: ndim
    INTEGER, INTENT(IN) :: kdimin

    COMPLEX(DP) :: aver
    INTEGER :: i, j

    aver = CMPLX(0D0, 0D0, DP)
    DO i = 1, ndim
      DO j = 1, ndim
        aver = aver + CONJG(vec(j))*a(j, i)*vec(i)
      END DO
    END DO
    avermatrix_c = aver

    RETURN
  END FUNCTION avermatrix_c

!-----calc_locwf--------------------------------------------------

  SUBROUTINE calc_locwf(q0, qnew)

! computes localized wavefunctions
! input is set of wavefunctions on 'q0'
! output are localized wavefunctions on 'q0'
! the array 'qnew' is used as workspace

    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(IN OUT) :: qnew(kdfull2, kstate)

    LOGICAL, PARAMETER :: ttest = .false.
    INTEGER :: ind, is, ishift, nb, nbeff
!------------------------------------------------------------------

! Compute matrix of radial moments and determine
! optimally localizing transformation.
! Results is unitary matrix 'vecs' communicated via
! COMMON /radmatrixr/.

    IF (ifspemoms) WRITE (6, '(a)') ' before localization'

    CALL spmomsr(q0, 6)

    CALL spmomsmatrixo(q0)
    CALL locgradstep(1, 0)
    CALL locgradstep(2, 0)

    DO nb = 1, nstate
      is = ispin(nrel2abs(nb))
      ishift = (is - 1)*nxyz ! store spin=2 IN upper BLOCK
      nbeff = nb - (is - 1)*ndim(1)
      CALL superpose_state(qnew(:, nb), vecsr(:, nbeff, is), q0, is)
    END DO

    DO nb = 1, nstate
      DO ind = 1, nxyz
        q0(ind, nb) = qnew(ind, nb)
      END DO
    END DO

    IF (ifspemoms) THEN
      WRITE (6, '(a)') ' after localization'
      CALL spmomsr(q0, 6)
    END IF

    RETURN
  END SUBROUTINE calc_locwf

END MODULE localize_rad
#endif
#endif

SUBROUTINE dummy9
  RETURN
END SUBROUTINE dummy9

