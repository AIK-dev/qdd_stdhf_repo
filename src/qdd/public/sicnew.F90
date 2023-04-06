
! ******************************

#ifdef REALSWITCH
SUBROUTINE calc_sicr(rho, aloc, q0)
#else
  SUBROUTINE calc_sic(rho, aloc, q0)
#endif

! ******************************

! Switchboard for the various brands of SIC.
! Cases are selected according to 'ifsicp' which is
! communicated via module 'params'.
!
! List parameters (handed through to calling routines):
! rho = local density
! aloc = local KS potential
! q0 = set of s.p. wavefunctions

    USE params
    USE kinetic
    IMPLICIT NONE
#ifdef REALSWITCH
    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
#else
    COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
#endif
    REAL(DP), INTENT(IN OUT) :: aloc(*), rho(*)

!------------------------------------------------------------

    IF (tfreezekspot .AND. tfs > 0) RETURN

#ifdef REALSWITCH
    SELECT CASE (ifsicp)
    CASE (2)
      CALL calc_adsicr(rho, aloc)
    CASE (3)
      CALL calc_slaterr(rho, aloc, q0)
    CASE (4)
      CALL calc_sicklir(rho, aloc, q0)
    CASE DEFAULT
      RETURN
    END SELECT
#else
    SELECT CASE (ifsicp)
    CASE (2)
      CALL calc_adsic(rho, aloc)
    CASE (3)
      CALL calc_slater(rho, aloc, q0)
    CASE (4)
      CALL calc_sickli(rho, aloc, q0)
    CASE DEFAULT
      RETURN
    END SELECT
#endif
    RETURN
#ifdef REALSWITCH
  END SUBROUTINE calc_sicr
#else
END SUBROUTINE calc_sic
#endif

! ******************************

#ifdef REALSWITCH
SUBROUTINE calc_adsicr(rho, aloc)
#else
  SUBROUTINE calc_adsic(rho, aloc)
#endif

! Computes and adds ADSIC to local KS potential.
!
! Input:
! rho = local density
! Input/Output:
! aloc = local KS potential (TO be modified by ADSIC)

    USE params
    USE kinetic
    USE coulsolv, ONLY: solv_poisson_f, solv_poisson_e, tcoulfalr

    IMPLICIT NONE

    REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
    REAL(DP), INTENT(IN) :: rho(2*kdfull2)


    REAL(DP), DIMENSION(:), ALLOCATABLE :: rhosp, chpdftsp, coulsum, couldif
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rho1, rho2
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rhospu, rhospd, chpdftspu, chpdftspd
    INTEGER :: ind, idx, npartdw, npartup, npartto, SIZE
    INTEGER :: indpri
    REAL(DP) :: enrearsave, enerpwsave, enrear1, enrear2, enpw1, enpw2, encadd
    REAL(DP) :: couldw, coulup, fac, facdw, facup, facdwh, facuph, factotal

    LOGICAL,PARAMETER :: ttest=.FALSE.

!-------------------------------------------------------------------------

    IF (.NOT. ALLOCATED(akv)) STOP ' ADSIC requires FFT'

    IF (ifsicp /= 2) THEN
      STOP ' CALC_SIC called with wrong option IFSICP'
    END IF

! set electron numbers to code-internal valiables
    npartto = nelect
    npartdw = nspdw
    npartup = nelect - nspdw

    SIZE = nxyz

    enrearsave = enrear
    enerpwsave = enerpw
    enrear1 = 0D0
    enrear2 = 0D0
    enpw1 = 0D0
    enpw2 = 0D0
    encadd = 0D0
    IF (numspin == 2) THEN

      ALLOCATE (rhospu(2*kdfull2), rhospd(2*kdfull2), &
                chpdftspu(2*kdfull2), chpdftspd(2*kdfull2))

! averaged spinup and spindown density

      IF (npartup > 0) THEN
        facuph = 0.5D0/npartup

        DO ind = 1, SIZE
          rhospu(ind) = rho(ind)*(1D0 + rho(ind + nxyz))*facuph
          rhospu(ind + nxyz) = 1D0
        END DO
      END IF

      IF (npartdw > 0) THEN
        facdwh = 0.5D0/npartdw

        DO ind = 1, SIZE
          rhospd(ind) = rho(ind)*(1D0 - rho(ind + nxyz))*facdwh
          rhospd(ind + nxyz) = -1D0
        END DO
      END IF

! DFT for averaged spinup and spindown density

      IF (npartup > 0) THEN
        CALL calc_lda(rhospu, chpdftspu)
        enrear1 = enrear
        enpw1 = enerpw
      END IF
      IF (npartdw > 0) THEN
        CALL calc_lda(rhospd, chpdftspd)
        enrear2 = enrear
        enpw2 = enerpw
      END IF

      enrear = enrearsave - enrear1*npartup - enrear2*npartdw
      enerpw = enerpwsave - enpw1*npartup - enpw2*npartdw
      IF(ttest) WRITE (6, *) ' enrear.s=', enrearsave, enrear1, enrear2, enrear
      IF(directenergy .AND. ttest) &
       WRITE (6, *) ' enerpw.s=', enerpwsave, enpw1, enpw2, enerpw


      IF (npartup > 0) THEN
        DO ind = 1, SIZE
          aloc(ind) = aloc(ind) - chpdftspu(ind)
        END DO
      END IF

      IF (npartdw > 0) THEN
        DO idx = nxyz + 1, 2*nxyz
          aloc(idx) = aloc(idx) - chpdftspd(idx)
        END DO
      END IF

      DEALLOCATE (rhospu, rhospd, chpdftspu, chpdftspd)

    ELSE

      ALLOCATE (chpdftsp(2*kdfull2))
      ALLOCATE (rhosp(2*kdfull2))

      factotal = 1D0/npartto


      DO ind = 1, SIZE
        rhosp(ind) = rho(ind)*factotal
        rhosp(ind + nxyz) = 1D0
      END DO


! DFT for averaged s.p. state

      CALL calc_lda(rhosp, chpdftsp)

! subtract from pure LDA potential

      DO ind = 1, SIZE
        aloc(ind) = aloc(ind) - chpdftsp(ind)
      END DO

      enrear = enrearsave - enrear*npartto
      DEALLOCATE (chpdftsp)
      DEALLOCATE (rhosp)

    END IF

! correct Coulomb potential by 1/N

    IF (numspin == 2) THEN
      ALLOCATE (coulsum(kdfull2))
      ALLOCATE (couldif(kdfull2))
      ALLOCATE (rho1(kdfull2))
      ALLOCATE (rho2(kdfull2))

! compute Coulomb

      DO ind = 1, SIZE
        rho2(ind) = rho(ind)
        rho1(ind) = rho(ind)*rho(ind + nxyz)
      END DO

      IF (tcoulfalr) THEN
        CALL solv_poisson_f(rho1, couldif, kdfull2)
      ELSE
        CALL solv_poisson_e(rho1, couldif, kdfull2)
      END IF
        IF (nion2 > 0) THEN
          coulsum = chpcoul
        ELSE
          IF (tcoulfalr) THEN
            CALL solv_poisson_f(rho2, coulsum, kdfull2)
          ELSE
            CALL solv_poisson_e(rho2, coulsum, kdfull2)
          END IF
        END IF
        facup = 1D0/npartup

        DO ind = 1, SIZE
          coulup = 0.5D0*(coulsum(ind) + couldif(ind))
          aloc(ind) = aloc(ind) - coulup*facup
        END DO

        IF (npartdw > 0) THEN

          facdw = 1D0/npartdw
          DO ind = 1, SIZE
            couldw = 0.5D0*(coulsum(ind) - couldif(ind))
            idx = ind + nxyz
            aloc(idx) = aloc(idx) - facdw*couldw
          END DO
        IF(directenergy) THEN
          encadd = (facup*SUM((rho2(1:kdfull2)+rho1(1:kdfull2))*&
                              (coulsum(1:kdfull2)+couldif(1:kdfull2))) + &
                    facdw*SUM((rho2(1:kdfull2)-rho1(1:kdfull2))*&
                              (coulsum(1:kdfull2)-couldif(1:kdfull2))) )* &
                   dvol/2D0/4D0
          IF(ttest) WRITE(*,*) ' enerpw,encadd=',enerpw,encadd,enerpw - encadd
          enerpw = enerpw - encadd
        END IF

        END IF

        DEALLOCATE (coulsum)
        DEALLOCATE (couldif)
        DEALLOCATE (rho1)
        DEALLOCATE (rho2)

      ELSE

! recalculate Coulomb part for jellium

        IF (nion2 == 0) THEN
          IF (tcoulfalr) THEN
            CALL solv_poisson_f(rho(1:kdfull2), chpcoul, kdfull2)
          ELSE
            CALL solv_poisson_e(rho(1:kdfull2), chpcoul, kdfull2)
          END IF
        END IF

        fac = 1D0/npartto

        IF(directenergy) THEN
          encadd = fac*SUM(rho(1:kdfull2)*chpcoul(1:kdfull2))*dvol/2D0
          IF(ttest) WRITE(*,*) ' enerpw,encadd=',enerpw,encadd,enerpw - encadd
          enerpw = enerpw - encadd
        END IF

        DO ind = 1, SIZE
          aloc(ind) = aloc(ind) - fac*chpcoul(ind)
        END DO

      END IF

      RETURN
#ifdef REALSWITCH
      END SUBROUTINE calc_adsicr
#else
    END SUBROUTINE calc_adsic
#endif

! ******************************

#ifdef REALSWITCH
    SUBROUTINE calc_slaterr(rho, aloc, q0)
#else
      SUBROUTINE calc_slater(rho, aloc, q0)
#endif

! ******************************

! Computes SIC-Slater corrective potentials (one for each s.p.state)

! Input:
! rho = local density
! aloc = local KS potential
! q0 = set of s.p. wavefunctions
! Output:
! set of s.p. corrective potential 'usicsp' via module 'params'

        USE params
        USE util, ONLY: prifld2
        USE kinetic
        IMPLICIT NONE

#ifdef REALSWITCH
        REAL(DP) :: q0(kdfull2, kstate)
#else
        COMPLEX(DP) :: q0(kdfull2, kstate)
#endif
        REAL(DP) :: aloc(2*kdfull2), rho(2*kdfull2)

        REAL(DP), DIMENSION(:), ALLOCATABLE :: usicsp, rhosp
        REAL(DP), DIMENSION(:), ALLOCATABLE :: rhospu, rhospd
        LOGICAL :: testprint
        INTEGER :: ind, idx, ishift, nb, npartdw, npartup, npartto
        REAL(DP) :: enrearsave, enerpwsave, enrear1, enrear2, enpw1, enpw2, encadd, reldw, relup

        DATA testprint/.false./

!-------------------------------------------------------------------

        IF (numspin .NE. 2) STOP ' SIC-Slater requires full spin'
        IF (ifsicp /= 3) STOP ' CALC_SLATER called with wrong option IFSICP'

        CALL act_part_num(npartup, npartdw, npartto)

        IF (testprint .AND. ifsicp == 3) THEN
          OPEN (11, FILE='testslater')
          CALL prifld2(11, aloc, 'pot before u')
          CALL prifld2(11, aloc(1 + nxyz), 'pot before d')
        END IF

        ALLOCATE (rhosp(2*kdfull2))
        ALLOCATE (usicsp(2*kdfull2))
        ALLOCATE (rhospu(2*kdfull2), rhospd(2*kdfull2))

        enrearsave = enrear
        enerpwsave = enerpw
        enrear1 = 0D0
        enrear2 = 0D0
        enpw1 = 0D0
        enpw2 = 0D0
        encadd = 0D0

! recombine to spin up and spin down densities

        IF (npartup > 0) THEN
          DO ind = 1, nxyz
            rhospu(ind) = MAX(rho(ind)*(1D0 + rho(ind + nxyz))*0.5D0, small)
          END DO
        END IF
        IF (npartdw > 0) THEN
          DO ind = 1, nxyz
            rhospd(ind) = MAX(rho(ind)*(1D0 - rho(ind + nxyz))*0.5D0, small)
          END DO
        END IF

! loop over s.p. states, compute and accumulate SIC-Slater potential

        DO nb = 1, nstate
          IF (occup(nb) > small) THEN
            ishift = (ispin(nrel2abs(nb)) - 1)*nxyz ! store spin=2 in upper BLOCK
#ifdef REALSWITCH
            CALL calc_sicspr(rhosp, usicsp, q0(1, nb), nb)
#else
            CALL calc_sicsp(rhosp, usicsp, q0(1, nb), nb)
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
                relup = occup(nb)*rhosp(ind)/MAX(rhospu(ind), small)
                aloc(ind) = aloc(ind) - relup*usicsp(ind)
              END DO
            ELSE IF (ispin(nrel2abs(nb)) == 2) THEN
              DO ind = 1, nxyz
                idx = ind + nxyz
                reldw = occup(nb)*rhosp(ind)/MAX(rhospd(ind), small)
                aloc(idx) = aloc(idx) - reldw*usicsp(idx)
              END DO
            END IF
          END IF
        END DO
        encadd = encadd/2D0
        enrear = enrearsave - enrear1 - enrear2
        IF (directenergy) THEN
          enerpw = enerpwsave - enpw1 - enpw2 - encadd
        END IF
        DEALLOCATE (rhosp)
        DEALLOCATE (usicsp)
        DEALLOCATE (rhospu, rhospd)

        RETURN
#ifdef REALSWITCH
      END SUBROUTINE calc_slaterr
#else
    END SUBROUTINE calc_slater
#endif

! ******************************

#ifdef REALSWITCH
    SUBROUTINE calc_sicklir(rho, aloc, q0)
#else
      SUBROUTINE calc_sickli(rho, aloc, q0)
#endif

! Computes SIC-KLI correctiion to KS potential

! Input:
! rho = local density
! q0 = set of s.p. wavefunctions
! Input/Output:
! aloc = local KS potential

        USE params
        USE util, ONLY: prifld
        USE kinetic
        IMPLICIT NONE
        REAL(DP), PARAMETER :: sgnkli = 1D0 ! -1D0


#ifdef REALSWITCH
        REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
#else
        COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
#endif
        REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2), rho(2*kdfull2)

        REAL(DP), DIMENSION(:), ALLOCATABLE :: usicsp, rhosp, rho1, rho2
        REAL(DP), DIMENSION(:), ALLOCATABLE :: rhospu, rhospd
        REAL(DP), ALLOCATABLE :: rhokli(:, :), uslater(:), ukli(:)
        LOGICAL, PARAMETER :: testprint = .false.
        INTEGER :: itmaxkli = 200
        LOGICAL:: converged = .false.

        INTEGER :: ind, idx, ishift, itkli, nb, nfdw, nfup, npartup, npartdw, npartto
        REAL(DP) :: addn, addo, correct, efermup, efermdw
        REAL(DP) :: enrearsave, enerpwsave, enrear1, enrear2, enpw1, enpw2, encadd, reldw, relup
        REAL(DP) :: acc, avuslatdw, avuslatup, sumsldw, sumslup, sumdw, sumup
!--------------------------------------------------------------------

        IF (numspin .NE. 2) STOP ' SIC-KLI requires full spin code'
        IF (ifsicp /= 4) STOP ' CALC_SICKLI called with wrong option IFSICP'

        ALLOCATE (rhokli(kdfull2, kstate), uslater(2*kdfull2), ukli(2*kdfull2))

        CALL act_part_num(npartup, npartdw, npartto)

! determine Fermi levels

        IF (temp /= 0D0) STOP 'KLI not yet compatible with temperature'
        nfup = 0
        nfdw = 0
        IF (sgnkli == -1D0) THEN
          efermup = -1000D0
          efermdw = -1000D0
          DO nb = 1, nstate
            IF (ispin(nrel2abs(nb)) == 1) THEN ! M : 1=up 2=down
              IF (amoy(nb) > efermup .AND. occup(nb) > 0.5D0) THEN
                nfup = nb
                efermup = amoy(nb)
              END IF
            ELSE
              IF (amoy(nb) > efermdw .AND. occup(nb) > 0.5D0) THEN
                nfdw = nb
                efermdw = amoy(nb)
              END IF
            END IF
          END DO
          IF (nfup == 0) STOP ' Fermi state spin-up not found'
          IF (nfdw == 0) STOP ' Fermi state spin-down not found'
        END IF

! check workspace

        enrearsave = enrear
        enerpwsave = enerpw
        enrear1 = 0D0
        enrear2 = 0D0
        enpw1 = 0D0
        enpw2 = 0D0
        encadd = 0D0
        ALLOCATE (rhospu(2*kdfull2), rhospd(2*kdfull2))

! spin up and spin down densities,

        IF (npartup > 0) THEN
          DO ind = 1, nxyz
            rhospu(ind) = MAX(rho(ind)*(1D0 + rho(ind + nxyz))*0.5D0, small)
          END DO
        END IF
        IF (npartdw > 0) THEN
          DO ind = 1, nxyz
            rhospd(ind) = MAX(rho(ind)*(1D0 - rho(ind + nxyz))*0.5D0, small)
          END DO
        END IF

! reset accumulator

        DO ind = 1, 2*nxyz
          uslater(ind) = 0D0
        END DO
        sumslup = 0D0
        sumsldw = 0D0

! loop over s.p. states
        ALLOCATE (rhosp(2*kdfull2))
        ALLOCATE (usicsp(2*kdfull2))

        DO nb = 1, nstate
          IF (occup(nb) > small) THEN
            ishift = (ispin(nrel2abs(nb)) - 1)*nxyz ! store spin=2 in upper BLOCK
#ifdef REALSWITCH
            CALL calc_sicspr(rhosp, usicsp, q0(1, nb), nb)
#else
            CALL calc_sicsp(rhosp, usicsp, q0(1, nb), nb)
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
              acc = 0D0
              DO ind = 1, nxyz
                relup = rhosp(ind)/MAX(rhospu(ind), small)
                rhokli(ind, nb) = relup
                uslater(ind) = uslater(ind) + relup*usicsp(ind)
                acc = rhosp(ind)*usicsp(ind) + acc
              END DO
              IF (nb == nfup) avuslatup = acc*dvol
              acc = acc*dvol*sgnkli
              DO ind = 1, nxyz
                uslater(ind) = uslater(ind) + acc*rhokli(ind, nb)
              END DO
              sumslup = sumslup + acc
            ELSE IF (ispin(nrel2abs(nb)) == 2) THEN
              acc = 0D0
              DO ind = 1, nxyz
                idx = ind + nxyz
                reldw = rhosp(ind)/MAX(rhospd(ind), small)
                rhokli(ind, nb) = reldw
                uslater(idx) = uslater(idx) + reldw*usicsp(idx)
                acc = rhosp(ind)*usicsp(idx) + acc
              END DO
              IF (nb == nfdw) avuslatdw = acc*dvol
              acc = acc*dvol*sgnkli
              DO ind = 1, nxyz
                idx = ind + nxyz
                uslater(idx) = uslater(idx) + acc*rhokli(ind, nb)
              END DO
              sumsldw = sumsldw + acc
            END IF
          END IF
        END DO
        encadd = encadd/2D0
        enrear = enrearsave - enrear1 - enrear2
        IF (directenergy) THEN
          enerpw = enerpwsave - enpw1 - enpw2 - encadd
        END IF

        DEALLOCATE (rhosp)
        DEALLOCATE (usicsp)

! the first loop finished - Slater potential now on 'uslater'

! the KLI iteration

        IF (testprint) WRITE (6, '(a,2(1pg12.4))') ' sumslp=', sumslup, sumsldw
        sumslup = 0D0
        sumsldw = 0D0
        DO ind = 1, 2*nxyz
          ukli(ind) = uslater(ind)
        END DO

        ALLOCATE (rho1(kdfull2))
        ALLOCATE (rho2(kdfull2))

        DO itkli = 1, itmaxkli
          DO ind = 1, nxyz
            rho1(ind) = 0D0
            rho2(ind) = 0D0
          END DO
          DO nb = 1, nstate
            IF (occup(nb) > small) THEN
              ishift = (ispin(nrel2abs(nb)) - 1)*nxyz
              acc = 0D0
              IF (ishift <= 0) THEN
                DO ind = 1, nxyz
                  acc = ukli(ind)*rhokli(ind, nb)*rhospu(ind) + acc
                END DO
                acc = acc*dvol*sgnkli
                DO ind = 1, nxyz
                  rho1(ind) = rhokli(ind, nb)*acc + rho1(ind)
                END DO
              ELSE
                DO ind = 1, nxyz
                  idx = ind + ishift
                  acc = ukli(idx)*rhokli(ind, nb)*rhospd(ind) + acc
                END DO
                acc = acc*dvol*sgnkli
                DO ind = 1, nxyz
                  rho2(ind) = rhokli(ind, nb)*acc + rho2(ind)
                END DO
              END IF
            END IF
          END DO
          addn = 0.3D0
          addo = 1D0 - addn
          sumup = 0D0
          sumdw = 0D0
          DO ind = 1, nxyz
            idx = ind + nxyz
            ukli(ind) = addo*ukli(ind) + addn*(uslater(ind) - rho1(ind))
            sumup = (ukli(ind) - (uslater(ind) - rho1(ind)))**2*rhospu(ind) + sumup
            ukli(idx) = addo*ukli(idx) + addn*(uslater(idx) - rho2(ind))
            sumdw = (ukli(idx) - (uslater(idx) - rho2(ind)))**2*rhospd(ind) + sumdw
          END DO
          sumup = dvol*sumup
          sumdw = dvol*sumdw
          IF (sgnkli == -1D0) THEN ! asymptotic correction
            acc = 0D0
            DO ind = 1, nxyz
              acc = rhospu(ind)*rhokli(ind, nfup)*ukli(ind) + acc
            END DO
            acc = acc*dvol
            correct = acc - avuslatup
            DO ind = 1, nxyz
              ukli(ind) = ukli(ind) - correct
            END DO
            acc = 0D0
            DO ind = 1, nxyz
              idx = ind + nxyz
              acc = rhospd(ind)*rhokli(ind, nfdw)*ukli(idx) + acc
            END DO
            acc = acc*dvol
            correct = acc - avuslatdw
            DO ind = nxyz + 1, 2*nxyz
              ukli(ind) = ukli(ind) - correct
            END DO
          END IF
          IF (testprint) WRITE (6, '(a,i4,3(1pg12.4))') ' itkli,sumup,sumdw=', itkli, &
            sumup, sumdw, epsoro
          IF (sumup + sumdw < epsoro**2*2D0) THEN
            converged = .true.
            EXIT
          END IF
          sumslup = sumup
          sumsldw = sumdw
        END DO

        IF (.NOT. converged) WRITE (6, '(a,2(1pg12.4))') &
          ' KLI not converged: errors=', sumup - sumslup, sumdw - sumsldw

        WRITE (6, '(a,i5,4(1pg12.4))') &
          ' KLI itkli,sums=', itkli, sumup, sumslup, sumdw, sumsldw

        DEALLOCATE (rho1)
        DEALLOCATE (rho2)
! KLI iteration terminated - modify local potential

        DO ind = 1, nxyz
          idx = ind + nxyz
          aloc(ind) = aloc(ind) - ukli(ind)
          aloc(idx) = aloc(idx) - ukli(idx)
        END DO
        IF (testprint) CALL prifld(ukli, ' SIC-KLI ')

        DEALLOCATE (rhokli, uslater, ukli)

        RETURN
#ifdef REALSWITCH
      END SUBROUTINE calc_sicklir
#else
    END SUBROUTINE calc_sickli
#endif

#ifdef REALSWITCH

! ******************************

    SUBROUTINE act_part_num(npartup, npartdw, npartto)

! Computes actual particle NUMBER for spin-up and spin-down.
!
! Input:
! occup = array of occupation numbers via module 'params'
! Output:
! npartup = nr. of spin up particles
! npartdw = nr. of spin down particles
! npartto = total nr. of particles

      USE params
      USE kinetic
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: npartup
      INTEGER, INTENT(OUT) :: npartdw
      INTEGER, INTENT(OUT) :: npartto

      INTEGER :: nb
      REAL(DP) :: partup, partdw

! compute NUMBER of states

      partup = 0D0
      partdw = 0D0
      DO nb = 1, nstate
        IF (ispin(nrel2abs(nb)) == 1) THEN ! M : 1=up 2=down
          partup = occup(nb) + partup
        ELSE
          partdw = occup(nb) + partdw
        END IF
      END DO


      npartup = NINT(partup)
      npartdw = NINT(partdw)

      npartto = npartup + npartdw

      RETURN
    END SUBROUTINE act_part_num
#endif

! ******************************

#ifdef REALSWITCH
    SUBROUTINE calc_sicspr(rhosp, usicsp, q0state, nb)
#else
      SUBROUTINE calc_sicsp(rhosp, usicsp, q0state, nb)
#endif
! ******************************

! computes SIC potential for state 'nb'.
! input is
! q0state = wavefunction for s.p. state
! nb = NUMBER of state
! output is
! usicsp = the s.p. SIC potential
! output via module 'params'
! enrear,enerpw = rearrangement energies for 'nb'
! encoulsp = Coulomb energy for 'nb'

        USE params
        USE util, ONLY: prifld2
        USE kinetic
        USE coulsolv, ONLY: solv_poisson_f, solv_poisson_e, tcoulfalr

        IMPLICIT NONE


#ifdef REALSWITCH
        REAL(DP), INTENT(IN) :: q0state(kdfull2)
#else
        COMPLEX(DP), INTENT(IN) :: q0state(kdfull2)
#endif
        REAL(DP), INTENT(IN OUT) :: rhosp(2*kdfull2)
        REAL(DP), INTENT(IN OUT) :: usicsp(2*kdfull2)

        REAL(DP), DIMENSION(:), ALLOCATABLE :: chpdftsp, couldif, rho1

        LOGICAL :: testprint = .false.
        INTEGER :: ind, idx, ishift, nb
        REAL(DP) :: enrearsave, enerpwsave, enrear1, enrear2, enpw1, enpw2, encadd, encsum
!--------------------------------------------------------------------

        IF (numspin .NE. 2) STOP ' CALC_SICSP requires fullspin code'

        ALLOCATE (chpdftsp(2*kdfull2))
        ALLOCATE (couldif(kdfull2))
        ALLOCATE (rho1(kdfull2))

        enrearsave = enrear
        enerpwsave = enerpw
        enrear1 = 0D0
        enrear2 = 0D0
        enpw1 = 0D0
        enpw2 = 0D0
        encadd = 0D0

        IF (occup(nb) > small) THEN
          ishift = (ispin(nrel2abs(nb)) - 1)*nxyz ! store spin=2 in upper block

! density for s.p. state

          DO ind = 1, nxyz
#ifdef REALSWITCH
            rhosp(ind) = q0state(ind)*q0state(ind)
#else
            rhosp(ind) = REAL(CONJG(q0state(ind))*q0state(ind), DP)
#endif
            rhosp(ind + nxyz) = 3 - 2*ispin(nrel2abs(nb)) ! M : 1 if spinup -1 if spin down
            rho1(ind) = rhosp(ind)
          END DO
! DFT for s.p. state

          CALL calc_lda(rhosp, chpdftsp) ! --> enrear,enerpw

          IF (ispin(nrel2abs(nb)) == 1) THEN
            DO ind = 1, nxyz
              usicsp(ind) = chpdftsp(ind)
            END DO
          ELSE
            DO ind = 1, nxyz
              idx = ind + nxyz
              usicsp(idx) = chpdftsp(idx)
            END DO
          END IF

! s.p. Coulomb and subtract from LDA

          DO ind = 1, nxyz
            rho1(ind) = rhosp(ind)
          END DO
          IF (tcoulfalr) THEN
            CALL solv_poisson_f(rho1, couldif, kdfull2)
          ELSE
            CALL solv_poisson_e(rho1, couldif, kdfull2)
          END IF

          IF (testprint) THEN
            WRITE (11, '(a)') ' ', ' '
            WRITE (11, '(a,i3)') '# NB=', nb
            CALL prifld2(11, rhosp, ' density')
            CALL prifld2(11, couldif, ' Coulomb')
            CALL prifld2(11, chpdftsp, ' P&W')
          END IF

          encsum = 0D0
          IF (ispin(nrel2abs(nb)) == 1) THEN
            DO ind = 1, nxyz
              usicsp(ind) = usicsp(ind) + couldif(ind)
              IF (directenergy) THEN
                encsum = encsum + rhosp(ind)*couldif(ind)
              END IF
            END DO
            IF (testprint) CALL prifld2(11, usicsp, ' SIC pot')
          ELSE
            DO ind = 1, nxyz
              idx = ind + nxyz
              usicsp(idx) = usicsp(idx) + couldif(ind)
              IF (directenergy) THEN
                encsum = encsum + rhosp(ind)*couldif(ind)
              END IF
            END DO
            IF (testprint) CALL prifld2(11, usicsp(nxyz + 1), ' SIC pot')
          END IF
          encoulsp = encsum*dvol/2D0
        ELSE
          usicsp(1:nxyz) = 0D0
          encoulsp = 0D0
        END IF
        DEALLOCATE (chpdftsp)
        DEALLOCATE (couldif)
        DEALLOCATE (rho1)

        RETURN
#ifdef REALSWITCH
      END SUBROUTINE calc_sicspr
#else
    END SUBROUTINE calc_sicsp
#endif

! ******************************

#ifdef REALSWITCH
    SUBROUTINE exchgr(q0, qex, nbel)
#else
      SUBROUTINE exchg(q0, qex, nbel)
#endif

! ******************************

! exact exchange:
! input for static version is set of wavefunctions on 'q0'
! input for the dynamic version are
! the wavefunctions 'q0' on which the exchange acts
! and the wf's 'psisavex' for the density matrix
! which are communicated via module PARAMS
! output are accumulated exchange wavefunctions on 'qex'
! only the exchange for the one state 'nbe' is handled
! in case of COMPLEX wavefunctions

        USE params
        USE coulsolv, ONLY: solv_poisson_f, solv_poisson_e, tcoulfalr

        IMPLICIT NONE

#ifdef REALSWITCH
        REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
        REAL(DP), INTENT(OUT) :: qex(kdfull2, kstate)
        INTEGER, INTENT(IN), OPTIONAL :: nbel
        INTEGER :: i
        REAL(DP) :: sump
#else
        COMPLEX(DP), INTENT(IN) :: q0(kdfull2)
        COMPLEX(DP), INTENT(OUT) :: qex(kdfull2)
        INTEGER, INTENT(IN) :: nbel
        REAL(DP), DIMENSION(:), ALLOCATABLE :: acli
        REAL(DP), DIMENSION(:), ALLOCATABLE :: rhi
        COMPLEX(DP) :: rhoc
#endif

! workspaces

        REAL(DP), DIMENSION(:), ALLOCATABLE :: rh
        REAL(DP), DIMENSION(:), ALLOCATABLE :: acl

        LOGICAL, PARAMETER :: ttest = .FALSE.
        INTEGER :: ind, nb2, nbe, nlow, nup

        IF (ifsicp /= 5) STOP ' in EXCHANGE: wrong option IFSICP'

        ALLOCATE (acl(2*kdfull2))
        ALLOCATE (rh(2*kdfull2))

#ifdef COMPLEXSWITCH
        ALLOCATE (rhi(2*kdfull2))
        ALLOCATE (acli(2*kdfull2))
        nbe = nbel
#endif

        qex = 0D0

! DOUBLE loop over states:
! outer loop for states on which exchange acts
! inner loop for states with which exchange acts

        dvol = dx*dy*dz
#ifdef REALSWITCH
        IF(PRESENT(nbel)) THEN
          nlow = nbel
          nup = nbel
        ELSE
          nlow = 1
          nup = nstate
        END IF
        DO nbe = 1, nstate
#endif
          DO nb2 = 1, nstate
            IF (ispin(nrel2abs(nbe)) == ispin(nrel2abs(nb2)) &
                .AND. occup(nrel2abs(nb2)) > 0.5D0) THEN
              IF (ttest) WRITE (*, *) ' NBE,NB2,SPINS,OCCUP:', &
                nbe, nb2, ispin(nrel2abs(nbe)), ispin(nrel2abs(nb2)), &
                occup(nrel2abs(nb2))

! compute transition density

              DO ind = 1, nxyz
#ifdef REALSWITCH
                rh(ind) = q0(ind, nb2)*q0(ind, nbe)
#else
                rhoc = CONJG(psisavex(ind, nb2))*q0(ind)
                rh(ind) = REAL(rhoc, DP)
                rhi(ind) = AIMAG(rhoc)
#endif
              END DO

! the Coulomb potential for the transition density
! (warning : counet inserts the esquar factor)

              IF (tcoulfalr) THEN
                CALL solv_poisson_f(rh, acl, kdfull2)
#ifdef COMPLEXSWITCH
                CALL solv_poisson_f(rhi, acli, kdfull2)
#endif
              ELSE
                CALL solv_poisson_e(rh, acl, kdfull2)
#ifdef COMPLEXSWITCH
                CALL solv_poisson_e(rhi, acli, kdfull2)
#endif
              END IF
! accumulate on wavefunction

#ifdef REALSWITCH
              IF (ttest) THEN
                sump = 0D0
                DO i = 1, nxyz
                  sump = q0(i, nbe)*q0(i, nb2)*acl(i) + sump
                END DO
                WRITE (6, '(a,2i5,1pg12.4)') &
                  ' EXCHANGE: nbe,nb2,overlap=', nbe, nb2, sump*dvol
              END IF
#endif

#ifdef REALSWITCH
              qex(:, nbe) = qex(:, nbe) - q0(:, nb2)*acl(1:kdfull2)
#else
              qex(:) = qex(:) - psisavex(:, nb2)*CMPLX(acl(:), acli(:), DP)
#endif

#ifdef REALSWITCH
              IF (ttest) THEN
                sump = 0D0
                DO i = 1, nxyz
                  sump = q0(i, nbe)*qex(i, nbe) + sump
                END DO
                WRITE (6, '(a,2i5,1pg12.4)') &
                  ' EXCHANGE: nbe,nb2,overlap=', nbe, nb2, sump*dvol
              END IF
#endif

            END IF
          END DO
#ifdef REALSWITCH
        END DO
#endif

        DEALLOCATE (acl)
        DEALLOCATE (rh)

#ifdef COMPLEXSWITCH
        DEALLOCATE (acli)
        DEALLOCATE (rhi)
#endif

        RETURN
#ifdef REALSWITCH
      END SUBROUTINE exchgr
#else
    END SUBROUTINE exchg
#endif

