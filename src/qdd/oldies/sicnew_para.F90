
! ******************************

#ifdef REALSWITCH
SUBROUTINE calc_sicr(rho, aloc, q0)
#else
  SUBROUTINE calc_sic(rho, aloc, q0)
#endif

! ******************************

! switchboard for the various brands of SIC

    USE params
    USE kinetic
    IMPLICIT REAL(DP) (A - H, O - Z)!
#ifdef REALSWITCH
    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
#else
    COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
#endif
!REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2),rho(2*kdfull2)
    REAL(DP), INTENT(IN OUT) :: aloc(*), rho(*)

!------------------------------------------------------------

    IF (ifreezekspot .AND. tfs > 0) RETURN

#ifdef REALSWITCH

    IF (ifsicp == 2) THEN
      CALL calc_adsicr(rho, aloc, q0)
    ELSE IF (ifsicp == 3) THEN
      CALL calc_slaterr(rho, aloc, q0)
    ELSE IF (ifsicp == 4) THEN
      CALL calc_sicklir(rho, aloc, q0)
    ELSE
      RETURN
    END IF
#else
    IF (ifsicp == 2) THEN
      CALL calc_adsic(rho, aloc, q0)
    ELSE IF (ifsicp == 3) THEN
      CALL calc_slater(rho, aloc, q0)
    ELSE IF (ifsicp == 4) THEN
      CALL calc_sickli(rho, aloc, q0)
    ELSE
      RETURN
    END IF
#endif
    RETURN
#ifdef REALSWITCH
  END SUBROUTINE calc_sicr
#else
END SUBROUTINE calc_sic
#endif

! ******************************

#ifdef REALSWITCH
SUBROUTINE calc_adsicr(rho, aloc, q0)
#else
  SUBROUTINE calc_adsic(rho, aloc, q0)
#endif

! ******************************

! computes local part of Hamiltonian
! 'time' <= 0 signal static iteration i.e. without laser field

    USE params
    USE kinetic
#if(netlib_fft|fftw_cpu)
    USE coulsolv
#endif
    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), PARAMETER :: sgnkli = -1D0

#if(mpi)
    INCLUDE 'mpif.h'
    INTEGER :: is(mpi_status_size)
#endif

#ifdef REALSWITCH
    REAL(DP) :: q0(kdfull2, kstate)
#else
    COMPLEX(DP) :: q0(kdfull2, kstate)
#endif
    REAL(DP) :: aloc(2*kdfull2), rho(2*kdfull2)
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rhosp, chpdftsp, coulsum, couldif
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rho1, rho2
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rhospu, rhospd, chpdftspu, chpdftspd
!EQUIVALENCE (coulsum,w1),(couldif,w2)
!EQUIVALENCE (rhosp,w1)
!EQUIVALENCE (chpdftsp,w3)
!EQUIVALENCE (rho2,rhosp)
    INTEGER :: SIZE

#if(mpi)
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rhonod1, rhonod2, rhotp
    REAL(DP), DIMENSION(:), ALLOCATABLE :: tp1, tp2
    REAL(DP), DIMENSION(:), ALLOCATABLE :: aloc1, aloc2
    INTEGER :: nod
#endif

!-------------------------------------------------------------------------

    IF (ifsicp /= 2) THEN
      STOP ' CALC_SIC called with wrong option IFSICP'
    END IF

    CALL act_part_num(npartup, npartdw, npartto)

#if(mpi)
    CALL mpi_comm_rank(mpi_comm_world, nod, icode)
    SIZE = lengnod(nod + 1)
    ALLOCATE (rhonod1(SIZE))
    ALLOCATE (rhonod2(SIZE))
    CALL pi_scatterv(rho, nxyz, rhonod1, SIZE, icode)
    CALL pi_scatterv(rho(nxyz + 1), nxyz, rhonod2, SIZE, icode)
    ALLOCATE (tp1(SIZE))
    ALLOCATE (tp2(SIZE))
#else
    SIZE = nxyz
#endif

    IF (numspin == 2) ALLOCATE (rhospu(2*kdfull2), rhospd(2*kdfull2), &
                                chpdftspu(2*kdfull2), chpdftspd(2*kdfull2))

    enrearsave = enrear
    enerpwsave = enerpw
    enrear1 = 0D0
    enrear2 = 0D0
    enpw1 = 0D0
    enpw2 = 0D0
    encadd = 0D0
    IF (numspin == 2) THEN

! averaged spinup and spindown density

      IF (npartup > 0) THEN
        facuph = 0.5/npartup

! DO ind=1,nxyz
        DO ind = 1, SIZE
#if(mpi)
          tp1(ind) = rhonod1(ind)*(1.0 + rhonod2(ind))*facuph
          tp2(ind) = 1.0
#else
          rhospu(ind) = rho(ind)*(1.0 + rho(ind + nxyz))*facuph
          rhospu(ind + nxyz) = 1.0
#endif
        END DO
      END IF
#if(mpi)
      CALL pi_allgatherv(tp1, SIZE, rhospu, nxyz, icode)
      CALL pi_allgatherv(tp2, SIZE, rhospu(nxyz + 1), nxyz, icode)
#endif

      IF (npartdw > 0) THEN
        facdwh = 0.5/npartdw

! DO ind=1,nxyz
        DO ind = 1, SIZE
#if(mpi)
          tp1(ind) = rhonod1(ind)*(1.0 - rhonod2(ind))*facuph
          tp2(ind) = -1.0
#else
          rhospd(ind) = rho(ind)*(1.0 - rho(ind + nxyz))*facdwh
          rhospd(ind + nxyz) = -1.0
#endif
        END DO
      END IF
#if(mpi)
      CALL pi_allgatherv(tp1, SIZE, rhospd, nxyz, icode)
      CALL pi_allgatherv(tp2, SIZE, rhospd(nxyz + 1), nxyz, icode)
#endif

! DFT for averaged spinup and spindown density

      IF (npartup > 0) THEN
        CALL calc_lda(rhospu, chpdftspu)
        enrear1 = enrear
      END IF
      IF (npartdw > 0) THEN
        CALL calc_lda(rhospd, chpdftspd)
        enrear2 = enrear
      END IF

      enrear = enrearsave - enrear1*npartup - enrear2*npartdw
! WRITE(6,*) ' enrear.s=',enrearsave, enrear1,enrear2,enrear

#if(mpi)
      ALLOCATE (aloc1(SIZE))
      ALLOCATE (aloc2(SIZE))
      CALL pi_scatterv(chpdftspu, nxyz, tp1, SIZE, icode)
      CALL pi_scatterv(chpdftspd(nxyz + 1), nxyz, tp2, SIZE, icode)
      CALL pi_scatterv(aloc, nxyz, aloc1, SIZE, icode)
      CALL pi_scatterv(aloc(nxyz + 1), nxyz, aloc2, SIZE, icode)
#endif

      IF (npartup > 0) THEN
! DO ind=1,nxyz
        DO ind = 1, SIZE
#if(mpi)
          aloc1(ind) = aloc1(ind) - tp1(ind)
#else
          aloc(ind) = aloc(ind) - chpdftspu(ind)
#endif
        END DO
#if(mpi)
        CALL pi_allgatherv(aloc1, SIZE, aloc, nxyz, icode)
#endif
      END IF

      IF (npartdw > 0) THEN
#if(mpi)
        DO idx = 1, SIZE
          aloc2(idx) = aloc2(idx) - tp2(idx)
        END DO
#else
        DO idx = nxyz + 1, 2*nxyz
          aloc(idx) = aloc(idx) - chpdftspd(idx)
        END DO
#endif
#if(mpi)
        CALL pi_allgatherv(aloc2, SIZE, aloc(nxyz + 1), nxyz, icode)
#endif
      END IF

      DEALLOCATE (rhospu, rhospd, chpdftspu, chpdftspd)

    ELSE

      ALLOCATE (chpdftsp(2*kdfull2))
      ALLOCATE (rhosp(2*kdfull2))

      factotal = 1.0/npartto

#if(mpi)
      ALLOCATE (tp1(SIZE))
      ALLOCATE (tp2(SIZE))
#endif

! DO ind=1,nxyz
      DO ind = 1, SIZE
#if(mpi)
        tp1(ind) = rhonod1(ind)*factotal
        tp2(ind) = 1.0
#else
        rhosp(ind) = rho(ind)*factotal
        rhosp(ind + nxyz) = 1.0
#endif
      END DO

#if(mpi)
      CALL pi_allgatherv(tp1, SIZE, rhosp, nxyz, icode)
      CALL pi_allgatherv(tp2, SIZE, rhosp(nxyz + 1), nxyz, icode)
#endif

! DFT for averaged s.p. state

      CALL calc_lda(rhosp, chpdftsp)

! subtract from pure LDA potential

#if(mpi)
      CALL pi_scatterv(chpdftsp, nxyz, tp1, SIZE, icode)
      CALL pi_scatterv(aloc, nxyz, aloc1, SIZE, icode)
#endif
! DO ind=1,nxyz
      DO ind = 1, SIZE
#if(mpi)
        aloc1(ind) = aloc1(ind) - tp1(ind)
#else
        aloc(ind) = aloc(ind) - chpdftsp(ind)
#endif
      END DO
#if(mpi)
      CALL pi_allgatherv(aloc1, SIZE, aloc, nxyz, icode)
#endif

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

! DO ind=1,nxyz
      DO ind = 1, SIZE
#if(mpi)
        tp2(ind) = rhonod1(ind)
        tp1(ind) = rhonod1(ind)*rhonod2(ind)
#else
        rho2(ind) = rho(ind)
        rho1(ind) = rho(ind)*rho(ind + nxyz)
#endif
      END DO
#if(mpi)
      CALL pi_allgatherv(tp1, SIZE, rho1, nxyz, icode)
      CALL pi_allgatherv(tp2, SIZE, rho2, nxyz, icode)
#endif

#if(gridfft)
      CALL falr(rho1(1), couldif, nx2, ny2, nz2, kdfull2)
      IF (idielec == 0 .AND. nion2 > 0) THEN
        coulsum = chpcoul
      ELSE
        CALL falr(rho2(1), coulsum, nx2, ny2, nz2, kdfull2)
      END IF
#endif
#if(findiff|numerov)
      CALL solv_fft(rho1(1), couldif, dx, dy, dz)
      CALL solv_fft(rho2(1), coulsum, dx, dy, dz)
#endif
      facup = 1.0/npartup

#if(mpi)
      CALL pi_scatterv(coulsum, nxyz, tp1, SIZE, icode)
      CALL pi_scatterv(couldif, nxyz, tp2, SIZE, icode)
      CALL pi_scatterv(aloc, nxyz, aloc1, SIZE, icode)
#endif
! DO ind=1,nxyz
      DO ind = 1, SIZE
#if(mpi)
        coulup = 0.5*(tp1(ind) + tp2(ind))
        aloc1(ind) = aloc1(ind) - coulup*facup
#else
        coulup = 0.5*(coulsum(ind) + couldif(ind))
        aloc(ind) = aloc(ind) - coulup*facup
#endif
      END DO
#if(mpi)
      CALL pi_allgatherv(aloc1, SIZE, aloc, nxyz, icode)
#endif

      IF (npartdw > 0) THEN

#if(mpi)
        CALL pi_scatterv(aloc(nxyz + 1), nxyz, aloc2, SIZE, icode)
#endif
        facdw = 1.0/npartdw
! DO ind=1,nxyz
        DO ind = 1, SIZE
#if(mpi)
          couldw = 0.5*(tp1(ind) - tp2(ind))
          aloc2(ind) = aloc2(ind) - facdw*couldw
#else
          couldw = 0.5*(coulsum(ind) - couldif(ind))
          idx = ind + nxyz
          aloc(idx) = aloc(idx) - facdw*couldw
#endif
        END DO
#if(mpi)
        CALL pi_allgatherv(aloc2, SIZE, aloc(nxyz + 1), nxyz, icode)
        DEALLOCATE (aloc1, aloc2, tp1, tp2)
#endif

      END IF

      DEALLOCATE (coulsum)
      DEALLOCATE (couldif)
      DEALLOCATE (rho1)
      DEALLOCATE (rho2)

    ELSE

! recalculate Coulomb part for jellium

#if(gridfft)
      IF (nion2 == 0) THEN
        CALL falr(rho(1), chpcoul, nx2, ny2, nz2, kdfull2)
      END IF
#endif
#if(findiff|numerov)
      STOP ' ADSIC not yet ready for nospin and finite diff.'
#endif

      fac = 1.0/npartto
#if(mpi)
      CALL pi_scatterv(chpcoul, nxyz, tp1, SIZE, icode)
      CALL pi_scatterv(aloc, nxyz, aloc1, SIZE, icode)
      WRITE (*, *) 'scattering of chpcoul and aloc done'
#endif
! DO ind=1,nxyz
      DO ind = 1, SIZE
#if(mpi)
        aloc1(ind) = aloc1(ind) - fac*tp1(ind)
#else
        aloc(ind) = aloc(ind) - fac*chpcoul(ind)
#endif
      END DO
#if(mpi)
      CALL pi_allgatherv(aloc1, SIZE, aloc, nxyz, icode)
      DEALLOCATE (aloc1, aloc2, tp1)
#endif

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

! computes local part of Hamiltonian
! 'time' <= 0 signal static iteration i.e. without laser field

    USE params
    USE kinetic
    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), PARAMETER :: sgnkli = -1D0
#if(mpi)
    INCLUDE 'mpif.h'
    INTEGER :: is(mpi_status_size)
#endif
#ifdef REALSWITCH
    REAL(DP) :: q0(kdfull2, kstate)
#else
    COMPLEX(DP) :: q0(kdfull2, kstate)
#endif
    REAL(DP) :: aloc(2*kdfull2), rho(2*kdfull2)

    REAL(DP), DIMENSION(:), ALLOCATABLE :: usicsp, rhosp
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rhospu, rhospd
    LOGICAL :: testprint
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

! IF(usew1) STOP ' IN CALC_SIC: workspace W1 already active '
! usew1 = .true.
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

! recombine TO spin up and spin down densities

    IF (npartup > 0) THEN
      DO ind = 1, nxyz
        rhospu(ind) = MAX(rho(ind)*(1.0 + rho(ind + nxyz))*0.5, small)
      END DO
    END IF
    IF (npartdw > 0) THEN
      DO ind = 1, nxyz
        rhospd(ind) = MAX(rho(ind)*(1.0 - rho(ind + nxyz))*0.5, small)
      END DO
    END IF

! loop over s.p. states, compute and accumulate SIC-Slater potential

    DO nb = 1, nstate
      IF (occup(nb) > small) THEN
        ishift = (ispin(nrel2abs(nb)) - 1)*nxyz ! store spin=2 IN upper BLOCK
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
    encadd = encadd/2.0
    enrear = enrearsave - enrear1 - enrear2
    IF (directenergy) THEN
      enerpw = enerpwsave - enpw1 - enpw2 - encadd
    END IF
! usew1 = .false.
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

! ******************************

! computes local part of Hamiltonian
! 'time' <= 0 signal static iteration i.e. without laser field

    USE params
    USE kinetic
    IMPLICIT REAL(DP) (A - H, O - Z)
    PARAMETER(sgnkli=-1.0)
#if(mpi)
    INCLUDE 'mpif.h'
    INTEGER :: is(mpi_status_size)
#endif

#ifdef REALSWITCH
    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
#else
    COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate)
#endif
    REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2), rho(2*kdfull2)

    REAL(DP), DIMENSION(:), ALLOCATABLE :: usicsp, rhosp, rho1, rho2
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rhospu, rhospd
    REAL(DP), ALLOCATABLE :: rhokli(:, :), uslater(:), ukli(:)
! EQUIVALENCE (rhosp,w1)
    LOGICAL, PARAMETER :: testprint = .false.
    INTEGER :: itmaxkli = 200

!--------------------------------------------------------------------

    IF (numspin .NE. 2) STOP ' SIC-KLI requires full spin code'
    IF (ifsicp /= 4) STOP ' CALC_SICKLI called with wrong option IFSICP'

    ALLOCATE (rhokli(kdfull2, kstate), uslater(2*kdfull2), ukli(2*kdfull2))

    CALL act_part_num(npartup, npartdw, npartto)

! determine Fermi levels

    IF (temp /= 0D0) STOP 'KLI not yet compatible with temperature'
    nfup = 0
    nfdw = 0
    IF (sgnkli == -1.0) THEN
      efermup = -1000D0
      efermdw = -1000D0
      DO nb = 1, nstate
        IF (ispin(nrel2abs(nb)) == 1) THEN ! M : 1=up 2=down
          IF (amoy(nb) > efermup .AND. occup(nb) > 0.5) THEN
            nfup = nb
            efermup = amoy(nb)
          END IF
        ELSE
          IF (amoy(nb) > efermdw .AND. occup(nb) > 0.5) THEN
            nfdw = nb
            efermdw = amoy(nb)
          END IF
        END IF
      END DO
      IF (nfup == 0) STOP ' Fermi state spin-up not found'
      IF (nfdw == 0) STOP ' Fermi state spin-down not found'
    END IF

! check workspace

! IF(usew1) STOP ' IN CALC_SIC: workspace W1 already active '
! usew1 = .true.

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
        rhospu(ind) = MAX(rho(ind)*(1.0 + rho(ind + nxyz))*0.5, small)
      END DO
    END IF
    IF (npartdw > 0) THEN
      DO ind = 1, nxyz
        rhospd(ind) = MAX(rho(ind)*(1.0 - rho(ind + nxyz))*0.5, small)
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
        ishift = (ispin(nrel2abs(nb)) - 1)*nxyz ! store spin=2 IN upper BLOCK
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
!old uslater(ind) = uslater(ind)-acc*rhokli(ind,nb)
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
!old uslater(idx) = uslater(idx)-acc*rhokli(ind,nb)
            uslater(idx) = uslater(idx) + acc*rhokli(ind, nb)
          END DO
          sumsldw = sumsldw + acc
        END IF
      END IF
    END DO
    encadd = encadd/2.0
    enrear = enrearsave - enrear1 - enrear2
    IF (directenergy) THEN
      enerpw = enerpwsave - enpw1 - enpw2 - encadd
    END IF

    DEALLOCATE (rhosp)
    DEALLOCATE (usicsp)
!ALLOCATE(rhospu(2*kdfull2),rhospd(2*kdfull2))

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
! w2(ind) = 0D0
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
! DO ind=1,nxyz
! w2(ind) = rhokli(ind,nb)*acc + w2(ind)
! END DO
          END IF
        END IF
      END DO
      addn = 0.3
      addo = 1.0 - addn
      sumup = 0D0
      sumdw = 0D0
      DO ind = 1, nxyz
        idx = ind + nxyz
!old ukli(ind)=uslater(ind)+w1(ind)
!old ukli(idx)=uslater(idx)+w2(ind)
        ukli(ind) = addo*ukli(ind) + addn*(uslater(ind) - rho1(ind))
        sumup = (ukli(ind) - (uslater(ind) - rho1(ind)))**2*rhospu(ind) + sumup
        ukli(idx) = addo*ukli(idx) + addn*(uslater(idx) - rho2(ind))
        sumdw = (ukli(idx) - (uslater(idx) - rho2(ind)))**2*rhospd(ind) + sumdw
!old sumup = rhospu(ind)*w1(ind) +sumup
!old sumdw = rhospd(ind)*w2(ind) +sumdw
      END DO
      sumup = dvol*sumup
      sumdw = dvol*sumdw
      IF (sgnkli == -1.0) THEN ! asymptotic correction
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
!old & sumup-sumslup,sumdw-sumsldw,epsoro
!old IF(abs(sumup-sumslup)+abs(sumdw-sumsldw).lt.epsoro) goto 99
      IF (sumup + sumdw < epsoro**2*2.0) GO TO 99
      sumslup = sumup
      sumsldw = sumdw
    END DO
    WRITE (6, '(a,2(1pg12.4))') &
      ' KLI not converged: errors=', sumup - sumslup, sumdw - sumsldw
99  CONTINUE
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

! ******************************

! computes actual particle NUMBER for spin-up and spin-down.

  USE params
  USE kinetic
  IMPLICIT REAL(DP) (A - H, O - Z)
#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
#endif

  REAL(DP) :: partup, partdw
  INTEGER :: npartup, npartdw, npartto

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

#if(mpi)
  CALL mpi_barrier(mpi_comm_world, mpi_ierror)
  CALL mpi_allreduce(partup, partupp, 1, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, ic)
  CALL mpi_allreduce(partdw, partdwp, 1, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, ic)
  CALL mpi_barrier(mpi_comm_world, mpi_ierror)
  partup = partupp
  partdw = partdwp
#endif

!mb npartup = nint(partup+0.000001) ! M : INTEGER ???
!mb npartdw = nint(partdw+0.000001)
  npartup = partup + 0.000001D0 ! M : INTEGER ???
  npartdw = partdw + 0.000001D0

  npartto = npartup + npartdw

  RETURN
END SUBROUTINE act_part_num
#endif

#ifdef COMPLEXSWITCH
SUBROUTINE calc_sicspc(rhosp, usicsp, q0state, nb)

! ******************************

! computes SIC potential for state 'nb'.
! input is
! q0state = wavefunction for s.p. state
! nb = NUMBER of state
! output is
! usicsp = the s.p. SIC potential
! output via COMMON
! enrear,enerpw = rearrangement energies for 'nb'
! encoulsp = Coulomb energy for 'nb'

!#include "all.inc"
  USE params
  USE kinetic
#if(netlib_fft|fftw_cpu)
  USE coulsolv
#endif
  IMPLICIT REAL(DP) (A - H, O - Z)

  COMPLEX(DP), INTENT(IN OUT) :: q0state(kdfull2)
  REAL(DP), INTENT(IN OUT) :: rhosp(2*kdfull2), usicsp(2*kdfull2)

! the usage of workspacss is tricky:
! 'rhosp' from the calling list may be equivalenced with 'w1', but
! occupies temporarily two workspaces,. i.e. 'w2=couldif'.
! The upper BLOCK is obsolete after 'calc_lda' as is the whole
! 'chpdftsp' after transfer TO 'uscisp'. The free slots are THEN
! used for 'rho1' and 'couldif' IN the Coulomb solver.

  REAL(DP), ALLOCATABLE :: chpdftsp(:)
  REAL(DP), ALLOCATABLE :: couldif(:)
  REAL(DP), ALLOCATABLE :: rho1(:)

  LOGICAL, PARAMETER :: testprint = .false.

!--------------------------------------------------------------------

  IF (numspin < 2) STOP ' CALC_SICSP requires fullspin code'

! ALLOCATE workspace

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
! WRITE(*,*) ' CALC_SICSPC: nb,occup(nb)',nb,occup(nb)

  IF (occup(nb) > small) THEN
    ishift = (ispin(nrel2abs(nb)) - 1)*nxyz ! store spin=2 IN upper BLOCK

! density for s.p. state

    DO ind = 1, nxyz
      rhosp(ind) = REAL(CONJG(q0state(ind))*q0state(ind))
      rhosp(ind + nxyz) = 3 - 2*ispin(nrel2abs(nb)) ! M : 1 IF spinup -1 IF spin down
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
#if(gridfft)
    CALL falr(rho1(1), couldif, nx2, ny2, nz2, kdfull2)
#endif
#if(findiff|numerov)
    CALL solv_fft(rho1(1), couldif, dx, dy, dz)
#endif
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
      encoulsp = encsum*dvol/2D0
      IF (testprint) CALL prifld2(11, usicsp, ' SIC pot')
    ELSE
      DO ind = 1, nxyz
        idx = ind + nxyz
        usicsp(idx) = usicsp(idx) + couldif(ind)
        IF (directenergy) THEN
          encsum = encsum + rhosp(ind)*couldif(ind)
        END IF
      END DO
      encoulsp = encsum*dvol/2D0
      IF (testprint) CALL prifld2(11, usicsp(nxyz + 1), ' SIC pot')
    END IF
  ELSE
    usicsp = 0D0
    encoulsp = 0D0
  END IF

!WRITE(*,*) ' SICSPC: encoulsp,enerpw=',encsum,enerpw

  DEALLOCATE (chpdftsp)
  DEALLOCATE (couldif)
  DEALLOCATE (rho1)

  RETURN
END SUBROUTINE calc_sicspc
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
! output via COMMON
! enrear,enerpw = rearrangement energies for 'nb'
! encoulsp = Coulomb energy for 'nb'

    USE params
    USE kinetic
#if(netlib_fft|fftw_cpu)
    USE coulsolv
#endif
    IMPLICIT REAL(DP) (A - H, O - Z)

#if(mpi)
    INCLUDE 'mpif.h'
    INTEGER :: is(mpi_status_size)
#endif

#ifdef REALSWITCH
    REAL(DP), INTENT(IN) :: q0state(kdfull2)
#else
    COMPLEX(DP), INTENT(IN) :: q0state(kdfull2)
#endif
    REAL(DP), INTENT(IN OUT) :: rhosp(2*kdfull2), usicsp(2*kdfull2)

    REAL(DP), DIMENSION(:), ALLOCATABLE :: chpdftsp, couldif, rho1
!DIMENSION chpdftsp(2*kdfull2)
!DIMENSION couldif(kdfull2)
!DIMENSION rho1(kdfull2)
!EQUIVALENCE (couldif,w2)
!EQUIVALENCE (chpdftsp,w3)
!EQUIVALENCE (rho1,w4)

    LOGICAL :: testprint = .false.

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
! WRITE(*,*) ' CALC_SICSPR: nb,occup(nb)',nb,occup(nb)

    IF (occup(nb) > small) THEN
      ishift = (ispin(nrel2abs(nb)) - 1)*nxyz ! store spin=2 IN upper BLOCK

! density for s.p. state

      DO ind = 1, nxyz
#ifdef REALSWITCH
        rhosp(ind) = q0state(ind)*q0state(ind)
#else
        rhosp(ind) = REAL(CONJG(q0state(ind))*q0state(ind))
#endif
        rhosp(ind + nxyz) = 3 - 2*ispin(nrel2abs(nb)) ! M : 1 IF spinup -1 IF spin down
        rho1(ind) = rhosp(ind)
      END DO
! WRITE(*,*) 'SICSP: nb,norm-dens=',nb,dvol*SUM(rhosp(1:nxyz))
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
#if(gridfft)
      CALL falr(rho1(1), couldif, nx2, ny2, nz2, kdfull2)
#endif
#if(findiff|numerov)
      CALL solv_fft(rho1(1), couldif, dx, dy, dz)
#endif
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
! WRITE(*,*) 'SICSP: encoulsp,enerpw=',encoulsp,enerpw

    RETURN
#ifdef REALSWITCH
  END SUBROUTINE calc_sicspr
#else
END SUBROUTINE calc_sicsp
#endif

! ******************************

#ifdef REALSWITCH
SUBROUTINE exchgr(q0, qex)
#else
  SUBROUTINE exchg(q0, qex, nbe)
#endif

! ******************************

! exact exchange:
! input for static version is set of wavefunctions on 'q0'
! input for the dynamic version are
! the wavefunctions 'q0' on which the exchange acts
! and the wf's 'psisavex' for the density matrix
! which are communicated via MODULE PARAMS
! output are accumulated exchange wavefunctions on 'qex'
! ONLY the exchange for the one state 'nbe' is handled
! IN CASE of COMPLEX wavefunctions

    USE params
!USE kinetic
#if(netlib_fft|fftw_cpu)
    USE coulsolv
#endif
    IMPLICIT REAL(DP) (A - H, O - Z)

#ifdef REALSWITCH
    REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(OUT) :: qex(kdfull2, kstate)
#else
    COMPLEX(DP), INTENT(IN) :: q0(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: qex(kdfull2)
    REAL(DP), DIMENSION(:), ALLOCATABLE :: acli
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rhi
    COMPLEX(DP) :: rhoc
#endif

! workspaces

    REAL(DP), DIMENSION(:), ALLOCATABLE :: rh
    REAL(DP), DIMENSION(:), ALLOCATABLE :: acl

    LOGICAL, PARAMETER :: ttest = .false.

    IF (ifsicp /= 5) STOP ' IN EXCHANGE: wrong option IFSICP'

    ALLOCATE (acl(2*kdfull2))
    ALLOCATE (rh(2*kdfull2))

#ifdef COMPLEXSWITCH
    ALLOCATE (rhi(2*kdfull2))
    ALLOCATE (acli(2*kdfull2))
#endif

    qex = 0D0

! DOUBLE loop over states:
! outer loop for states on which exchange acts
! inner loop for states with which exchange acts

    dvol = dx*dy*dz
#ifdef REALSWITCH
!WRITE(6,'(a,10i5)') ' ispin:',ispin
!WRITE(6,'(a,10f8.2)') ' occup:',occup
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
            rh(ind) = REAL(rhoc)
            rhi(ind) = AIMAG(rhoc)
#endif
          END DO

! the Coulomb potential for the transition density
! (warning : counet inserts the esquar factor)

#if(gridfft)
          CALL falr(rh, acl, nx2, ny2, nz2, kdfull2)
#ifdef COMPLEXSWITCH
          CALL falr(rhi, acli, nx2, ny2, nz2, kdfull2)
#endif
#endif
#if(findiff|numerov)
          CALL solv_fft(rh, acl, dx, dy, dz)
#ifdef COMPLEXSWITCH
          CALL solv_fft(rhi, acli, dx, dy, dz)
#endif
#endif

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
          qex(:, nbe) = qex(:, nbe) - q0(:, nb2)*acl(:)
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

