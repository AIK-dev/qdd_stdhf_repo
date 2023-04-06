
!-----statit---------------------------------------------------------

SUBROUTINE statit(psir, rho, aloc)

! master routine for static iteration

  USE params
  USE util, ONLY: inttostring, pricm, printfield, printCubeFile, prifld, prifldz, tmodbuf
  USE coulsolv, ONLY: solv_poisson_f, solv_poisson_e, tcoulfalr
!USE coulsolv_e, ONLY:solv_poisson


  IMPLICIT NONE
  REAL(DP), INTENT(IN OUT) :: psir(kdfull2, kstate)
  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)

  LOGICAL, PARAMETER :: tcpu = .true.
  LOGICAL, PARAMETER :: tp_prints = .false.
  INTEGER :: is
  INTEGER :: i, ifsicpsav, iter1, j, nbe, nbeabs
  INTEGER :: ii, jj

  REAL(DP) :: time_init
  REAL(DP) :: xcm, ycm, zcm
!REAL(DP),ALLOCATABLE :: qaux(:,:)

  CALL calcrhor(rho, psir)

  WRITE (*, *) 'for charge=', SUM(rho)*dvol

  IF (tcoulfalr) THEN
    CALL solv_poisson_f(rho, chpcoul, kdfull2)
  ELSE
    CALL solv_poisson_e(rho, chpcoul, kdfull2)
  END IF

  CALL prifld(chpcoul, 'coulomb pot')

! number of pre-iterations for static solution with IFSICP=6.
! This parameter is to  be set here "by hand" such that it can
! be communicated by 'all.inc'

  WRITE (*, *) ' IN STATIT'

  itersicp6 = 2*istinf

  IF (tstat) THEN
    itersicp6 = 0
! CALL sstep(psir,akv,aloc,0) ! also for GSlat, DSIC
    CALL resume(psir, outnam)
    CALL calcpseudo() ! initialize pseudo-potentials
    CALL static_mfield(rho, aloc, psir, 0)
    CALL pricm(rho)
    CALL infor(rho, aloc, 0)
  END IF

! initial computation of mean field, protocol prints, headers

  CALL calcpseudo() ! initialize pseudo-potentials
  CALL static_mfield(rho, aloc, psir, 0)
  CALL pricm(rho)

  IF (myn == 0) WRITE (6, '(a,3e17.7)') 'rhoCM: ', &
    rvectmp(1), rvectmp(2), rvectmp(3)

  CALL infor(rho, aloc, 0)
  IF (myn == 0) THEN
    CALL prifld(rho, 'density ')
    CALL prifld(aloc, 'potential ')

    IF (nion2 /= 0) CALL prifld(potion, 'potential_io')
    WRITE (7, '(f8.4,a,4f12.4)') 0.0, ' initial moments', (qe(j), j=1, 4)
    WRITE (6, '(f8.4,a,4f12.4)') 0.0, ' initial moments', (qe(j), j=1, 4)
    WRITE (7, *)

    WRITE (6, '(a)') '+++ start of static iteration +++'
    WRITE (7, '(a)') '+++ start of static iteration +++'

    IF (dpolx*dpolx + dpoly*dpoly*dpolz*dpolz > 0D0) WRITE (7, '(a,3f8.4)') &
      ' static dipole potential: dpolx,dpoly,dpolz=', dpolx, dpoly, dpolz
    WRITE (7, *)
    WRITE (6, *) 'ismax=', ismax
  END IF


! the static iteration starts here

  ifsicpsav = ifsicp ! save for later recycling
  sumvar = 100D0 ! initial value for terminator

  DO iter1 = 1, ismax

    IF (tcpu) CALL cpu_time(time_init)

    IF (ifsicpsav == 4) THEN ! switch safe pre-iterations for KLI
      IF (iter1 < 40) THEN
        ifsicp = 3
      ELSE
        ifsicp = ifsicpsav
      END IF
    END IF


    CALL sstep(psir, aloc, iter1)
    CALL static_mfield(rho, aloc, psir, iter1)
    CALL pricm(rho)
    IF (myn == 0) WRITE (6, '(a,3e17.7)') 'rhoCM: ', &
      rvectmp(1), rvectmp(2), rvectmp(3)

    IF(istinf > 0) THEN; IF (MOD(iter1, istinf) == 0) THEN
      CALL infor(rho, aloc, iter1)

      IF (myn == 0) THEN
        WRITE (7, '(a,i5)') 'iter= ', iter1
        WRITE (7, '(a,f12.4,a,2(/5x,5f13.5))') 'binding energy=', binerg, &
          ', moments: monop.,dip,quad=', qe(1), qe(2), qe(3), qe(4), &
          qe(5), qe(6), qe(7), qe(8), qe(9), qe(10)
        IF (numspin == 2) WRITE (7, '(a,3f10.4)') 'spindipole', se(1), se(2), se(3)
      END IF
      IF (sumvar2 < epsoro) EXIT
    END IF; END IF

    IF (tmodbuf(iter1,isaves)) &
        CALL rsave(psir, iter1, outnam)

    FLUSH (6)

  END DO ! END iteration loop

  IF (myn == 0) THEN
    WRITE (7, *) ' static iteration terminated with ', iter1, ' iterations'

    CALL prifld(rho, 'density ')
    CALL prifld(aloc, 'potential ')
    CALL prifld(chpcoul, 'Coul-potent.')
    CALL prifldz(rho, 'density ')
    CALL prifldz(aloc, 'potential ')
    CALL prifldz(chpcoul, 'Coul-potent.')
  END IF

! compute and PRINT localization

  IF (iflocaliz) THEN
    CALL localizer(rho, psir, 0)
  END IF

! final protocol on FILE 'pstat.<NAME>''

  CALL pri_pstat(psir, rho)

! SAVE REAL wavefunctions for further applications

! CALL printSurfPot(592)

  IF (tp_prints .AND. (myn == 0 .OR. knode == 1)) THEN
    CALL printfield(491, aloc, 'tp.aloc')
    CALL printfield(492, rho, 'tp.density')
    CALL printCubeFile(493, rho)
    CALL printfield(496, chpcoul, 'tp.coulomb')
    CALL printfield(497, potion, 'tp.potion')
  END IF

  IF (tplotorbitals) THEN
    DO nbe = 1, nstate
      nbeabs = nrel2abs(nbe)
      IF (nbeabs < 1000) THEN
        OPEN (522, STATUS='unknown', FILE='pOrbitals.'//trim(adjustl(inttostring(nbeabs)))//'.'//outnam)
      ELSE
        STOP 'ERROR: Too many states for tplotorbitals'
      END IF

      WRITE (522, '(a,i3)') '# state nr: ', nbeabs
      WRITE (522, '(a,f12.5)') '# occupation: ', occup(nbe)
      WRITE (522, '(a,f12.5)') '# s.p. energy: ', amoy(nbe)
      CALL printfield(522, psir(1, nbe), 'tp.psir')
      WRITE (522, *) ! separate blocks for gnuplot
      WRITE (522, *) !
      CLOSE (522)
    END DO
  END IF

  CALL calcrhor(rho, psir)

  IF (jdensitydiff /= 0) THEN
    OPEN (590, STATUS='unknown', FILE='densdiff')
    DO i = 1, kdfull2*2
      WRITE (590, *) rho(i)
    END DO
    CLOSE (590)
  END IF

  CALL rsave(psir, iter1, outnam)
  tstat = .TRUE. ! prepare for reading in time step
    IF (itmax == 0 .AND. isaves > 0) THEN
      WRITE (*, *) ' CALL RSAVE 1. CASE'
      CALL infor(rho, aloc, -1)
      STOP ' terminate with static iteration '
    END IF

! file for groundstate at the end of static iteration

! optionally print KS potential along axes
! use 'jforce' as switch

  IF (jforce > 0 .AND. myn == 0) CALL prifld(aloc, 'potential ')

  IF (icooltyp == 1 .AND. itmax > 0) THEN
    WRITE (*, *) ' CALL RSAVE 2. CASE'
    CALL rsave(psir, iter1, outnam)
  END IF


! check 1p-h transition matrix elements to  unoccupied state:
! attention: nstate must be bigger than number of electrons
  IF (iftransme .AND. (nstate > nelect)) CALL transel(psir)



  CALL pricm(rho)
  IF (myn == 0) WRITE (6, '(a,3e17.7)') 'rhoCM: ', &
    rvectmp(1), rvectmp(2), rvectmp(3)

  xcm = rvectmp(1)
  ycm = rvectmp(2)
  zcm = rvectmp(3)


  IF (itmax <= 0) STOP ' terminate with static iteration '

  RETURN
END SUBROUTINE statit
!-----static_mfield------------------------------------------------

SUBROUTINE static_mfield(rho, aloc, psir, iter1)

! The Kohn-Sham mean field potential for given set of s.p. wavefunctions.

! Input:
! psir = REAL wavefunctions
! iter1 = nr. of static itzeration
! Output:
! rho = electron density
! aloc = local mean-field potential

  USE params


  IMPLICIT NONE

  REAL(DP), INTENT(OUT) :: rho(2*kdfull2)
  REAL(DP), INTENT(OUT) :: aloc(2*kdfull2)
  REAL(DP), INTENT(IN) :: psir(kdfull2, kstate)
  INTEGER, INTENT(IN) :: iter1

  REAL(8), PARAMETER :: addrho = 1D0

!----------------------------------------------------------------

  IF (addrho < 1D0) aloc = rho
  CALL calcrhor(rho, psir)
  IF (addrho < 1D0) rho = addrho*rho + (1D0 - addrho)*aloc

  CALL coul_mfield(rho)


  CALL calclocal(rho, aloc) ! LDA part of the potential

  IF (ifsicp > 0 .AND. ifsicp < 6) THEN
    CALL calc_sicr(rho, aloc, psir)
  END IF
  RETURN

END SUBROUTINE static_mfield

!-----sstep----------------------------------------------------------
! FFT version

SUBROUTINE sstep(q0, aloc, iter)

! Performs one static step for all wavefunctions and for given
! mean fields.
! The step involves: action of H->psi, some analysis, and damping.

! Optionally the mean field Hamiltonian is diagonalized IN the
! space of given s.p. states.

! Input/Output:
! q0 = set of s.p. wavefunctions to  be iterated
! Input:
! aloc = local Kohn-Sham field
! (akv = kinetic energy IN momentum space, via MODULE 'kinetic')
! iter = iteration number
!
! This version for FFT.

  USE params
  USE util, ONLY: wfovlp, project, tmodbuf
  USE kinetic


  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
  INTEGER, INTENT(IN) :: iter

  INTEGER :: iactsp, ii, nbc, nbcs, nbes

  INTEGER :: i, ishift
  INTEGER :: ktridig !=(kstate+kstate*kstate)/2
  INTEGER :: nbe, nph
  INTEGER :: ncount_init, ncount_rate, ncount_max, ncount_step, ncount_orth, ncount_syst, ncount_fin
  REAL(DP) :: time_init, time_step, time_orth, time_cpu, time_fin
  REAL(DP) :: sum0, sumk, sume, sum2
  REAL(DP) :: wf0, wfstep
  REAL(DP), ALLOCATABLE :: hmatr(:, :)
  REAL(DP), ALLOCATABLE :: heigen(:)
  REAL(DP), ALLOCATABLE :: vect(:, :)
  REAL(DP), ALLOCATABLE :: psistate(:)
  INTEGER, ALLOCATABLE :: npoi(:, :)
  INTEGER :: ntridig(2), nstsp(2)
  LOGICAL, PARAMETER :: tprham = .TRUE.

  REAL(DP):: espbef, espaft, sumvarold = 1D10
  INTEGER :: ni

  LOGICAL :: tocc = .TRUE. ! allow for re-occupation according to s.p. energy
  LOGICAL :: tcpu = .TRUE.


  LOGICAL, PARAMETER :: tproj = .false.
! workspaces
  REAL(DP)::vol
  REAL(DP), ALLOCATABLE :: qex(:, :)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: q1
  REAL(DP), DIMENSION(:), ALLOCATABLE :: w4
  COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: psipr
  COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: q2
  LOGICAL :: thamdiag, talert = .FALSE.
  INTERFACE
    SUBROUTINE exchgr(q0, qex, nbel)
        USE params
        USE coulsolv, ONLY: solv_poisson_f, solv_poisson_e, tcoulfalr
        REAL(DP), INTENT(IN) :: q0(kdfull2, kstate)
        REAL(DP), INTENT(OUT) :: qex(kdfull2, kstate)
        INTEGER, INTENT(IN), OPTIONAL :: nbel
    END SUBROUTINE exchgr
  END INTERFACE

!-------------------------------------------------------------------------


!  WRITE(*,*) 'reoccupations: nhstate,tocc=',nhstate,tocc

  thamdiag = tmodbuf(iter,ifhamdiag) .AND. iter > 49

  nph = 3 - numspin

! set timer

  IF (tcpu) THEN
    CALL cpu_time(time_init)
    CALL system_clock(ncount_init, ncount_rate, ncount_max)
  END IF

! exact exchange, to  be completed before wavefunctions are modified

  IF (ifsicp == 5) THEN
    ALLOCATE (qex(kdfull2, kstate))
    CALL exchgr(q0, qex)
  END IF

  ALLOCATE (q1(kdfull2))
  ALLOCATE (q2(kdfull2))


  IF (ifhamdiag > 0) THEN; IF (MOD(iter, ifhamdiag) == 0) THEN
    ktridig = (kstate + kstate*kstate)/2
    ALLOCATE (hmatr(ktridig, 2), heigen(kstate), vect(kstate, kstate))
    ALLOCATE (npoi(kstate, 2))
    ntridig(1) = 0
    ntridig(2) = 0
    nstsp(1) = 0
    nstsp(2) = 0
  END IF; END IF


!-----------------------------------------------------------------------
  sumvar = 0D0
  DO nbe = 1, nstate
    ishift = (ispin(nbe) - 1)*nxyz ! store spin=2 in upper block

! action of the potential in coordinate space
! plus non local part of ps
! plus optionally exchange part

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

    IF (ifsicp == 5) THEN
      q1 = q1 + qex(:, nbe)
    END IF

! subtract SIC potential for state NBE

! optionally compute expectation value of potential energy

    epotsp(nbe) = wfovlp(q0(:, nbe), q1) + amoy(nbe)

    ALLOCATE (psipr(kdfull2))

! action of the kinetic energy in momentum space
    CALL rftf(q0(1, nbe), psipr)

! FFT of V*psi
    CALL rftf(q1, q2)

! compose to  h|psi>
    q2 = psipr*akv + q2

! Optionally compute expectation value of kinetic energy
! and variance of mean field Hamiltonian.
! This is done in Fourier space to  SAVE FFT's.
! Variance 'evarsp2' excludes non-diagonal elements within
! occupied space.


    IF (tmodbuf(iter, istinf) .AND. ifsicp /= 6) THEN
      ALLOCATE (w4(kdfull2))
      CALL rfftback(q2, w4)
      CALL project(w4, w4, ispin(nbe), q0)
      evarsp2(nbe) = SQRT(wfovlp(w4, w4))
      DEALLOCATE (w4)
    ELSE
      evarsp2(nbe) = 0D0
    END IF

    sum0 = 0D0
    sumk = 0D0
    sume = 0D0
    sum2 = 0D0
    DO i = 1, nxyz
      vol = REAL(psipr(i), DP)*REAL(psipr(i), DP) + AIMAG(psipr(i))*AIMAG(psipr(i))
      sum0 = vol + sum0
      sumk = vol*akv(i) + sumk
      sume = REAL(q2(i), DP)*REAL(psipr(i), DP) + AIMAG(q2(i))*AIMAG(psipr(i)) &
             + sume
      sum2 = REAL(q2(i), DP)*REAL(q2(i), DP) + AIMAG(q2(i))*AIMAG(q2(i)) + sum2
    END DO

    ekinsp(nbe) = sumk/sum0
    sume = sume/sum0
    sum2 = sum2/sum0
    amoy(nbe) = sume + amoy(nbe)
    evarsp(nbe) = SQRT(MAX(sum2 - sume**2, small))
    sumvar = sumvar + occup(nbe)*evarsp(nbe)**2

    IF (ifhamdiag > 10000) THEN; IF (MOD(iter, ifhamdiag) == 0) THEN
! accumulate mean-field Hamiltonian within occupied states,
! for later diagonalization

      IF (iter > 0) THEN
        CALL rfftback(q2, q1)
        iactsp = ispin(nbe)
        nstsp(iactsp) = 1 + nstsp(iactsp)
        npoi(nstsp(iactsp), iactsp) = nbe
        DO nbc = 1, nbe
          IF (iactsp == ispin(nbc)) THEN
            ntridig(iactsp) = 1 + ntridig(iactsp)
            hmatr(ntridig(iactsp), iactsp) = wfovlp(q0(:, nbc), q1)
            IF (nbc == nbe) hmatr(ntridig(iactsp), iactsp) = &
              hmatr(ntridig(iactsp), iactsp) + amoy(nbe)
            IF (tprham) WRITE (6, '(a,2i5,1pg13.5)') ' nbe,nbc,hmatr=', nbe, nbc, &
              hmatr(ntridig(iactsp), iactsp)
          END IF
        END DO
      END IF

    END IF; END IF

! perform the damped gradient step and orthogonalize the new basis

    IF (idyniter /= 0 .AND. iter > idyniter) e0dmp = MAX(ABS(amoy(nbe)), 0.5D0)

    IF (tproj) THEN
      IF (e0dmp > small) THEN
        q2 = -q2*epswf/(akv + e0dmp)
      ELSE
        q2 = -epswf*q2
      END IF
      CALL rfftback(q2, q1)
      CALL project(q1, q1, ispin(nbe), q0)
      q0(:, nbe) = q0(:, nbe) + q1
    ELSE
      IF (e0dmp > small) THEN
        psipr = psipr - q2*epswf/(akv + e0dmp)
      ELSE
        psipr = psipr - epswf*q2
      END IF
      CALL rfftback(psipr, q0(1, nbe))
    END IF
    DEALLOCATE (psipr)
  END DO ! END loop over states
  sumvar = SQRT(sumvar/REAL(nelect, DP))

  DEALLOCATE (q1)

  IF (ifsicp == 5) DEALLOCATE (qex)

  IF (tcpu) THEN
    CALL cpu_time(time_step)
    CALL system_clock(ncount_step, ncount_rate, ncount_max)
  END IF

  IF (ifhamdiag > 10000) THEN; IF (MOD(iter, ifhamdiag) == 0) THEN

! diagonalize mean-field Hamiltonianx
! and transform occupied states correspondingly

    IF (tprham) THEN
      WRITE (6, '(a,2i5)') ' nstsp=', nstsp
      WRITE (6, '(a)') 'npoi:'
      DO iactsp = 1, 2
        WRITE (6, '(10i5)') (npoi(nbcs, iactsp), nbcs=1, nstsp(iactsp))
      END DO
    END IF

#if(!omp)
    ALLOCATE (psistate(kstate))
#endif
    DO iactsp = 1, 2
      CALL givens(hmatr(1, iactsp), heigen, vect, nstsp(iactsp), nstsp(iactsp), kstate)
      IF (tprham) WRITE (6, '(a/20(1pg13.5))') ' eigenvalues:', &
        (heigen(nbe), nbe=1, nstsp(iactsp))
      IF (tprham) WRITE (6, '(a/20(1pg13.5))') ' amoy before:', &
        (amoy(nbes), nbes=1, nstsp(iactsp))
      DO nbes = 1, nstsp(iactsp)
        nbe = npoi(nbes, iactsp)
        amoy(nbe) = heigen(nbes)
      END DO
      IF (tprham) WRITE (6, '(a/20(1pg13.5))') ' amoy after: ', &
        (amoy(nbe), nbe=1, nstsp(iactsp))
#if(omp)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii,psistate,nbes,nbe,nbcs,nbc)
      ALLOCATE (psistate(kstate))
!$OMP DO SCHEDULE(STATIC)
#endif
      DO ii = 1, nxyz
        psistate = 0D0
        DO nbes = 1, nstsp(iactsp)
          nbe = npoi(nbes, iactsp)
          amoy(nbe) = heigen(nbes)
          DO nbcs = 1, nstsp(iactsp)
            nbc = npoi(nbcs, iactsp)
            psistate(nbe) = psistate(nbe) + q0(ii, nbc)*vect(nbcs, nbes)
          END DO
        END DO
        DO nbes = 1, nstsp(iactsp)
          nbe = npoi(nbes, iactsp)
          q0(ii, nbe) = psistate(nbe)
        END DO
      END DO
#if(omp)
!$OMP END DO
      DEALLOCATE (psistate)
!$OMP END PARALLEL
#endif
    END DO
    DEALLOCATE (hmatr, heigen, vect)
    DEALLOCATE (npoi)
#if(!omp)
    DEALLOCATE (psistate)
#endif

  END IF; END IF

! Schmidt ortho-normalisation

  IF (tcpu) THEN
    CALL cpu_time(time_orth)
    CALL system_clock(ncount_orth, ncount_rate, ncount_max)
  END IF

  CALL schmidt(q0)
! invoke Hamiltonian diagonalization only if actual variance
! is factor 'variance_gain' lower than last value before diagonalization.
  IF (thamdiag .OR. talert) THEN
    IF (sumvar < sumvarold*variance_gain) THEN
      sumvarold = sumvar
      thamdiag = .TRUE.
      talert = .FALSE.
      IF (tprham) WRITE (*, *) 'Hamiltonian diag. at iter=', iter, sumvar
    ELSE
      thamdiag = .FALSE.
      talert = .TRUE.
    END IF
  END IF
  IF (thamdiag) CALL hamdiag(q0, aloc)

! readjust occupations numbers to  actual s.p. energies

  IF (tocc .AND. nelect < nstate*nph .AND. iter > 3*istinf) CALL reocc()

! protocol of SSTEP on standard output files

  IF (tcpu) THEN
    CALL cpu_time(time_fin)
    time_cpu = time_fin - time_init
    CALL system_clock(ncount_fin, ncount_rate, ncount_max)
    ncount_syst = ncount_fin - ncount_init
  END IF
  WRITE (6, '(a,i5,6(f10.4))') &
    'iter,up/down,CPU=', iter, se(4), se(5), time_cpu, ncount_syst*1D-4
  WRITE (7, '(a,i5,6(f10.4))') &
    'iter,up/down,CPU=', iter, se(4), se(5), time_cpu, ncount_syst*1D-4
  IF (tcpu) THEN
    WRITE (7, '(a,6(f9.3))') ' times: step,diag,orth=', &
      time_step - time_init, time_orth - time_step, time_fin - time_orth, &
      (ncount_step - ncount_init)*1D-4, (ncount_orth - ncount_step)*1D-4, &
      (ncount_fin - ncount_orth)*1D-4
  END IF

  RETURN
END SUBROUTINE sstep
! FFT switch finished


!-----infor--------------------------------------------------------

SUBROUTINE infor(rho, aloc, i)

! Information subroutine for static steps.
! Computes observables (energies, radii, ...) and prints to  standard output.
!
! Input:
! rho = local density
! i = iteration number

  USE params
  USE util, ONLY: printfieldx, printfieldy, printfieldz, cleanfile, emoms


  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
  INTEGER, INTENT(IN) :: i

  INTEGER :: ind, nb
  REAL(DP) :: eshell, enonlc, ensav, ekin, ehilf
  REAL(DP), PARAMETER :: alpha_ar = 10.6D0
  REAL(DP), SAVE :: energyold = 0D0
  INTEGER, PARAMETER :: idebug = 0 ! set to  1 for testprint KS potentials
  REAL(DP) :: denscrit(5) ! array for density and potential criteria


  REAL(DP) :: en(kstate)

  REAL(DP), EXTERNAL :: energ_ions

! compute s.p. energies

  eshell = 0D0
  esh1 = 0D0

  sumvar = 0D0
  sumvar2 = 0D0
  espnb = 0D0
  enonlc = 0D0
  DO nb = 1, nstate

    ekin = ekinsp(nb)
    epot = epotsp(nb)
    ehilf = epot
    epot = epot + ekin
    espnb = espnb + (epotsp(nb) + ekinsp(nb))*occup(nb)
    epot = epot + ekin
    en(nb) = epot*occup(nb)
    eshell = eshell + en(nb)
    esh1 = esh1 + ekin*occup(nb)
    sumvar = sumvar + occup(nb)*evarsp(nb)**2
    sumvar2 = sumvar2 + occup(nb)*evarsp2(nb)**2
    enonlc = enonlc + enonlo(nb)*occup(nb)
  END DO
  ecorr = energ_ions()

  DO nb = 1, nstate
    ekin = ekinsp(nb)
    IF (numspin == 2) THEN
      WRITE (6, '(a,i3,a,i3,3f9.5,2(1pg12.4))') &
        'level:', nrel2abs(nb), ' spin,occup,ekin,esp,var=', &
        3 - 2*ispin(nrel2abs(nb)), occup(nb), ekin, amoy(nb), evarsp(nb), evarsp2(nb)
    ELSE
      WRITE (6, '(a,i3,a,3f9.5,1pg12.4)') &
        'level:', nrel2abs(nb), ' occup,ekin,esp,var=', &
        occup(nb), ekin, amoy(nb), evarsp(nb)
    END IF
  END DO

  sumvar = SQRT(sumvar/REAL(nelect, DP))
  sumvar2 = SQRT(sumvar2/REAL(nelect, DP))
  IF (myn == 0) WRITE (6, '(a,2(1pg12.4))') ' total variance =', sumvar, sumvar2

  ecback = 0D0
  ecrho = 0D0
  DO ind = 1, nxyz
    IF (nion2 /= 0) THEN
      ecback = ecback - rho(ind)*potion(ind)
      ecrho = ecrho + rho(ind)*(chpcoul(ind) - potion(ind))
    ELSE IF (nion2 == 0) THEN
      ecback = ecback - rhojel(ind)*chpcoul(ind)
      ecrho = ecrho + rho(ind)*chpcoul(ind)
    END IF
  END DO
  ecback = ecback*dvol/2D0
  ecrho = ecrho*dvol/2D0

  esh1 = esh1
  eshell = eshell/2D0 !(=t+v/2)

  energy = espnb/2D0 + esh1/2D0 + enrear + ecback + ecorr + enonlc/2D0 - ecrhoimage
  IF (directenergy) &
    energ2 = esh1 + enerpw + ecrho + ecback + ecorr + enonlc - ecrhoimage
  binerg = energy
    CALL emoms(rho)
    WRITE (6, *) 'sp pot. energy =', espnb - esh1
    WRITE (6, *) 'sp kin. energy =', esh1
    WRITE (6, *) 'tot sp energy =', espnb
    WRITE (6, *) 'e_coul:ion-ion =', ecorr
    WRITE (6, *) 'e_coul:el-ion =', 2.*ecback
    WRITE (6, *) 'e_coul:el-el =', ecrho - ecback
    WRITE (6, *) 'e_coul:total =', ecback + ecrho + ecorr
    WRITE (6, *) 'rearge. energy =', enrear
    WRITE (6, *) 'nonlocal energy =', enonlc
    WRITE (6, *) 'binding energy =', energy
    IF (directenergy) THEN
      WRITE (6, *) 'binding energy2 =', energ2
      WRITE (6, *) 'potential energ =', enerpw
      ensav = energ2
      energ2 = energy
      energy = ensav
    END IF
    WRITE (6, '(a,i5,a,f12.6)') 'iter= ', i, ' binding energy', binerg
    WRITE (6, '(a,3(1pg13.5))') 'quadrupole:', qe(5), qe(6), qe(7)
    WRITE (6, '(a)') ' '

    IF (istinf > 0) THEN; IF (MOD(i, istinf) == 0) THEN
      CALL cleanfile(17)
      OPEN (17, POSITION='append', FILE='infosp.'//outnam)
      WRITE (17, '(a,i5,4(1pg13.5))') 'iteration,energy,variances=', &
        i, energy, (energy - energyold)/jinfo, sumvar, sumvar2
      FLUSH (17)
      IF (i == 0) THEN
        OPEN (199, FILE='sconver.'//outnam)
        WRITE (199, *) '# protocol of static convergence'
        WRITE (199, *) '# col 1: iteration number'
        WRITE (199, *) '# col 2: energy difference'
        WRITE (199, *) '# col 3: variance of s.p. energies'
        WRITE (199, *) '# col 4: variance of s.p. energies projected'
        WRITE (199, *) '# cols 5-9 show density and potential criteria as follows:'
        WRITE (199, *) '# col 5: sqrt of average (rho-rho_old)**2'
        WRITE (199, *) '# col 6: sqrt of average (aloc-aloc_old)**2 both spins'
        WRITE (199, *) '# col 7: radius'
        WRITE (199, *) '# col 8: cartesian quadrupole along z-axis'
        WRITE (199, *) '# col 9: rms spin polarization field'
        WRITE (199, *) '# iter diff_energ variance variance_2 density_criteria'
        CLOSE (199)
        OPEN (199, FILE='sspenergies.'//outnam)
        WRITE (199, *) '# protocol of s.p. energies'
        WRITE (199, *) '# iteration s.p. energies'
        CLOSE (199)
        OPEN (199, FILE='sspoccup.'//outnam)
        WRITE (199, *) '# protocol of s.p. occupations'
        WRITE (199, *) '# iteration s.p. occupations'
        CLOSE (199)
        OPEN (199, FILE='sspvariances.'// outnam)
        WRITE (199, *) '# protocol of s.p. variances'
        WRITE (199, *) '# iteration s.p. variances'
        CLOSE (199)
      ELSE
        OPEN (199, POSITION='append', FILE='sconver.'//outnam)
        CALL density_criterion(i, rho, aloc, denscrit)
        IF (MINVAL(ABS(denscrit)) .GE. 0D0) WRITE (199, '(i6,15(1pg13.5))') &
          i, (energy - energyold)/jinfo, sumvar, sumvar2, denscrit
        CLOSE (199)
        OPEN (199, POSITION='append', FILE='sspenergies.'//outnam)
        WRITE (199, '(i6,250(1pg13.5))') i, amoy(1:nstate)
        CLOSE (199)
        OPEN (199, POSITION='append', FILE='sspoccup.'//outnam)
        WRITE (199, '(i6,250(1pg13.5))') i, occup(1:nstate)
        CLOSE (199)
        OPEN (199, POSITION='append', FILE='sspvariances.'//outnam)
        WRITE (199, '(i6,250(1pg13.5))') i, evarsp2(1:nstate)
        CLOSE (199)
      END IF
      energyold = energy
    END IF; END IF


  IF (myn == 0) WRITE (6, '(a,f8.4,a,f7.2,a,3f9.3)') &
    'IN info: t=', tfs, ' moments: monop.=', qe(1), ' dip.: ', qe(2), qe(3), qe(4)

  IF (idebug == 1) THEN

    WRITE (6, *) 'Printing KS potential...'

    WRITE (6, *) 'Printing electron density...'
    CALL printfieldx(2000 + i, rho, 0D0, 0D0)
    CALL printfieldz(2100 + i, rho, 0D0, 0D0)

    WRITE (6, *) 'Printing ion potential...'
    CALL printfieldx(3000 + i, potion, 0D0, 0D0)
    CALL printfieldz(3100 + i, potion, 0D0, 0D0)

  END IF

! protocol of dipoles and polarizability

  IF (myn == 0 .AND. dpolx*dpolx + dpoly*dpoly + dpolz*dpolz > 0D0) THEN
    WRITE (7, '(a,3f8.4)') ' EXTERNAL dipole fields: dpolx,dpoly,dpolz=', &
      dpolx, dpoly, dpolz
    IF (dpolx > 0D0) WRITE (7, '(a,1pg13.5)') &
      ' dipole polarizability IN x=', -qe(1)*qe(2)/dpolx
    IF (dpoly > 0D0) WRITE (7, '(a,1pg13.5)') &
      ' dipole polarizability IN y=', -qe(1)*qe(3)/dpoly
    IF (dpolz > 0D0) WRITE (7, '(a,1pg13.5)') &
      ' dipole polarizability IN z=', -qe(1)*qe(4)/dpolz
    FLUSH (7)
  END IF



  RETURN
END SUBROUTINE infor

!-----pri_pstat----------------------------------------------------

SUBROUTINE pri_pstat(psi, rho)

! PRINT short protocol of final static state on FILE 'pstat.*'
!
! Input:
! psi = set of static s.p. wavefunctions
! rho = local density

  USE params
  USE coulsolv


  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: psi(kdfull2, kstate)
  REAL(DP), INTENT(IN) :: rho(kdfull2)

  LOGICAL, PARAMETER :: mumpri = .false.
  REAL(DP) :: enonlc, omegam, rms
  INTEGER :: nb
! for FALR Coulomb solver
  REAL(DP) :: p00, p10, p11r, p11i, p20, p21r, p21i, p22r, p22i, &
              p30, p31r, p31i, p32r, p32i, p33r, p33i, &
              p40, p41r, p41i, p42r, p42i, p43r, p43i, p44r, p44i, pr2


    OPEN (42, POSITION='append', FILE='pstat.'//outnam)
    IF (ifsicp == 6) THEN
      WRITE (42, '(a,i3,a,i5)') 'final protocol of static for IFSICP=', ifsicp, &
        ', pre-iterations with Slater=', itersicp6
    ELSE
      WRITE (42, '(a,i3)') 'final protocol of static for IFSICP=', ifsicp
    END IF

  IF (tcoulfalr) CALL mulmws(rho, p00, p10, p11r, p11i, p20, p21r, p21i, &
                             p22r, p22i, p30, p31r, p31i, p32r, p32i, p33r, p33i, &
                             p40, p41r, p41i, p42r, p42i, p43r, p43i, p44r, p44i, pr2, mumpri)

  DO nb = 1, nstate
    IF (numspin == 2) THEN
      WRITE (42, '(a,i3,a,i3,3f9.5,1pg12.4)') &
        'level:', nrel2abs(nb), ' spin,occup,ekin,esp,variance =', &
        3 - 2*ispin(nrel2abs(nb)), occup(nb), ekinsp(nb), amoy(nb), evarsp(nb)
    ELSE
      WRITE (42, '(a,i3,a,3f9.5,1pg12.4)') &
        'level:', nrel2abs(nb), ' occup,ekin,esp,variance=', &
        occup(nb), ekinsp(nb), amoy(nb), evarsp(nb)
    END IF
  END DO

  IF (myn == 0) THEN
    IF (temp == 0D0) THEN
      IF (directenergy) THEN
        WRITE (42, '(a,2f13.7)') 'binding energy =', energ2, energy
      ELSE
        WRITE (42, '(a,f13.7)') 'binding energy =', energy
      END IF
    ELSE
      IF (directenergy) THEN
        energy = energ2
      END IF
      WRITE (42, '(a,4f13.7)') 'energies:E,TS,F,T =', &
        energy, temp*entrop, energy - temp*entrop, temp
    END IF
    WRITE (42, '(a,1pg12.4)') 'total variance =', sumvar
    WRITE (42, '(a,4f11.5)') 'sp pot, sp kin, rearr, nonlocal=', &
      espnb - esh1, esh1, enrear, enonlc
      WRITE (42, '(a,4f11.5)') 'e_coul: i-i , e-i , e-e , total=', &
        ecorr, 2D0*ecback, ecrho - ecback, ecback + ecrho + ecorr
    WRITE (42, '(a,f7.2)') 'mon.:', qe(1)
    WRITE (42, '(a,3f11.5)') 'dip.IN :', dpolx, dpoly, dpolz
    WRITE (42, '(a,3f11.5)') 'dip.OUT :', qe(2), qe(3), qe(4)
    WRITE (42, '(a)') 'quadrupole moments:'
    WRITE (42, '(a,3f11.4)') 'xx,yy,zz:', qe(5), qe(6), qe(7)
    WRITE (42, '(a,3f11.4)') 'xy,zx,zy:', qe(8), qe(9), qe(10)
    rms = SQRT(qe(5) + qe(6) + qe(7))
    IF (tcoulfalr) THEN
      WRITE (42, '(a,3f11.4)') 'q20,30,4:', p20, p30, p40
      WRITE (42, '(a,3f11.4)') 'renorm. :', &
        p20/(qe(1)*rms**2), p30/(qe(1)*rms**3), p40/(qe(1)*rms**4)
    END IF
    WRITE (42, '(a,3f11.4)') 'spindip.:', se(1), se(2), se(3)
    WRITE (42, '(a,f11.4,a,2(1pg13.5))') 'rms radius:', rms, &
      ', average density,k_F:', &
      qe(1)*3D0/(4D0*PI*SQRT(1.6666667D0)*rms**3), &
      (2.25D0*PI)**0.33333D0/rms
  END IF

  IF (ifspemoms) CALL spmomsr(psi, 42)

  IF (myn == 0) THEN
    WRITE (42, '(1x)')
    CLOSE (42)
  END IF

  RETURN
END SUBROUTINE pri_pstat

!-----transel-------------------------------------------------transel--

SUBROUTINE transel(psir)

! Dipole transition matrix elements between all actual s.p. states.
! result is printed on file 'mte_xyz'.
!
! Input:
! psir = set of stationary s.p. wavefunctions

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: psir(kdfull2, kstate)

  INTEGER :: i, ind, ix, iy, iz, nb1, nb2
  REAL(DP) :: x1, y1, z1
  REAL(DP) :: te(kstate, kstate, 3)

  OPEN (34, STATUS='unknown', FILE='mte-xyz')
  WRITE (34, *) '# spin n1 n2 dipol(1:3) e_ph'
  DO nb1 = 1, nstate
    DO nb2 = 1, nstate
      te(nb1, nb2, 1) = 0D0
      te(nb1, nb2, 2) = 0D0
      te(nb1, nb2, 3) = 0D0

      IF (ispin(nb1) == ispin(nb2) .AND. occup(nb2) < 0.1D0 &
          .AND. occup(nb1) > 0.9D0) THEN
        ind = 0
        DO iz = minz, maxz
          z1 = (iz - nzsh)*dz
          DO iy = miny, maxy
            y1 = (iy - nysh)*dy
            DO ix = minx, maxx
              ind = ind + 1
              IF ((ix /= nx2) .AND. (iy /= ny2) .AND. (iz /= nz2)) THEN
                x1 = (ix - nxsh)*dx

                te(nb1, nb2, 1) = te(nb1, nb2, 1) + psir(ind, nb1)*x1*psir(ind, nb2)
                te(nb1, nb2, 2) = te(nb1, nb2, 2) + psir(ind, nb1)*y1*psir(ind, nb2)
                te(nb1, nb2, 3) = te(nb1, nb2, 3) + psir(ind, nb1)*z1*psir(ind, nb2)

              END IF
            END DO
          END DO
        END DO
      END IF
      te(nb1, nb2, 1) = ABS(te(nb1, nb2, 1)*dvol)**2D0
      te(nb1, nb2, 2) = ABS(te(nb1, nb2, 2)*dvol)**2D0
      te(nb1, nb2, 3) = ABS(te(nb1, nb2, 3)*dvol)**2D0

      IF (te(nb1, nb2, 1) > 0D0 .OR. te(nb1, nb2, 2) > 0D0 &
          .OR. te(nb1, nb2, 3) > 0D0) THEN
        WRITE (34, '(3i3,3(1pg13.5),0pf12.5)') ispin(nb1), nb1, nb2, &
          (te(nb1, nb2, i), i=1, 3), amoy(nb2) - amoy(nb1)
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE transel

! ******************************

SUBROUTINE reocc()

! ******************************

! Readjust occupations numbers to  actual s.p. energies.
! I/O through module 'params'.
!

  USE params
  USE util, ONLY: pair, sort_energy
  IMPLICIT NONE

  INTEGER :: is, nbe, numpart
  INTEGER :: nelecs(2), nstspin(2, 2)
  REAL(DP) :: epstmp, occo, partnm, gp
  REAL(DP) :: occold(kstate), ocwork(kstate)
  REAL(DP) :: ph(ksttot) ! degeneracy of wavefunction, for
  REAL(DP), ALLOCATABLE :: eord(:)
  INTEGER, ALLOCATABLE :: isort(:)
  REAL(DP), SAVE :: efermsav(2) = (/0D0, 0D0/)
  LOGICAL, PARAMETER :: tprintp = .FALSE.
  LOGICAL :: tdown

!--------------------------------------------------------------------

  DO nbe = 1, nstate
    occold(nbe) = occup(nbe)
    ocwork(nbe) = occold(nbe)
  END DO
  epstmp = 1D-8
  occo = 1D0 - occmix

  IF (temp == 0D0) THEN

    ALLOCATE (eord(nstate), isort(nstate))
    eord(1:nstate) = amoy(1:nstate)
    occup(1:nstate) = 0D0
    CALL sort_energy(eord, isort, nstate)
    IF (numspin == 2) THEN
      nelecs(1) = nelect - nspdw
      nelecs(2) = nspdw
    ELSE
      nelecs(1) = nelect
    END IF

    DO is = 1, numspin
      numpart = 0
      DO nbe = 1, nstate
        IF (ispin(isort(nbe)) == is) THEN
          occup(isort(nbe)) = 1D0
          numpart = 1 + numpart
          IF (numpart == nelecs(is)) EXIT
        END IF
      END DO
    END DO

    DEALLOCATE (eord, isort)

  ELSE IF (numspin == 2) THEN
    nelecs(1) = nelect - nspdw
    nelecs(2) = nspdw
    ph = 1D0
    nstspin(1, 1) = 1
    nstspin(2, 1) = 0
    tdown = .FALSE.
    DO nbe = 1, nstate
      IF (ispin(nbe) .NE. 1) THEN
        tdown = .TRUE.
      ELSE
        IF (.NOT. tdown) nstspin(2, 1) = nbe
        IF (nbe > nstspin(2, 1)) STOP "non-continguous spin"
      END IF
      nstspin(1, 2) = nstspin(2, 1) + 1
      nstspin(2, 2) = nstate
    END DO

    DO is = 1, numspin
      eferm = efermsav(is)
      IF (tprintp) WRITE (*, *) 'REOCC: is,etc=', &
        is, nstspin(1:2, is), nstspin(2, is) - nstspin(1, is) + 1, temp
      CALL pair(amoy(nstspin(1, is):nstspin(2, is)), &
                ocwork(nstspin(1, is):nstspin(2, is)), &
                ph(nstspin(1, is):nstspin(2, is)), &
                nelecs(is), nstspin(2, is) - nstspin(1, is) + 1, &
                gp, eferm, temp, partnm, 200, 4, epstmp, -1, ksttot)
      efermsav(is) = eferm
      IF (tprintp) THEN
        WRITE (6, '(a,2i5,4(1pg13.5))') 'REOCC basic:', nelecs(is), nstate, &
          eferm, temp, partnm, occmix
        WRITE (6, '(a//20(5(1pg13.5)/))') 'REOCC energ:', amoy(nstspin(1, is):nstspin(2, is))
        WRITE (6, '(a//20(5(1pg13.5)/))') 'REOCC phase:', ph(nstspin(1, is):nstspin(2, is))
        WRITE (6, '(a//20(5(1pg13.5)/))') 'REOCC occup:', ocwork(nstspin(1, is):nstspin(2, is))
      END IF
      DO nbe = 1, nstate
        IF (ispin(nbe) == is) THEN
          occup(nbe) = occmix*ocwork(nbe)*ph(nbe) + occo*occold(nbe)
        END IF
      END DO
    END DO
  ELSE
    eferm = 0D0 ! restart search of eferm from scratch
    IF (tprintp) THEN
      WRITE (6, '(a,2i5,3(1pg13.5))') 'REOCC basic:', nelecs(1), nstate, &
        eferm, temp, partnm
      WRITE (6, '(a,200(1pg13.5))') 'REOCC energ:', amoy(1:nstate)
    END IF
    ph = 2D0
    CALL pair(amoy, ocwork, ph, nelect, nstate, gp, eferm, temp, partnm, &
              60, 4, epstmp, -1, ksttot)
    IF (tprintp) WRITE (6, '(a,200(1pg13.5))') 'REOCC occup:', ocwork(1:nstate)
    DO nbe = 1, nstate
      occup(nbe) = occmix*ocwork(nbe)*ph(nbe) + occo*occold(nbe)
    END DO
  END IF
  RETURN
END SUBROUTINE reocc

!-----rhpsi -------------------------------------------------------------

SUBROUTINE rhpsi(qact, aloc, nbe, itpri)

! action of mean-field Hamiltonian on A REAL s.p. wavefunction:
! qact = wavefunction on which H acts and resulting w.f.
! aloc = local potential for the actual spin component
! ak = kinetic energies in momentum space
! nbe = number of state
! itpri = switch for computing s.p. energies (for ABS(itpri)=1)
! <0 switches to  subtract mean-value of s.p. energy
!

  USE params
  USE util, ONLY: wfovlp
  USE kinetic

  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: qact(kdfull2)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN) :: akv(kdfull2)
  INTEGER, INTENT(IN) :: nbe
  INTEGER, INTENT(IN) :: itpri
  REAL(DP), ALLOCATABLE :: qex(:)

! workspaces
  REAL(DP), ALLOCATABLE :: q1fine(:), q2fine(:), qactfine(:), q1(:), q2(:)
  COMPLEX(DP), ALLOCATABLE :: qc(:)
  REAL(DP), ALLOCATABLE :: qarray(:, :, :), qarrayfine(:, :, :)
  REAL(DP) :: wnorm
  LOGICAL :: tpri
  LOGICAL, PARAMETER :: tsubmean = .TRUE.
  LOGICAL, PARAMETER :: ttest = .FALSE.
  INTEGER :: i, is, na
  REAL(DP) :: cf

!----------------------------------------------------------------------

  tpri = ABS(itpri) == 1

  ALLOCATE (q1(kdfull2), q2(kdfull2))
  q1 = 0D0
  q2 = 0D0
! ACTION of kinetic energy

  ALLOCATE (qc(kdfull2))
  CALL rftf(qact, qc)
  DO i = 1, nxyz
    qc(i) = akv(i)*qc(i)
  END DO
  CALL rfftback(qc, q2)
  DEALLOCATE (qc)

! action of potential and non-local PsP (optionally)

  IF (ipsptyp == 1) THEN
    CALL nonlocalr(qact, q1)

    IF (tpri) enonlo(nbe) = wfovlp(qact, q1)
    DO i = 1, nxyz
      q1(i) = q1(i) + qact(i)*aloc(i)
    END DO
  ELSE
    DO i = 1, nxyz
      q1(i) = qact(i)*aloc(i)
    END DO
  END IF

  IF (ifsicp == 5) THEN
    ALLOCATE (qex(kdfull2))
    CALL exchgr(qact, qex, nbe)
    q1 = q1 + qex
    DEALLOCATE (qex)
  END IF

! subtract SIC potential for state NBE

  IF (tpri) THEN
    q2 = q1 + q2
    spvariance(nbe) = SQRT(MAX(REAL(wfovlp(q2, q2), DP) - ABS(wfovlp(qact, q2))**2, 1D-99))
    is = ispin(nrel2abs(nbe))
    IF (ttest) WRITE (*, '(a,2i4,5(1pg13.5))') &
      ' HPSI: nbe,is,esp,var=', nbe, is, amoy(nbe), spvariance(nbe), &
      ekinsp(nbe), epotsp(nbe), REAL(wfovlp(q2, q2), DP)
    FLUSH (6)
  ELSE
    q2 = q1 + q2
  END IF

  qact = q2

  DEALLOCATE (q1, q2)

  RETURN
END SUBROUTINE rhpsi

!-----density_criterion---------------------------------------------------------

SUBROUTINE density_criterion(iter, rho, aloc, denscrit)

! change of density and potential as convergence criterion
!
! Input:
! iter = iteration number at which routine is called
! rho = density
! aloc = local potentials for spin-up and spin-down
!
! Output:
! denscrit = array of criteria
! 1 = rms change of density 'rho' (total density)
! 2 = rms change of potential'aloc' (averaged spin up and down)
! 3 = change of rms radius
! 4 = change of quadrupole about z-axis
! 5 = rms change of spin polarization field

  USE params

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iter
  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
  REAL(DP), INTENT(OUT) :: denscrit(5)

  REAL(DP) :: diffnorm, diffiter
  INTEGER, SAVE :: iter_old, nmin, nmax
  LOGICAL :: tfirst = .TRUE.
  REAL(DP), ALLOCATABLE, SAVE :: rho_old(:), aloc_old(:), qe_old(:)

  IF (tfirst) THEN
    ALLOCATE (rho_old(2*kdfull2))
    ALLOCATE (aloc_old(2*kdfull2))
    ALLOCATE (qe_old(5:7))
    rho_old = rho
    aloc_old = aloc
    iter_old = iter
    qe_old(5:7) = qe(5:7)
    denscrit = -1D0 ! to  indicate undefined staus
    tfirst = .FALSE.
    RETURN
  END IF

  diffiter = 1D0/(iter - iter_old)
  diffnorm = diffiter/kdfull2

! total density
  denscrit(1) = SQRT(SUM((rho(1:kdfull2) - rho_old(1:kdfull2))**2)*diffnorm)

! local potential
  nmax = 2*kdfull2
  denscrit(2) = SQRT(SUM((aloc(1:nmax) - aloc_old(1:nmax))**2)*diffnorm/2D0)

! rms radius
  denscrit(3) = (SQRT(SUM(qe(5:7))) - SQRT(SUM(qe_old(5:7))))*diffiter

! cartesiian quadrupole along z axis
  denscrit(4) = ((2D0*qe(7) - qe(6) - qe(5)) - (2D0*qe_old(7) - qe_old(6) - qe_old(5))) &
                *diffiter

! rms spin polarization
  nmin = kdfull2 + 1
  nmax = 2*kdfull2
  denscrit(5) = SQRT(SUM((rho(nmin:nmax) - rho_old(nmin:nmax))**2)*diffnorm)

  rho_old = rho
  aloc_old = aloc
  iter_old = iter
  qe_old(5:7) = qe(5:7)

  RETURN

END SUBROUTINE density_criterion
