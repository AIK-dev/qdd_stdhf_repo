!-----dyn_propag---------------------------------------------------------

SUBROUTINE dyn_propag(psi, rho, aloc)

!
! Master routine for dynamic propagation of electrons and ions.
!
! Input/Ouput:
! psi = set of s.p. wavefunctions to be propagated
! rho = local density
! aloc = local KS potential
! Other variables (e.g. ionic coordinates) are communicated
! via MODULE 'params'.

  USE params
  USE kinetic
  USE util, ONLY: timer, safeopen, testcurrent, rhointxy, rhointyz, rhointxz, tmodbuf
  USE RTA_module

#if(fsic)
  USE twost, ONLY: tnearest
#endif

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
  REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)

  COMPLEX(DP), ALLOCATABLE :: psiw(:, :)

  INTEGER :: it, ion
  REAL(DP):: totalprob, totalovlp
  REAL(DP), EXTERNAL:: energ_ions ! declared inion_md
  REAL(DP), EXTERNAL:: enerkin_ions ! declared inion_md

  LOGICAL, PARAMETER :: ttimestop = .FALSE.
  REAL(DP) :: t1, t2, acctime(2) = (/0D0, 0D0/)
  REAL(DP) :: spenergbef(kstate),spenergaft(kstate),Eoutloc!mv
  INTEGER :: it1, it2 , isp
  complex(dp)::Eoutstd(kstate),Eoutmat(kstate,kstate)!mv

  INTERFACE
    SUBROUTINE tstep(q0, aloc, rho, it)
      USE params
      COMPLEX(DP), INTENT(IN OUT), TARGET :: q0(kdfull2, kstate)
      REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
      REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
      INTEGER, INTENT(IN) :: it
    END SUBROUTINE
  END INTERFACE

!--- here starts true propagation --------------
  Eoutstd= cmplx(0.d0,0.d0,8)!mv
  Eoutmat= cmplx(0.d0,0.d0,8)!mv
  Eoutloc=0;
  energref=energy

  CALL dyn_mfield(rho, aloc, psi, 0D0, 0)
  IF (ifexpevol) ALLOCATE (psiw(kdfull2, kstate))
  FLUSH (7)

  WRITE (*, *) 'before loop: cpx,y,z:', cpx(1:nion), cpy(1:nion), cpz(1:nion)
  CALL testcurrent(psi, rho, 0)

  time = 0
  tfs = 0
  DO it = irest, itmax ! time-loop

    tfs = it*dt1*0.0484D0 !/(2.0*ame)

    IF (ttimestop) CALL system_clock(it1)

    CALL print_densdiff(rho, it) 
    IF (ttimestop) THEN
      CALL system_clock(it2)
      WRITE (6, '(a,2(1pg13.5))') &
        'TIMES print_densdiff:', (it2 - it1)*systime_factor
      acctime(2) = acctime(2) + (it2 - it1)*systime_factor
    END IF

    IF (jrtaint .NE. 0) THEN
!  CALL testcurrent(psi, rho, it)
      IF (MOD(it, jrtaint) .EQ. 0 .AND. it > 0 .AND. (it > irest)) THEN 
        CALL rta(psi, aloc, rho, it)
      END IF
    END IF

#if(extended)
    IF (jattach /= 0 .AND. it == irest) THEN
      CALL init_occ_target()
      WRITE (*, *) 'nstate_target, after init_occ_target:', nstate_target
      ALLOCATE (psi_target(kdfull2, nstate_target))
      CALL init_psitarget()
    END IF
#endif

    IF (it > irest) THEN

! pure electronic dynamics

      IF (nelect > 0) THEN

! propagation of the wfs
        IF (ttimestop) CALL system_clock(it1)

        WRITE (*, *) 'propagation of the wfs'
        IF (ifexpevol) THEN
          CALL tstep_exp(psi, aloc, rho, it, psiw, .FALSE.) ! expon. evolution
        ELSE
#if(extended)
! IF(ifcnevol == 1) THEN
! WRITE(*,*) 'Crank CALL '
! CALL CrankNicolson_exp(psi,aloc,rho,it) ! yet experimental
! ELSE
          CALL tstep(psi, aloc, rho, it) ! T-V splitting
! END IF
#else
          CALL tstep(psi, aloc, rho, it) ! T-V splitting
#endif
        END IF

        IF (jlaser > 0) THEN;  IF (MOD(it, jlaser) == 0) THEN
          WRITE (7, '(a,2f9.3,3(1pg12.4))') ' tfs,tRy,foft1,foft2,foftXUV=', &
            tfs, dataold(:4)
          WRITE (38, '(1f15.5,14(1pg12.4))') tfs, datalaser(:14)
        END IF; END IF
       
        IF (nabsorb > 0) then 
          if(jstochasticint>0)CALL get_spenerg(psi,aloc,spenergbef) !!JH!mv
          CALL absbc(psi, rho)
          if(jstochasticint>0)CALL get_spenerg(psi,aloc,spenergaft) !!JH!mv
        endif
        if(jstochasticint>0) then
        
          Eoutstd= Eoutstd + spenergbef - spenergaft!!JH tableau
          Eoutloc=0.d0
          do isp=1,nstate; Eoutmat(isp,isp)=Eoutstd(isp);Eoutloc=Eoutloc+Eoutmat(isp,isp)*occup(isp); enddo           
          write(165,'(f8.5,4e14.6)') tfs,dot_product(real(Eoutstd(1:nstate)),occup(1:nstate)),energy,Eoutloc,energy+eoutloc
          call flush(165)
          if (mod(it,jstochasticint).eq.0.and.it>0.and.(it>irest)) then
              write(*,*) 'call stochastic'
              call stochastic(psi,rho,aloc,eoutmat,it)!JH
          endif
        endif
        IF (ttimestop) THEN
          CALL system_clock(it2)
          WRITE (6, '(a,2(1pg13.5))') &
            'TIMES electronic step:', (it2 - it1)*systime_factor
          acctime(2) = acctime(2) + (it2 - it1)*systime_factor
        END IF

! protocol of densities

        IF (jdensity1d > 0) THEN; IF (MOD(it, jdensity1d) == 0) THEN
          CALL rhointxy(rho, it)
          CALL rhointxz(rho, it)
          CALL rhointyz(rho, it)
        END IF; END IF
        IF (tmodbuf(it,jstinf)) CALL testcurrent(psi, rho, it)

#if(raregas)
      ELSE
        IF (isurf /= 0 .AND. NE > 0) CALL valence_step(rho, dt1, .true.)
#endif
      END IF

      IF (ionmdtyp == 1 .OR. &
          (ionmdtyp == 2 .AND. tmodbuf(it,modionstep))) THEN
        IF (ttimestop) CALL system_clock(it1)

! ionic propagation

        IF (ionmdtyp == 1) CALL itstep(rho, it, psi)
        IF (ionmdtyp == 2) CALL itstepv(rho, it, psi)
        ekion = enerkin_ions()
        IF (icooltyp > 0) CALL reset_ions()

! new total mean field 
#if(extended)
!('calclocal' still necessary ?)
#endif

        IF (nelect > 0) THEN
          CALL calcpseudo()
          CALL calclocal(rho, aloc) 
          IF (ifsicp > 0) CALL calc_sic(rho, aloc, psi)
          IF (ipsptyp == 1) THEN
          DO ion = 1, nion
#if(extended)
          IF (tinterpol) THEN
            CALL calc_projFine(cx(ion), cy(ion), cz(ion), &
                               cx(ion), cy(ion), cz(ion), ion)
          ELSE
            CALL calc_proj(cx(ion), cy(ion), cz(ion), &
                           cx(ion), cy(ion), cz(ion), ion)
          ENDIF
#else
          CALL calc_proj(cx(ion), cy(ion), cz(ion), &
                         cx(ion), cy(ion), cz(ion), ion)
#endif
          END DO
          END IF
        END IF
        IF (ttimestop) THEN
          CALL system_clock(it2)
          WRITE (6, '(a,2(1pg13.5))') &
            'TIMES ionic step:', (it2 - it1)*systime_factor
          acctime(2) = acctime(2) + (it2 - it1)*systime_factor
        END IF
      END IF

    END IF

! ******** compute and WRITE observables: ********

    CALL timer(2)

    IF (it > irest) THEN
      IF (ttimestop) CALL system_clock(it1)
      IF (myn == 0) THEN
        tfs = it*dt1*0.0484D0
        IF (nion2 > 0) CALL analyze_ions(it)
#if(raregas)
        IF (isurf > 0) CALL analyze_surf(it)
#endif
      END IF
      IF (nelect > 0) THEN
        CALL analyze_elect(psi, rho, aloc, it)
      ELSE IF (MOD(it, 100) == 0) THEN
        ecorr = energ_ions()
        etot = ecorr + ekion
      END IF
      IF (ttimestop) THEN
        CALL system_clock(it2)
        WRITE (6, '(a,2(1pg13.5))') 'TIMES analysis:', (it2 - it1)*systime_factor
        acctime(2) = acctime(2) + (it2 - it1)*systime_factor
      END IF

#if(extended)
      ! The calculation of an electron attachment on a water molecule should
      ! proceed as follows:
      ! 1) Perform a static calculation for a water molecule alone, with a
      ! 'deocc' high enough, so that all bound states, occupied and empty,
      ! are calculated.
      ! 2) Create the FILE 'occ_spe_target' containing the binding energy and
      ! the total energy as a first line, and THEN the isopin, occupation
      ! NUMBER and s.p. energie of each state.
      ! 3) Perform a dynamical calculation with itmax=1 to generate a SAVE FILE
      ! which should be moved to the parent directory and renamed 'wfs_target'.
      ! 4) Perform a static and dynamical calculation with 'iscatterelectron=1'
      ! and a small deocc (e.g., 'deocc=0.1'), so that ONLY occupied states
      ! of the water molecule are considered. The incoming electron is added
      ! just before the dynamics as the last (occupied) state, and thus
      ! labeled by nstate.
      IF (jattach /= 0 .AND. it > irest) THEN
      IF (MOD(it, jattach) == 0) THEN
        CALL attach_prob(totalprob, totalovlp, psi)
        totintegprob = totintegprob + dt1*0.0484D0*jattach*totalprob
        WRITE (6, '(a,e12.5,1x,i8,3(1x,1pg13.5))') &
          'after ATTACHEMENT:', &
          tfs, nmatch, totalprob, totintegprob, totalovlp
        CALL safeopen(809, it, jattach, 'pattach')
        WRITE (809, '(e12.5,1x,i8,3(1x,1pg13.5))') &
          tfs, nmatch, totalprob, totintegprob, totalovlp
        FLUSH (809)
      END IF
      END IF
#endif

      IF (isaved > 0 .AND. it /= 0) THEN
      IF (MOD(it, isaved) == 0 .OR. it == itmax) THEN
      IF (irest /= 0 .AND. ABS(it - irest) <= 2) THEN
        ! DO nothing
      ELSE
        CALL SAVE(psi, it, outnam)
      END IF
      END IF
      END IF
#if(extended)
      IF (trequest > 0D0) CALL savings(psi, it)
#endif
    END IF

  END DO

! ******************** END of dynamic loop ****************************

  IF (ifexpevol) DEALLOCATE (psiw)

  IF (ttimestop) WRITE (6, '(a,2(1pg13.5))') 'TIMES ACCUMULATED:', acctime

  RETURN

END SUBROUTINE dyn_propag

!------------------------------------------------------------

SUBROUTINE init_dynwf(psi)

! Initializates for dynamic evolution.
!
! Output:
! psi = set of s.p. wavefunctions
! Other variables though module 'params'

  USE params
  USE util, ONLY: stateoverl, dipole_qp, prifld, testcurrent
#if(fsic)
  USE twost
#endif
  USE orthmat

  IMPLICIT NONE

! initializes dynamical wavefunctions from static solution

  COMPLEX(DP), INTENT(OUT) :: psi(kdfull2, kstate)

  INTEGER :: ifreset
  REAL(DP) ::acc, palph, pbeta, phexe, sqtest, xion
  REAL(DP), ALLOCATABLE :: rhosav(:), alocsav(:), rho(:), aloc(:)
  LOGICAL, PARAMETER :: ttestrho = .FALSE.
!----------------------------------------------------

#if(fsic)
  IF (ifsicp .GE. 8) CALL expdabvol_rotate_init ! MV initialise ExpDabOld
#endif

  IF (ifsicp .EQ. 5 .AND. .NOT. ifexpevol) &
    STOP ' exact exchange requires exponential evolution'
  IF (ifsicp .EQ. 5 .OR. jstateoverlap > 0) ALLOCATE (psisavex(kdfull2, kstate))

! USE of saved REAL wavefunctions, READ and copy to COMPLEX

  CALL restart2(psi, outnam, .true.)
#if(extended)
  IF (nhstate > 0 .AND. npstate == 0) CALL hstate()
#endif

  WRITE (*, *) 'ttestrho=', ttestrho
  IF (ttestrho) THEN
    ALLOCATE (rhosav(2*kdfull2), alocsav(2*kdfull2))
    ALLOCATE (rho(2*kdfull2), aloc(2*kdfull2))
    rhosav = rho
    alocsav = aloc
    CALL dyn_mfield(rho, aloc, psi, 0D0, 0)
    WRITE (*, *) 'INIT-DYN diff rho,aloc:', SUM(ABS(rho - rhosav)), &
      SUM(ABS(aloc - alocsav))
    CALL testcurrent(psi, rho, 0)
    DEALLOCATE (rhosav, alocsav)
    DEALLOCATE (rho, aloc)
  END IF

  IF (ifsicp .EQ. 5 .OR. jstateoverlap > 0) psisavex = psi

#if(raregas)
! boost the wavefunctions in case of ionic c.m. boost
  IF (ekin0pp > 0D0) CALL init_velocel(psi)
#endif

#if(extended)
! optionally rotate a 1ph state, only for finite phase angle
  IF (ABS(phangle) + ABS(phphase) > small) CALL phstate(psi)
! optionally a 2. 1ph rotation
  IF (ABS(phangle2) + ABS(phphase2) > small) THEN
    nhstate = nhstate2
    npstate = npstate2
    phangle = phangle2
    phphase = phphase2
    CALL phstate(psi)
  END IF

! atomic projectile wavefunctions
  IF (eproj > 0D0) CALL init_projwf(psi)

! add scattering electron

  IF (iscatterelectron /= 0) CALL init_scattel(psi)

#endif


! OPTIONAL initial excitations, either boost or rotation

  IF (irotat == 0) THEN
    IF (ABS(shiftinix) + ABS(shiftiniy) + ABS(shiftiniz) > small) THEN
      CALL dipole_qp(psi, shiftinix, shiftiniy, shiftiniz, centfx, centfy, centfz)
    ELSE
      CALL boost(psi) !dipole boost of electr. cloud
    END IF
  ELSE
    IF (nion2 /= 0) THEN
! rotate ionic coordinates
      CALL mrote
      CALL calcpseudo()
    ELSE IF (nion2 == 0) THEN
! rotate jelliumdensity:
      palph = falph
      pbeta = fbeta
      phexe = fhexe
      sqtest = 0D0
      xion = 1D0*nion
      CALL jelbak(xion, palph, pbeta, phexe, sqtest, 1)
      acc = dvol*SUM(rhojel)
      WRITE (6, *) 'jellium-density after rotation:', acc
      WRITE (7, *) 'jellium-density after rotation:', acc
      WRITE (7, *) ' '
    END IF
  END IF

  IF (jstateoverlap > 0) THEN
    CALL stateoverl(psi, psisavex)
! DEALLOCATE(psisavex)
  END IF

  RETURN
END SUBROUTINE init_dynwf

!-----init_velocel-------------------------------------------------

SUBROUTINE init_velocel(psi)

! Boosts all electronic wavefunctions 'psi' by the same velocity.
! The absolute value is associated with the boost kinetic
! energy 'ekin0pp' and the direction is given by 'vxn0', 'vyn0','vzn0'.
! This routine is similar to 'boost', but it is used in
! connection with ionic initialization to produce a cluster
! where electrons and ions move with the same velocity.

  USE params
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)

  INTEGER :: ind, ix, iy, iz, nbe
  COMPLEX(DP) :: cfac
  REAL(DP) :: arg, rnorm, v0, x1, y1, z1
!------------------------------------------------------------------

  v0 = SQRT(2D0*ekin0pp/(amu(np(nion))*1836.0D0*ame))
  rnorm = vxn0**2 + vyn0**2 + vzn0**2
  rnorm = SQRT(rnorm)

  IF (rnorm == 0) STOP 'Velocity vector not normalizable'

  vxn0 = vxn0/rnorm*v0*ame
  vyn0 = vyn0/rnorm*v0*ame
  vzn0 = vzn0/rnorm*v0*ame

  ind = 0
  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        ind = ind + 1
        arg = x1*vxn0 + y1*vyn0 + z1*vzn0
        cfac = CMPLX(COS(arg), SIN(arg), DP)
        DO nbe = 1, nstate
          psi(ind, nbe) = cfac*psi(ind, nbe)
        END DO
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE init_velocel

!-----tstep---------------------------------------------------------

SUBROUTINE tstep(q0, aloc, rho, it)

! One electronic time step by TV splitting method.
!
! Input:
! it = time step incalling routine
! Input/Output:
! q0 = set of s.p. wavefunctions
! rho = local density (spin-up and spin-down)
! aloc = local KS potential

  USE params
  USE kinetic
  USE util, ONLY: stimer, tmodbuf
#if(fsic)
  USE twost, ONLY: tnearest, eval_unitrot
!USE twost_util, ONLY:
#endif

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT), TARGET :: q0(kdfull2, kstate)
  REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
  REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
  INTEGER, INTENT(IN) :: it

  COMPLEX(DP), DIMENSION(:, :), ALLOCATABLE :: q1, q2
  COMPLEX(DP), DIMENSION(:, :), ALLOCATABLE :: qwork
  COMPLEX(DP), DIMENSION(:, :, :), POINTER :: q3D

  LOGICAL :: tenerg
  INTEGER :: ind, ishift, ithr, itsub, nb, nlocact, nup, ntm, nta, ita
  INTEGER :: ncount_init, ncount_rate, ncount_max, ncount_syst, ncount_fin
  REAL(DP) :: ri, dt, pr, time_init, time_fin, time_cpu

  LOGICAL, PARAMETER :: ttimestop = .FALSE.
  REAL(DP) :: t1, t2
! acctime is accumulator for time usage
! acctime(1) = preparation of local propagator
! acctime(2) = local propagator
! acctime(3) = nonlocal propagator
! acctime(4) = kinetic propagator
! acctime(5) = mean field
! acctime(6) = TSTEP
! acctime(7) = nr. of calls to TSTEP
  REAL(DP), SAVE :: acctime(1:10) = &
                    (/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/)

  INTEGER :: it1, it2, it3, it4, it5

#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)

  CALL mpi_comm_rank(mpi_comm_world, myn, mpi_ierror)
#else
  myn = 0
#endif

  acctime(1:5) = 0D0
  ALLOCATE (q1(2*kdfull2, 0:nthr))
  ALLOCATE (q2(kdfull2, 0:nthr))

  CALL cpu_time(time_init)
  CALL system_clock(ncount_init, ncount_rate, ncount_max)
#if(extended)
  IF (ntref > 0 .AND. it > ntref) nabsorb = 0 ! is that the correct place?
#endif

! 'itsub' indicates the number of subiteration before
! next analyzing step (routine 'info').
! itsub = MOD(it, ipasinf) + 1

  ri = -dt1*0.5D0
  dt = dt1*0.5D0
  nlocact = numspin*nxyz

! half time step incoordinate space
! local phase field on workspace 'q1'

  IF (ttimestop) THEN
    acctime(1:5) = 0D0
    CALL system_clock(it1)
  END IF

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ind) SCHEDULE(STATIC)
  DO ind = 1, nlocact
#if(extended)
! pr=-dt*aloc(ind)
! q1(ind,0)=CMPLX(COS(pr),SIN(pr),DP)
#endif
    q1(ind, 0) = EXP(CMPLX(0D0, -dt*aloc(ind), DP))
  END DO
!$OMP END PARALLEL DO

  IF (ttimestop) CALL system_clock(it2)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nb,ishift) SCHEDULE(STATIC)
  DO nb = 1, nstate
    ishift = (ispin(nrel2abs(nb)) - 1)*nxyz
#if(extended)
! CALL cmult3d(q0(1,nb),q1(1+ishift))
#endif
    q0(:, nb) = q1(ishift + 1:ishift + kdfull2, 0)*q0(:, nb)
  END DO
!$OMP END PARALLEL DO

  IF (ttimestop) THEN
    CALL system_clock(it3)
    acctime(1) = acctime(1) + (it2 - it1)*systime_factor
    acctime(2) = acctime(2) + (it3 - it2)*systime_factor
    CALL system_clock(it1)
  END IF

! half non-local step

  IF (ipsptyp == 1 .AND. tnonlocany) THEN
    tenerg = tmodbuf(it,jinfo) .OR. tmodbuf(it,jstinf) .OR. &
             tmodbuf(it,jenergy)
#if(omp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nb,tenerg) SCHEDULE(STATIC)
    DO nb = 1, nstate
      ithr = OMP_GET_THREAD_NUM()
      CALL nonlocstep(q0(1, nb), q1(1, ithr), q2(1, ithr), dt, tenerg, nb, 6) ! 4
    END DO
!$OMP END PARALLEL DO
#else
    DO nb = 1, nstate
      CALL nonlocstep(q0(1, nb), q1, q2, dt, tenerg, nb, 6) ! 4
    END DO
#endif
  END IF

  IF (ttimestop) CALL system_clock(it2)

! one full time step for the kinetic energy

  ithr = 0
#if(dynomp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(q3D,nb,ishift,ithr) SCHEDULE(STATIC)
#endif
  DO nb = 1, nstate
#if(dynomp)
    ithr = OMP_GET_THREAD_NUM()
    IF(nb==1) THEN
      ntm = OMP_GET_MAX_THREADS() 
      nta = OMP_GET_NUM_THREADS()
      IF(ntm.NE.nta) WRITE(*,*) &
        'DYNAMIC: threads not exhausted: nb,ithr,maxthreads,actthreads=',&
         nb,ithr,ntm,nta
    END IF
#else
    ithr = 0
#endif
#if(fftw_cpu)
    q3D(1:kxbox, 1:kybox, 1:kzbox) => q0(1:kdfull2, nb)
    CALL kinprop(q3D)
    NULLIFY (q3D)
#else
    CALL kinprop(q0(1, nb), q1(1, ithr))
#endif
  END DO
#if(dynomp)
!$OMP END PARALLEL DO
#endif
  IF (ttimestop) CALL system_clock(it3)

  FLUSH  (7)

! half non-local step

  IF (ipsptyp == 1 .AND. tnonlocany) THEN
    tenerg = tmodbuf(it,jinfo) .OR. tmodbuf(it,jstinf) .OR. &
             tmodbuf(it,jenergy)

#if(omp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nb,tenerg) SCHEDULE(STATIC)
    DO nb = 1, nstate
      ithr = OMP_GET_THREAD_NUM()
      CALL nonlocstep(q0(1, nb), q1(1, ithr), q2(1, ithr), dt, tenerg, nb, 6) ! 4
    END DO
!$OMP END PARALLEL DO
#else
    DO nb = 1, nstate
      CALL nonlocstep(q0(1, nb), q1, q2, dt, tenerg, nb, 6) ! 4
    END DO
#endif
  END IF
  IF (ttimestop) CALL system_clock(it4)

!DEALLOCATE(q1)
  DEALLOCATE (q2)

#if(fsic)
  IF (tnearest .AND. ifsicp .GE. 8) THEN
    ALLOCATE (qwork(kdfull2, kstate))
    qwork = q0
    CALL eval_unitrot(q0, qwork)
    DEALLOCATE (qwork)
  END IF
#endif

! New density and local potential (this is already the density
! at the end of the step, because it is unchanged by unitary
! potential step) propagation of substrate dipoles is done in
! 'dyn_mfield'.

  CALL dyn_mfield(rho, aloc, q0, dt, it)
  IF (ttimestop) THEN
    CALL system_clock(it5)
    acctime(3) = acctime(3) + (it2 - it1)*systime_factor
    acctime(4) = acctime(4) + (it3 - it2)*systime_factor
    acctime(3) = acctime(3) + (it4 - it3)*systime_factor
    acctime(5) = acctime(5) + (it5 - it4)*systime_factor
    CALL system_clock(it1)
  END IF

#if(extended)
! re-open workspace
!ALLOCATE(q1(2*kdfull2,0:0))
#endif

! half time step incoordinate space:

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ind) SCHEDULE(STATIC)
  DO ind = 1, nlocact !nup
#if(extended)
! pr=-dt*aloc(ind)
! q1(ind,0)=CMPLX(COS(pr),SIN(pr),DP)
#endif
    q1(ind, 0) = EXP(CMPLX(0D0, -dt*aloc(ind), DP))
  END DO
!$OMP END PARALLEL DO

  IF (ttimestop) CALL system_clock(it2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nb,ishift) SCHEDULE(STATIC)
  DO nb = 1, nstate
    ishift = (ispin(nrel2abs(nb)) - 1)*nxyz
    q0(:, nb) = q1(ishift + 1:ishift + kdfull2, 0)*q0(:, nb)
  END DO
!$OMP END PARALLEL DO

  IF (ttimestop) THEN
    CALL system_clock(it3)
    acctime(1) = acctime(1) + (it2 - it1)*systime_factor
    acctime(2) = acctime(2) + (it3 - it2)*systime_factor
    WRITE (6, '(a,10(1pg11.3))') &
      'systime1: prep_loc,loc,nonloc,kin,mf=', acctime(1:5)
  END IF

!DEALLOCATE(q1)

  CALL cpu_time(time_fin)
  time_cpu = time_fin - time_init
  CALL system_clock(ncount_fin, ncount_rate, ncount_max)
  ncount_syst = ncount_fin - ncount_init
  IF (myn == 0) THEN
    acctime(6) = acctime(6) + ncount_syst*systime_factor
    acctime(7) = acctime(7) + 1D0
    WRITE (6, '(a,3(1pg13.5))') ' CPU time inTSTEP', time_cpu, &
      ncount_syst*systime_factor, acctime(6)/acctime(7)
    WRITE (7, '(a,3(1pg13.5))') ' CPU time inTSTEP', time_cpu, &
      ncount_syst*systime_factor, acctime(6)/acctime(7)
    FLUSH  (6)
    FLUSH  (7)
  END IF

#if(extended)
  IF (izforcecorr /= -1) CALL checkzeroforce(rho, aloc)

  IF (tmodbuf(it,jescmask) .OR. tmodbuf(it,jescmaskorb)) CALL escmask(it)
#endif

  FLUSH  (6)
  FLUSH  (7)

  IF (ttimestop) WRITE (6, '(a,2(1pg13.5))') 'times accumulated:', acctime

  RETURN
END SUBROUTINE tstep

!-----dyn_mfield---------------------------------------------------

SUBROUTINE dyn_mfield(rho, aloc, psi, dt, it)

! Computation of the mean field.

! Input/Output:
! psi = set of COMPLEX s.p. wavefunctions (changed in case of SIC)
! Input:
! dt = actual time step
! it = current iteration  
! Output:
! rho = electron density
! aloc = local mean-field potential

  USE params
#if(fsic)
  USE twost
#endif

  IMPLICIT NONE

  REAL(DP), INTENT(OUT) :: rho(2*kdfull2)
  REAL(DP), INTENT(OUT) :: aloc(2*kdfull2)
  COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
  REAL(DP), INTENT(IN) :: dt
  INTEGER, INTENT(IN) :: it

! timer variables
  LOGICAL, PARAMETER :: ttimestop = .FALSE.
  INTEGER :: it1, it2, it3, it4, it5

!----------------------------------------------------------------

  IF (ttimestop) CALL system_clock(it1)
  CALL calcrho(rho, psi)
  IF (ttimestop) CALL system_clock(it2)
  CALL coul_mfield(rho)
  IF (ttimestop) CALL system_clock(it3)

  WRITE (*, *) 'DYN_MFIELD: rho', dvol*SUM(rho(1:kdfull2)) !test

#if(raregas)
  IF (isurf /= 0 .AND. NE > 0) THEN
    IF (dt == 0D0) THEN
      CALL valence_step(rho, dt, it, .false.)
    ELSE
      CALL valence_step(rho, dt, it, .true.)
    END IF
  END IF
#endif

! LDA part of local potential
  CALL calclocal(rho, aloc)
  IF (ttimestop) CALL system_clock(it4)

! SIC
  IF (ifsicp > 0 .AND. ifsicp <= 6) THEN
    CALL calc_sic(rho, aloc, psi)
#if(fsic)
  ELSE IF (ifsicp >= 7) THEN
    CALL calc_utwfc(psi, psiut, NINT(tfs/(dt1*0.0484D0)))
    IF (ifsicp == 7) THEN
      ifsicp = 3
      CALL calc_sic(rho, aloc, psiut)
      ifsicp = 7
    ELSE IF (ifsicp .GE. 8) THEN
      CALL calc_fullsic(psiut, qnewut)
    END IF
#endif
  ELSE IF (ifsicp == 6) THEN
    STOP ' that kind of SIC not valid for dynamics'
  END IF

  IF (ttimestop) THEN
    CALL system_clock(it5)
    WRITE (6, '(a,i4,10(1pg11.3))') 'systime2: iteration,rho,coul,KS,SIC=', &
      it,(it2 - it1)*systime_factor, (it3 - it2)*systime_factor, &
      (it4 - it3)*systime_factor, (it5 - it4)*systime_factor
  END IF

  RETURN
END SUBROUTINE dyn_mfield

!-----boost--------------------------------------------------------

SUBROUTINE boost(q0)

! Boosts electronic wavefunctions 'q0' by 'centfx','centfy','centfz'.
!
! Type of action determined by 'tspindip'.
! tspindip==0 --> c.m. boost, the same for all electrons.
! tspindip==1 --> spin-dipole boost, spin-up and spin-down
! in opposite directions

  USE params
  IMPLICIT NONE

  INTEGER :: ind, ix, iy, iz, nbe
  REAL(DP) :: aspin, actsx, actsy, actsz
  REAL(DP) :: x1, y1, z1
  COMPLEX(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)

!--------------------------------------------------------------------

  DO nbe = 1, nstate

    IF (tspindip) THEN
      aspin = 3 - 2*ispin(nbe)
      actsz = 0.5D0*centfz*aspin
      actsy = 0.5D0*centfy*aspin
      actsx = 0.5D0*centfx*aspin
      WRITE (6, *) ' state nr.', nbe, ' : boost with spin=', aspin
      WRITE (7, *) ' state nr.', nbe, ' : boost with spin=', aspin
    ELSE
      actsz = centfz
      actsy = centfy
      actsx = centfx
    END IF
    ind = 0
    DO iz = minz, maxz
      z1 = (iz - nzsh)*dz*actsz
      DO iy = miny, maxy
        y1 = (iy - nysh)*dy*actsy
        DO ix = minx, maxx
          x1 = (ix - nxsh)*dx*actsx
          ind = ind + 1
          q0(ind, nbe) = EXP(CMPLX(0D0, x1 + y1 + z1, DP))*q0(ind, nbe)
        END DO
      END DO
    END DO

  END DO

  RETURN
END SUBROUTINE boost

!-----info ------------------------------------------------------------

SUBROUTINE info(psi, rho, aloc, it, lprint)

! Print information on energy observables: single particle energies,
! kinetic energy, total energy, ionic energy, ...
!
! Input:
! psi = set of COMPLEX s.p. wavefunctions
! rho = electron density
! aloc = local mean-field potential
! it = time step in calling routine
! lprint = additional switch for printing s.p. energies

  USE params
  USE kinetic, ONLY: akv, calc_ekin
  USE util, ONLY: wfovlp, safeopen, project, tmodbuf
  IMPLICIT NONE

#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
  REAL(DP) :: enonlcp, esh1p, eshellp
#endif

  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
  REAL(DP), INTENT(IN) :: rho(2*kdfull2)
  REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
  INTEGER, INTENT(IN) :: it
  LOGICAL, INTENT(IN) :: lprint

  INTEGER :: ind, ion, ishift, iss, nb, nbe
  REAL(DP) :: ekin, ehilf, eshell, enonlc
  REAL(DP) :: ek, tinfs, xm
  REAL(DP) :: en(kstate)
  COMPLEX(DP), ALLOCATABLE :: qtmp(:)
  COMPLEX(DP), ALLOCATABLE :: psitmp(:, :)
  REAL(DP), ALLOCATABLE :: current(:, :)
  REAL(DP) :: estar(2), estarETF(2)
!  COMPLEX(DP), ALLOCATABLE :: psiaux(:)

  LOGICAL :: topenf
  LOGICAL, PARAMETER :: ttest = .FALSE.
  LOGICAL, PARAMETER :: ttesthpsi = .FALSE.
#if(extended)
  LOGICAL, PARAMETER :: tvar_project = .FALSE. ! s.p variance projected (expensive)
#endif

  REAL(DP), EXTERNAL :: energ_ions
  REAL(DP), EXTERNAL:: enerkin_ions ! declared inion_md
#if(raregas)
  REAL(DP), PARAMETER :: alpha_ar = 10.6D0 ! for VdW
  INTEGER :: ico
#endif

! variables and parameter for timer
  INTEGER :: itimeinfo(11),irateinfo
  LOGICAL,PARAMETER :: ttimeinfo=.FALSE.     ! detailed timing information

  IF(ttimeinfo) THEN
    CALL system_clock(itimeinfo(1), irateinfo)
    CALL system_clock(itimeinfo(1))
  END IF

#if(mpi)
  CALL mpi_comm_rank(mpi_comm_world, myn, mpi_ierror)
#else
  myn = 0
#endif
  tinfs = it*dt1*0.0484D0/2.0D0/ame
  tfs = it*dt1*0.0484D0
#if(raregas)
  IF (idielec == 1) THEN
    CALL energ_dielec(rho)
  ELSE
    ecrhoimage = 0D0
  END IF
#endif

! compute single-particle energies and related quantities

  IF (ifsicp == 5) psisavex = psi

  eshell = 0D0
  esh1 = 0D0
  enonlc = 0D0

#if(dynomp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ishift) SCHEDULE(STATIC) 
#endif
  DO nb = 1, nstate
    spnorm(nb) = wfovlp(psi(:, nb), psi(:, nb))
    CALL calc_ekin(psi(:, nb), ekinsp(nb))
    ishift = (ispin(nrel2abs(nb)) - 1)*nxyz
    CALL calc_epot(psi(:, nb), aloc(ishift + 1), epotsp(nb), enonlo(nb), nb)
  END DO
#if(dynomp)
!$OMP END PARALLEL DO
#endif
  IF(ttimeinfo) CALL system_clock(itimeinfo(2))

  DO nb = 1, nstate
    ekin = ekinsp(nb)
    ehilf = epotsp(nb)
    ehilf = ehilf + ekin
    amoy(nb) = ehilf
    epot = ehilf + ekin
    en(nb) = epot*occup(nb)
    eshell = eshell + en(nb)
    esh1 = esh1 + ekin*occup(nb)
    enonlc = enonlc + enonlo(nb)*occup(nb)

#if(!mpi)
    WRITE (6, '(2(a,i2),5(a,f9.5))') 'level:', nrel2abs(nb), &
      ' ispin=', ispin(nrel2abs(nb)), ' occup=', occup(nb), &
      ' ekin=' &
      , ekin, ' epot=', epotsp(nb), ' esp=', amoy(nb), ' enonlo=', enonlo(nb)
!    IF (ttesthpsi) THEN
!      ALLOCATE (psiaux(kdfull2))
!      psiaux = psi(1:nxyz, nb)
!      CALL hpsi(psiaux, aloc(ishift + 1), nb, 1)
!      DEALLOCATE (psiaux)
!    END IF
#endif
  END DO

  IF(ttimeinfo) CALL system_clock(itimeinfo(3))

  IF(lprint) CALL spmoms(psi, 6)

  IF(ttimeinfo) CALL system_clock(itimeinfo(4))

  IF(jstinf > 0) THEN
    tstinf = MOD(it, jstinf) == 0 .AND. lprint
  ELSE
    tstinf = .FALSE.
  END IF
  IF (tstinf .AND. ifsicp == 5) psisavex = psi
  IF (tstinf .OR. ttesthpsi) THEN
    ALLOCATE (qtmp(kdfull2))
    DO nb = 1, nstate
      qtmp = psi(:, nb)
      ishift = (ispin(nrel2abs(nb)) - 1)*nxyz
      CALL hpsi(qtmp, aloc(ishift + 1), nb, 1)
#if(extended)
#if(!mpi)
      IF (tvar_project) THEN
        CALL project(qtmp, qtmp, ispin(nb), psi)
        spvariancep(nb) = SQRT(REAL(wfovlp(qtmp, qtmp), DP))
      END IF
#endif
#endif
    END DO
    DEALLOCATE (qtmp)


#if(!mpi)
    CALL safeopen(77, it, jstinf, 'pspenergies')
    WRITE (77, '(1f15.6,500f12.6)') tfs, (amoy(nb), nb=1, nstate)
    FLUSH  (77)
    CALL safeopen(76, it, jstinf, 'pspvariances')
    WRITE (76, '(1f15.6,500f12.6)') tfs, (spvariance(nb), nb=1, nstate)
    FLUSH  (76)
#if(extended)
    IF (tvar_project) THEN
      CALL safeopen(75, it, jstinf, 'pspvariancesp')
      WRITE (75, '(1f15.6,500f12.6)') tfs, (spvariancep(nb), nb=1, nstate)
      FLUSH  (75)
    END IF
#endif
#endif
  END IF
  IF(ttimeinfo) CALL system_clock(itimeinfo(5))

#if(extended)
  IF (tmodbuf(it, jstboostinv)) THEN
    IF (.NOT. ALLOCATED(akv)) STOP "jstboostinv NOT YET FOR FINITE DIFFERENCE"
    ALLOCATE (current(kdfull2, 3))
    ALLOCATE (qtmp(kdfull2))
    CALL calc_current(current, psi)
    DO ind = 1, kdfull2
      current(ind, :) = current(ind, :)/MAX(rho(ind), 1D-5)
    END DO
    DO nb = 1, nstate
      qtmp = psi(:, nb)
      ishift = (ispin(nrel2abs(nb)) - 1)*nxyz
      CALL hpsi_boostinv(qtmp, aloc(ishift + 1), current, rho, nbe)
      spenergybi(nb) = REAL(wfovlp(psi(:, nb), qtmp), DP)
      spvariancebi(nb) = SQRT(REAL(wfovlp(qtmp, qtmp), DP) - spenergybi(nb)**2)
    END DO
    IF (lprint) THEN
      CALL safeopen(91, it, jstboostinv, 'pspenergybi')
      WRITE (91, '(1f15.6,500f12.6)') tfs, (spenergybi(nb), nb=1, nstate)
      CLOSE (91)
      CALL safeopen(92, it, jstboostinv, 'pspvariancesbi')
      WRITE (92, '(1f15.6,500f12.6)') tfs, (spvariancebi(nb), nb=1, nstate)
      CLOSE (92)
    END IF
    DEALLOCATE (current, qtmp)
  END IF
#endif

#if(mpi)
  CALL prispe_parallele(6, it)
  IF (ttest) WRITE (*, *) ' INFO: before allreduce. myn=', myn
  CALL mpi_barrier(mpi_comm_world, mpi_ierror)
  CALL mpi_allreduce(eshell, eshellp, 1, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, mpi_ierror)
  CALL mpi_allreduce(esh1, esh1p, 1, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, mpi_ierror)
  CALL mpi_allreduce(enonlc, enonlcp, 1, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, mpi_ierror)
  CALL mpi_barrier(mpi_comm_world, mpi_ierror)
  IF (ttest) WRITE (*, *) ' INFO: after allreduce. myn=', myn
  esh1 = esh1p
  eshell = eshellp
  enonlc = enonlcp
#endif

  IF (myn == 0 .AND. jenergy > 0 .AND. lprint) THEN;
  IF (MOD(it, jenergy) == 0) THEN
    DO iss = 2, 1, -1
      CALL calc_estar(psi, iss, estar(iss), estarETF(iss))
    END DO

    IF (lprint) THEN
      WRITE (7, *) 'ekintot=', esh1
      WRITE (7, *) 'estar(1)=', estar(1), ' estar(2)=', estar(2)
    END IF
    WRITE (*, *) 'estar(1:2)=', estar, 'estarETF(1:2)=', estarETF
  END IF; END IF

  IF(ttimeinfo) CALL system_clock(itimeinfo(6))

! rearrangement and background Coulomb energy
  ecback = 0D0
  ecrho = 0D0
  IF (nion2 /= 0) THEN
    DO ind = 1, nxyz
      ecback = ecback - rho(ind)*potion(ind)
      ecrho = ecrho + rho(ind)*(chpcoul(ind) - potion(ind))
    END DO
  ELSE ! jellium CASE
    DO ind = 1, nxyz
      ecback = ecback - rhojel(ind)*chpcoul(ind)
      ecrho = ecrho + rho(ind)*chpcoul(ind)
    END DO
  END IF
  ecback = ecback*dvol/2.0D0
  ecrho = ecrho*dvol/2.0D0

#if(raregas)
  IF (nc > 0 .AND. ivdw == 1) THEN
    esub = 0D0
    IF (nion2 /= 0) THEN
      DO ind = 1, nxyz
        esub = esub + rho(ind)*potvdw(ind)
      END DO
    END IF
    esub = esub*dvol
    evdw = esub
    DO iss = 1, nc
      DO ico = 1, 3
        evdw = evdw - 0.5D0*e2*alpha_ar*frho(iss, ico)*frho(iss, ico)/(REAL(nelect, DP))
      END DO
    END DO
    eshell = eshell - esub
  END IF
#endif

  esh1 = esh1
  eshell = eshell/2.0D0 !(=t+v/2)

! ionic contributions to the energy

  ecorr = energ_ions()

!  ekion = 0D0 ! kinetic energy of Na cores
#if(raregas)
  ekinion = 0D0 ! kinetic energy of GSM cores
  ekinel = 0D0 ! kinetic energy of GSM shells
  ekinkat = 0D0 ! kinetic energy of cations
#endif
  IF (ionmdtyp > 0) THEN
    ekion = enerkin_ions()
#if(extended)
!    WRITE(*,*) 'EKION before/after:',xm,ekion
!    DO ion = 1, nion
!      ek = cpx(ion)*cpx(ion) + cpy(ion)*cpy(ion) + cpz(ion)*cpz(ion)
!      xm = 1836.0D0*amu(np(ion))*ame
!      ek = ek/2D0/xm
!      ekion = ekion + ek
!    END DO
#endif
#if(raregas)
    DO ion = 1, nc
      ek = pxc(ion)*pxc(ion) + pyc(ion)*pyc(ion) + pzc(ion)*pzc(ion)
      xm = 1836.0D0*mion*ame
      ek = ek/2.0D0/xm
      ekinion = ekinion + ek
    END DO
    DO ion = 1, NE
      ek = pxe(ion)*pxe(ion) + pye(ion)*pye(ion) + pze(ion)*pze(ion)
      xm = 1836.0D0*me*ame
      ek = ek/2.0D0/xm
      ekinel = ekinel + ek
    END DO
    DO ion = 1, nk
      ek = pxk(ion)*pxk(ion) + pyk(ion)*pyk(ion) + pzk(ion)*pzk(ion)
      xm = 1836.0D0*mkat*ame
      ek = ek/2.0D0/xm
      ekinkat = ekinkat + ek
    END DO
#endif
  END IF

! final composition and PRINT

  energy = eshell + enrear + ecback + ecorr + enonlc/2D0 + ecrhoimage
#if(raregas)
  IF (ivdw == 1) energy = energy + evdw
  ecoul = ecback + ecrho + ecorr + ecrhoimage
  etot = energy + ekion + ekinion + ekinel + ekinkat
#else
  etot = energy + ekion
#endif

  energ2 = esh1 + enerpw + ecrho + ecback + ecorr + enonlc - ecrhoimage
  IF (myn == 0 .AND. jenergy > 0 .AND. lprint) THEN;
  IF (MOD(it, jenergy) == 0) THEN
    CALL safeopen(163, it, jenergy, 'penergies')
    WRITE(*,*) 'EKION after safeopen:',ekion
    WRITE (163, '(1f12.6,25(1pg16.8))') tfs, &
      etot, &
      energy, &
      eshell*2.-esh1, &
      enrear, &
      ekion, &
      ecorr, &
      2D0*ecback, &
      enonlc, &
      2D0*ecback + ecorr + enonlc, &
      datalaser(2), &
      estar, &
#if(raregas)
      estarETF, &
      enii, &
      enig, &
      engg, &
      ekinion, &
      ekinel, &
      ekinkat &
#else
      estarETF
#endif
    FLUSH  (163)
    CLOSE (163)
  END IF; END IF

#if(mpi)
  IF (myn == 0) THEN
#endif
    WRITE (6, '(a)') ' '
    WRITE (6, *) 'tot sp energ = ', eshell*2 - esh1
    WRITE (6, *) 'rearge. energ = ', enrear
    WRITE (6, *) 'e_coul:ion-ion = ', ecorr
    WRITE (6, *) 'e_coul:el-ion = ', 2*ecback
    WRITE (6, *) 'extern. energy = ', 2*ecback + ecorr
    WRITE (6, *) 'Hartree energy = ', ecrho - ecback - ecrhoimage
    WRITE (6, *) 'nonlocal energy = ', enonlc
    WRITE (6, *) 'sim.ann.energy = ', 2*ecback + ecorr + enonlc
    WRITE (6, *) 'laser energy = ', datalaser(2)
    WRITE (6, *) 'internal exc. energy per spin = ', estar(1), estar(2)
#if(raregas)
    IF (idielec == 1) WRITE (6, *) 'image energy = ', ecrhoimage
    IF (ivdw == 1) WRITE (6, *) 'vdw energy = ', evdw
    IF (isurf /= 0) THEN
      WRITE (6, *) 'adsorb.energy = ', etot - enerinfty
    END IF
#endif
    WRITE (6, *) 'binding energy = ', energy
    WRITE (6, *) 'total energy = ', etot

    IF (tmodbuf(it,jinfo) .and. lprint) THEN
      INQUIRE (17, OPENED=topenf)
      IF (.NOT. topenf) OPEN (17, POSITION='append', FILE='infosp.'//outnam)
      WRITE (17, '(a,i5,f8.3,1pg13.5)') &
       ' time_step,time,energy=', it, tinfs, etot
      CLOSE (17)
    END IF

#if(mpi)
  END IF
#endif

  IF(ttimeinfo) THEN
    CALL system_clock(itimeinfo(7))
    WRITE(*,'(a,i12,6(1pg13.5))') 'INFO timer:', irateinfo,&
        ((itimeinfo(2)-itimeinfo(1))*1D0/irateinfo),&
        ((itimeinfo(3)-itimeinfo(2))*1D0/irateinfo),&
        ((itimeinfo(4)-itimeinfo(3))*1D0/irateinfo),&
        ((itimeinfo(5)-itimeinfo(4))*1D0/irateinfo),&
        ((itimeinfo(6)-itimeinfo(5))*1D0/irateinfo),&
        ((itimeinfo(7)-itimeinfo(6))*1D0/irateinfo)
  END IF
  IF (myn == 0) THEN
    WRITE (6, '(a,f8.4,a,f7.2,a,3f10.5)') &
      'in info: t=', tfs, ' moments: monop.=', qe(1), &
      ' dip.: ', qe(2), qe(3), qe(4)
    WRITE (6, '(a,f8.4,a,3f10.5)') &
      'in info: t=', tfs, ' moments: quads=', qe(5), qe(6), qe(7)
  END IF

  tstinf = .false.



#if(extended)
! OPTIONALLY scan 2Ph matrix elements (activate and set parameters by hand)
!  IF (it > 0 .AND. MOD(it, 50) == 0) CALL V2ph(psi, 1.3D0, 0D0, 1D0)
#endif

  RETURN
END SUBROUTINE info

!-----calc_epot------------------------------------------------------------

SUBROUTINE calc_epot(psin, alocact, epotout, enonlocout, nb)

! Calculates potential energy for state 'nb'.

! Input:
! psi = wavefunction of state 'nb'
! alocact = local mean-field potential
! nb = nr. of s.p. state
!
! Output:
! enonlocout = controbution from non-local part
! epot = contribution from local part,
! via module 'params'

  USE params
  USE util, ONLY: wfovlp
#if(fsic)
  USE twost
  USE twostr, ONLY: ndims
#endif

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN) :: psin(kdfull2)
  REAL(DP), INTENT(IN) :: alocact(2*kdfull2)
  REAL(DP), INTENT(OUT) :: enonlocout, epotout
  INTEGER, INTENT(IN) :: nb

  INTEGER :: i
  REAL(DP) :: sumnon, accepot
  COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: psi2, qex
  INTEGER :: is, na, nbe, nae
  COMPLEX(DP) :: cf

  LOGICAL, PARAMETER :: ttest = .false.

!------------------------------------------------------------------

  ALLOCATE (psi2(kdfull2))

  accepot = 0D0

  IF (ttest) WRITE (*, *) &
    ' in CALC_EPOT 1. average field=', SUM(alocact(1:nxyz))

! non local part of ps

  IF (ipsptyp == 1) THEN
    CALL nonlocalc(psin, psi2, 0)
    sumnon = 0D0
    DO i = 1, nxyz
      sumnon = sumnon + REAL(psin(i), DP)*REAL(psi2(i), DP) &
               + AIMAG(psin(i))*AIMAG(psi2(i))
    END DO
    enonlocout = sumnon*dvol
    DO i = 1, nxyz
      psi2(i) = psi2(i) + alocact(i)*psin(i)
    END DO
  ELSE
    DO i = 1, nxyz
      psi2(i) = alocact(i)*psin(i)
    END DO
  END IF

  IF (ttest) WRITE (*, *) ' in CALC_EPOT 2'
! subtract SIC potential for state NB
#if(fsic)
  IF (ifsicp .GE. 8) THEN
    is = ispin(nrel2abs(nb))
    DO na = 1, ndims(is)
      nbe = nb - (is - 1)*ndims(1)
      nae = na + (is - 1)*ndims(1)
      WRITE (*, *) ' nb,is,na,nbe,vecs=', nb, is, na, nbe, vecs(nbe, na, is)
      cf = CONJG(vecs(nbe, na, is))
      WRITE (*, *) ' nb,is,na,nbe,cf=', nb, is, na, nbe, cf
      DO i = 1, nxyz
        psi2(i) = psi2(i) - qnewut(i, nae)*cf
      END DO
    END DO
  END IF
#endif

  IF (ifsicp == 5) THEN
    ALLOCATE (qex(kdfull2))
    CALL exchg(psin, qex, nb)
    psi2 = psi2 + qex
    DEALLOCATE (qex)
  END IF

  IF (ttest) WRITE (*, *) ' in CALC_EPOT 3'
  accepot = REAL(wfovlp(psin, psi2)) ! develop real overlap
  epotout = accepot

  IF (ttest) WRITE (*, *) ' EPOT: nb,epot=', nb, accepot

  DEALLOCATE (psi2)

  RETURN
END SUBROUTINE calc_epot

!-----calc_estar------------------------------------------------

SUBROUTINE calc_estar(psin, iss, excit, excitETF)

! Compute kinetic excitation energy of the electron cloud
! as E* = Ekin - Ekin(Thomas-Fermi).
!
! Input:
! psin = set of s.p. wavefunctions to be analyzed
! iss = spin for which analysis is done
!
! Output:
! excit = excitation energy with mere TF subtraction
! excitETF = excitation energy with full ETF subtraction

  USE params
  USE kinetic
  IMPLICIT NONE
#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
#endif

  COMPLEX(DP), INTENT(IN) :: psin(kdfull2, kstate)
  INTEGER, INTENT(IN) :: iss
  REAL(DP), INTENT(OUT) :: excit, excitETF

  REAL(DP), DIMENSION(:), ALLOCATABLE :: arho, gradrho
  COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: gradrhok, gradrhoksav
#if(mpi)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: arhop
  REAL(DP) :: ekintotp
#endif
  REAL(DP) :: factgrad

  LOGICAL, PARAMETER :: extendedTF = .TRUE. ! switch to ETF
  INTEGER :: i, j, nb
  REAL(DP) :: anorm, ekintot
  REAL(DP), PARAMETER :: rholimitrel = 1D-5 ! relative lower limit of density (avoid overflow)
  REAL(DP), PARAMETER :: rholimitlog = 1D-3 ! lower limit in LOG(arho)
  REAL(DP) :: rholimit

!------------------------------------------------------------

  ALLOCATE (arho(kdfull2))
#if(mpi)
  ALLOCATE (arhop(kdfull2))
#endif

#if(mpi)
  CALL mpi_barrier(mpi_comm_world, mpi_ierror)
#endif

  arho = 0D0
  ekintot = 0D0
  DO nb = 1, nstate
    j = nrel2abs(nb)
    IF (ispin(j) .NE. iss) THEN
      CYCLE
    ELSE
      ekintot = ekintot + occup(nb)*ekinsp(nb)
      DO i = 1, kdfull2
        arho(i) = arho(i) + &
                  occup(nb)*(REAL(psin(i, nb), DP)**2 + AIMAG(psin(i, nb))**2)
      END DO
    END IF
  END DO

  rholimit = MAXVAL(arho)*rholimitrel

#if(mpi)
  CALL mpi_barrier(mpi_comm_world, mpi_ierror)
  CALL mpi_allreduce(arho, arhop, kdfull2, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, mpi_ierror)
  CALL mpi_allreduce(ekintot, ekintotp, 1, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, mpi_ierror)
  CALL mpi_barrier(mpi_comm_world, mpi_ierror)
#endif

#if(mpi)
  arho = arhop
  DEALLOCATE (arhop)
#endif

  excit = 0D0
  DO i = 1, kdfull2
    excit = excit + arho(i)**(5D0/3D0)
  END DO
  excit = excit*dvol*h2m*0.6D0*(6.0D0*pi**2)**(2.D0/3.D0)

  IF (extendedTF) THEN

    ALLOCATE (gradrho(kdfull2), gradrhok(kdfull2), gradrhoksav(kdfull2))
    factgrad = dvol*h2m/18.0D0
#if(extended)
! factgrad=dvol*h2m/36.0D0
#endif
    excitETF = excit

! x derivative
    DO i = 1, kdfull2
      gradrho(i) = LOG(MAX(arho(i), rholimit*rholimitlog))
    END DO
    CALL rftf(gradrho, gradrhoksav)
#if(netlib_fft|fftw_cpu)
    CALL gradient(gradrhoksav, gradrhok, 1)
    CALL rfftback(gradrhok, gradrho)
#endif
    DO i = 1, kdfull2
      IF (arho(i) .gt. rholimit) &
        excitETF = excitETF + factgrad*(gradrho(i))**2*arho(i)
    END DO

! y derivative
#if(netlib_fft|fftw_cpu)
    CALL gradient(gradrhoksav, gradrhok, 2)
    CALL rfftback(gradrhok, gradrho)
#endif
    DO i = 1, kdfull2
      IF (arho(i) .gt. rholimit) &
        excitETF = excitETF + factgrad*(gradrho(i))**2*arho(i)
    END DO

! z derivative
#if(netlib_fft|fftw_cpu)
    CALL gradient(gradrhoksav, gradrhok, 3)
    CALL rfftback(gradrhok, gradrho)
#endif
    DO i = 1, kdfull2
      IF (arho(i) .gt. rholimit) &
        excitETF = excitETF + factgrad*(gradrho(i))**2*arho(i)
    END DO

  END IF

  anorm = SUM(arho)*dvol

  DEALLOCATE (arho)
  DEALLOCATE (gradrho)
  DEALLOCATE (gradrhok, gradrhoksav)

  WRITE (7, '(a,7(1pg13.5))') ' norm,TF,ETF=', anorm, excit, excitETF

#if(mpi)
  excitETF = ekintotp - excitETF
  excit = ekintotp - excit
#else
  excitETF = ekintot - excitETF
  excit = ekintot - excit
#endif

  WRITE (7, '(a,2i5,3(1pg13.5))') &
    ' CALC_ESTAR: myn,iss,ekintot,excit,excitETF', &
    myn, iss, ekintot, excit, excitETF

  RETURN
END SUBROUTINE calc_estar

!-----mrote ------------------------------------------------------------

SUBROUTINE mrote

! Rotates ionic configuration by 'phirot' around axis 'irotat:
! irotat=1 -> rotate around x-axis
! irotat=2 -> rotate around y-axis
! irotat=3 -> rotate around z-axis
! irotat=4 -> rotate around diagonal axis
! The angle 'phirot' is to be given in degree.
! Configuration and parameters are communicated via module 'params'

  USE params
  USE util, ONLY: rotatevec3D
  IMPLICIT NONE

  INTEGER :: ion
  REAL(DP) :: vecin(3), vecout(3), vecalpha(3)

!------------------------------------------------------------------

! determine vector of rotation angles (in radian)
  IF (irotat == 1) THEN
    vecalpha(1) = phirot*PI/180D0
    vecalpha(2) = 0D0
    vecalpha(3) = 0D0
  ELSE IF (irotat == 2) THEN
    vecalpha(1) = 0D0
    vecalpha(2) = phirot*PI/180D0
    vecalpha(3) = 0D0
  ELSE IF (irotat == 3) THEN
    vecalpha(1) = 0D0
    vecalpha(2) = 0D0
    vecalpha(3) = phirot*PI/180D0
  ELSE IF (irotat == 4) THEN
    vecalpha(1) = phirot*PI/180D0
    vecalpha(2) = vecalpha(1)
    vecalpha(3) = vecalpha(1)
  END IF

  OPEN (79, FILE='rotated-config.res')
  WRITE (79, '(a)') ' ion x y z'

  DO ion = 1, nion

    vecin(1) = cx(ion)
    vecin(2) = cy(ion)
    vecin(3) = cz(ion)

    CALL rotatevec3D(vecin, vecout, vecalpha)

    cx(ion) = vecout(1)
    cy(ion) = vecout(2)
    cz(ion) = vecout(3)

    WRITE (7, '(/a,i4)') 'ion no:', ion
    WRITE (7, '(a,f12.4)') 'cxnew=', cx(ion)
    WRITE (7, '(a,f12.4)') 'cynew=', cy(ion)
    WRITE (7, '(a,f12.4)') 'cznew=', cz(ion)
    WRITE (79, '(3f12.4)') cx(ion), cy(ion), cz(ion)

  END DO

  RETURN
END

!-----instit----------------------------------------------------

SUBROUTINE instit(psi)

! Calculation of electronic angular momentum.
!
! Input:
! psi = set of s.p. wavefunctions
!
! Output 'ajx','ajy','ajz' via module 'params'

  USE params
  USE kinetic
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)

  INTEGER :: i, ind, ix, iy, iz, nb, occ
  REAL(DP) :: ajpx, ajpy, ajpz, dkx, dky, dkz, x1, y1, z1
  REAL(DP) :: test
  COMPLEX(DP), ALLOCATABLE :: q2(:)
  COMPLEX(DP), ALLOCATABLE :: jtx(:), jty(:), jtz(:)
  COMPLEX(DP) :: jalpha

!------------------------------------------------------------------

  ALLOCATE (q2(kdfull2), jtx(kdfull2), jty(kdfull2), jtz(kdfull2))

  dkx = pi/(dx*REAL(nx, DP))
  dky = pi/(dy*REAL(ny, DP))
  dkz = pi/(dz*REAL(nz, DP))

  DO i = 1, kdfull2
    jtx(i) = CMPLX(0D0, 0D0, DP)
    jty(i) = CMPLX(0D0, 0D0, DP)
    jtz(i) = CMPLX(0D0, 0D0, DP)
  END DO

  DO nb = 1, nstate
    occ = occup(nb)

#if(netlib_fft|fftw_cpu)
    CALL fftf(psi(1, nb), q2)

    DO ind = 1, kdfull2
      q2(ind) = q2(ind)*akx(ind)
    END DO

    CALL fftback(q2, q2)
#endif

    DO ind = 1, kdfull2
      test = eye/2D0*(CONJG(psi(ind, nb))*q2(ind) - psi(ind, nb)*CONJG(q2(ind)))
      jalpha = test
      jtx(ind) = jtx(ind) - occ*jalpha
    END DO

#if(netlib_fft|fftw_cpu)
    CALL fftf(psi(1, nb), q2)

    DO ind = 1, kdfull2
      q2(ind) = q2(ind)*aky(ind)
    END DO

    CALL fftback(q2, q2)
#endif
    DO ind = 1, kdfull2
      test = eye/2D0*(CONJG(psi(ind, nb))*q2(ind) - psi(ind, nb)*CONJG(q2(ind)))
      jalpha = test
      jty(ind) = jty(ind) - occ*jalpha
    END DO

#if(netlib_fft|fftw_cpu)
    CALL fftf(psi(1, nb), q2)

    DO ind = 1, kdfull2
      q2(ind) = q2(ind)*akz(ind)
    END DO

    CALL fftback(q2, q2)
#endif

    DO ind = 1, kdfull2
      test = eye/2D0*(CONJG(psi(ind, nb))*q2(ind) - psi(ind, nb)*CONJG(q2(ind)))
      jalpha = test
      jtz(ind) = jtz(ind) - occ*jalpha
    END DO

  END DO
! END loop over state

  ajx = 0D0
  ajy = 0D0
  ajz = 0D0

  ind = 0
  DO iz = minz, maxz
    z1 = (iz - nzsh)*dz
    DO iy = miny, maxy
      y1 = (iy - nysh)*dy
      DO ix = minx, maxx
        x1 = (ix - nxsh)*dx
        ind = ind + 1
        ajpx = y1*jtz(ind) - z1*jty(ind)
        ajpy = z1*jtx(ind) - x1*jtz(ind)
        ajpz = x1*jty(ind) - y1*jtx(ind)
        ajx = ajx + ajpx
        ajy = ajy + ajpy
        ajz = ajz + ajpz
      END DO
    END DO
  END DO
  ajx = ajx*dvol
  ajy = ajy*dvol
  ajz = ajz*dvol
  WRITE (6, '(a,3f12.4)') 'moments', ajx, ajy, ajz
  WRITE (6, *)

  DEALLOCATE (q2, jtx, jty, jtz)

  RETURN
END SUBROUTINE instit

!-----nonlocstep---------------------------------------------------

SUBROUTINE nonlocstep(qact, q1, q2, ri, tenerg, nb, norder)

! Computes one time step for the non-local potentials
! using exponential evolution.
!
! Input:
! ri = size of time step
! tenerg = (LOGICAL) switch to accumulate non-local energy
! nb = number of state which is propagated
! norder = order of step (up to 6, 4 or 6 recommended)
!
! Input/Output:
! qact = array for actual wavefunction to be propagated
! q1,q2 = auxiliary wavefunction arrays
!
! !! Note: This routine should be replaced by truly exponential
! propagation within the non-local plaquettes.
!
  USE params
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN OUT) :: qact(kdfull2)
  COMPLEX(DP), INTENT(IN OUT) :: q1(kdfull2)
  COMPLEX(DP), INTENT(IN OUT) :: q2(kdfull2)
  REAL(DP), INTENT(IN) :: ri
  LOGICAL, INTENT(IN) :: tenerg
  INTEGER, INTENT(IN) :: nb
  INTEGER, INTENT(IN) :: norder

  INTEGER :: i, ind
  REAL(DP) :: ri2, rfac, sumadd
  COMPLEX(DP) :: rii, cfac

!---------------------------------------------------------------------

  CALL nonlocalc(qact, q1, 0)
  IF (tenerg) THEN ! add nonloc.pot energy
    sumadd = 0D0
    DO i = 1, nxyz
      sumadd = REAL(qact(i), DP)*REAL(q1(i), DP) + &
               AIMAG(qact(i))*AIMAG(q1(i)) + sumadd
    END DO
    enonlo(nb) = sumadd*dvol
    epotsp(nb) = sumadd*dvol + epotsp(nb)
  END IF

  cfac = -ri*eye
  rii = cfac
  DO ind = 1, nxyz
    qact(ind) = qact(ind) + rii*q1(ind)
  END DO

  CALL nonlocalc(q1, q2, 0)
  rii = rii*cfac/2D0
  DO ind = 1, nxyz
    qact(ind) = qact(ind) + rii*q2(ind)
  END DO
  IF (norder <= 2) RETURN

  CALL nonlocalc(q2, q1, 0)
  rii = rii*cfac/3D0
  DO ind = 1, nxyz
    qact(ind) = qact(ind) + rii*q1(ind)
  END DO
  IF (norder <= 3) RETURN

  CALL nonlocalc(q1, q2, 0)
  rii = rii*cfac/4D0
  DO ind = 1, nxyz
    qact(ind) = qact(ind) + rii*q2(ind)
  END DO
  IF (norder <= 4) RETURN

  CALL nonlocalc(q2, q1, 0)
  cfac = -ri2*ri2*ri/30D0*eye
  rii = rii*cfac/5D0
  DO ind = 1, nxyz
    qact(ind) = qact(ind) + rii*q1(ind)
  END DO
  IF (norder <= 5) RETURN

  CALL nonlocalc(q1, q2, 0)
  rii = rii*cfac/6D0
  DO ind = 1, nxyz
    qact(ind) = qact(ind) + rii*q2(ind)
  END DO

  RETURN
END SUBROUTINE nonlocstep

!-----init_dynprotocol-------------------------------------------------

#if(raregas)
SUBROUTINE init_dynprotocol(rho, aloc, psi)
#else
  SUBROUTINE init_dynprotocol()
#endif

! Initializes files for protocol of dynamics.
!
! Header lines with 'H:' in first two columns serve as
! communicators to later off-line spectral analysis.

    USE params
    IMPLICIT NONE

#if(raregas)
    REAL(DP), INTENT(IN) :: rho(2*kdfull2)
    REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
    COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)

    REAL(DP) :: dist
#endif
    LOGICAL :: topenf
    INTEGER :: i, ion
    CHARACTER(LEN=1) :: ext
    REAL(DP), EXTERNAL :: getxval
    REAL(DP), EXTERNAL :: getyval
    REAL(DP), EXTERNAL :: getzval
    REAL(DP) :: scal1,scal2     ! auxiliary for phase corrected PES
    REAL(DP) :: xx,yy,zz,rr

    IF (irest <= 0) THEN ! WRITE FILE headers ONLY for fresh dynamics

#if(extended)
      IF (myn == 0 .AND. jdip /= 0 .AND. eproj /= 0) THEN
        OPEN (88, STATUS='unknown', FORM='FORMATTED', FILE='pprojdip.'//outnam)
        WRITE (88, '(a)') ' & '
        WRITE (88, '(a)') 'x:time'
        WRITE (88, '(a)') 'y:dipole-momenta of projectile'
        WRITE (88, '(a)') 'z:dipole-momenta of TARGET'
        WRITE (88, '(a)') 'n: time in fs dipole-momenta: x y z'
        WRITE (88, '(a)') 'H: X Yl Yd Yq Zl Zd Zq'
        CLOSE (88)
      END IF
#endif

      IF (myn == 0 .AND. jdip /= 0) THEN
        OPEN (8, STATUS='unknown', FORM='FORMATTED', FILE='pdip.'//outnam)
        WRITE (8, '(a)') '# protocol of electronic dipole moments, units a_0'
        WRITE (8, '(a)') '# time[fs] dipole_x dipole_y dipole_z'
        WRITE (8, '(a)') 'H: X YL YD YP'
        CLOSE (8)
      END IF

      IF (myn == 0 .AND. jdiporb /= 0) THEN
        OPEN (810, STATUS='unknown', FORM='FORMATTED', FILE='pdiporb.x.'//outnam)
        WRITE (810, '(a)') 'protocol of s.p. moments: time,x-dipole of orbitals'
        CLOSE (810)
        OPEN (811, STATUS='unknown', FORM='FORMATTED', FILE='pdiporb.y.'//outnam)
        WRITE (811, '(a)') 'protocol of s.p. moments: time,y-dipole of orbitals'
        CLOSE (811)
        OPEN (812, STATUS='unknown', FORM='FORMATTED', FILE='pdiporb.z.'//outnam)
        WRITE (812, '(a)') 'protocol of s.p. moments: time,z-dipole of orbitals'
        CLOSE (812)
      END IF

      IF (jmp /= 0 .AND. myn == 0) THEN
        OPEN (803, STATUS='unknown', FILE='pMP.'//outnam)
        WRITE (803, '(a)') '# nr. of meas.points, nr. of state'
        WRITE (803, '(a,2i6)') '# ', nmps, nstate_all
        DO i = 1, nmps
          xx = getxval(imps(i))
          yy = getyval(imps(i))
          zz = getzval(imps(i))
          rr = SQRT(xx*xx+yy*yy+zz*zz)
          scal1 = (xx*e1x + yy*e1y + zz*e1z)/rr
          scal2 = (xx*e2x + yy*e2y + zz*e2z)/rr
          WRITE (803, '(a,i6,a,3f12.1,3i4,a,2(1pg14.6))') &
            '# Point ', i, ' at : ', &
            xx, yy, zz, nint(xx/dx), nint(yy/dy), nint(zz/dz), &
            ' projections:',scal1,scal2
        END DO
        CLOSE (803)
      END IF

      IF (jnorms /= 0) THEN
        OPEN (806, STATUS='unknown', FORM='FORMATTED', FILE='pescOrb.'//outnam)
        WRITE (806, '(a)') '# tfs, 1-|phi_1|, 1-|phi_2|, ... , 1-|phi_N|, Sum_i (1-|phi_i|)'
        CLOSE (806)
        OPEN (808, STATUS='unknown', FORM='FORMATTED', FILE='pproba.'//outnam)
        WRITE (808, '(a)') '# probability to find a certain charge state'
        WRITE (808, '(a,i6,a,3f12.1)') '#time[fs] charge_state probability'
        CLOSE (808)
      END IF

      IF (jlaser > 0) THEN
        OPEN (38, STATUS='unknown', FILE='plaser.'//outnam)
        WRITE (38, *) ' & '
        WRITE (38, *) ' tfs,power,elaser,fpulseinteg1,fpulseinteg2,fXUVinteg,'// &
          ' Ex,Ey,Ez of 1st standard pulse,'// &
          ' same for 2nd pulse, same for attopulse train'
        FLUSH  (38)
      END IF

      IF (jgeomel /= 0) THEN
        OPEN (608, STATUS='unknown', FORM='FORMATTED', FILE='pgeomel.'//outnam)
        WRITE (608, *) '# geometry parameters for electron cloud:'
        WRITE (608, *) '# time center-of-mass (3 cols) rms '// &
          ' Cartesian quadrupole tensor (3x3 cols)'
        FLUSH  (608)
      END IF

      IF (jquad /= 0) THEN
        OPEN (9, STATUS='unknown', FORM='FORMATTED', FILE='pquad.'//outnam)
        WRITE (9, '(a)') '# protocol of Cartesian quadrupole momenta, units a_0^2'
        WRITE (9, '(a)') '# time[fs] <x^2> <y^2> <z^2> <xy> <yz> <zx>'
        WRITE (9, '(a)') 'H: X YL YD YP'
        CLOSE (9)
      END IF

      IF (jinfo /= 0) THEN
        INQUIRE (17, OPENED=topenf)
        IF (.NOT. topenf) OPEN (17, POSITION='append', FILE='infosp.'//outnam)
        WRITE (17, '(a)') ' & '
        IF (.NOT. topenf) CLOSE (17)
      END IF

      IF (jspdp /= 0) THEN
        OPEN (78, STATUS='unknown', FORM='FORMATTED', FILE='pspdip.'//outnam)
        WRITE (78, '(a)') '# protocol of spin dipole moments, units a_0 '
        WRITE (78, '(a)') '# time[fs] x-component y-component z-component'
        WRITE (78, '(a)') 'H: X YL YD YP'
        CLOSE (78)
      END IF

      IF (jang /= 0) THEN
        OPEN (68, STATUS='unknown', FORM='FORMATTED', FILE='pangmo.'//outnam)
        WRITE (68, '(a)') '# protocol of orbital angular momenta, dimensionless '
        WRITE (68, '(a)') '# time[fs] L_x L_y L_z'
        WRITE (68, '(a)') 'H: X YL YD YP'
        CLOSE (68)
      END IF

      IF (jenergy /= 0) THEN
        OPEN (163, STATUS='unknown', FORM='FORMATTED', FILE='penergies.'//outnam)
        WRITE (163, *) 'col 1: time (fs)'
        WRITE (163, *) 'col 2: total energy [Ry]'
        WRITE (163, *) 'col 3: electronic binding energy [Ry]'
        WRITE (163, *) 'col 4: total s.p. energy [Ry]'
        WRITE (163, *) 'col 5: rearrangement energy [Ry]'
        WRITE (163, *) 'col 6: ionic kinetic energy [Ry]'
        WRITE (163, *) 'col 7: ionic potential energy [Ry]'
        WRITE (163, *) 'col 8: el-ion Coulomb energy [Ry]'
        WRITE (163, *) 'col 9: energy from non-local PsP [Ry]'
        WRITE (163, *) 'col 10: total energy from ionic background [Ry]'
        WRITE (163, *) 'col 11: energy absorbed from laser [Ry]'
        WRITE (163, *) 'col 12/13: internal exc. energy (spin up/down)'
        WRITE (163, *) 'col 14/15: internal exc. energy ETF (spin up/down)'
#if(raregas)
        WRITE (163, *) 'col 16: ion-ion potential energy [Ry]'
        WRITE (163, *) 'col 17: ion-surf pot'
        WRITE (163, *) 'col 18: intra-surf. pot'
        WRITE (163, *) 'col 19: kin. cores,col 20: kin. en. shells'
        WRITE (163, *) 'col 21: kin. cations'
#endif
        CLOSE (163)
      END IF

      IF (jesc /= 0) THEN
        OPEN (23, STATUS='unknown', FORM='FORMATTED', FILE='pescel.'//outnam)
        WRITE (23, '(a)') '# protocol of ionization, dimensionless'
        WRITE (23, '(a)') '# time[fs] ioniz./N_elect ionization ioniz.from absorption'
        CLOSE (23)
      END IF

      IF (jelf /= 0) THEN
        OPEN (33, STATUS='unknown', FORM='FORMATTED', FILE='pelf.'//outnam)
        WRITE (33, '(a)') '# TD electron localization FUNCTION along axes'
        CLOSE (33)
        OPEN (33, STATUS='unknown', FORM='FORMATTED', FILE='pelf2Dxy.'//outnam)
        WRITE (33, '(a)') '# TD electron localization FUNCTION in z=0 plane'
        CLOSE (33)
        OPEN (33, STATUS='unknown', FORM='FORMATTED', FILE='pelf2Dyz.'//outnam)
        WRITE (33, '(a)') '# TD electron localization FUNCTION in x=0 plane'
        CLOSE (33)
        OPEN (33, STATUS='unknown', FORM='FORMATTED', FILE='pelf2Dxz.'//outnam)
        WRITE (33, '(a)') '# TD electron localization FUNCTION in y=0 plane'
        CLOSE (33)
      END IF

      IF (jstinf /= 0) THEN
        OPEN (77, STATUS='unknown', FORM='FORMATTED', FILE='pspenergies.'//outnam)
        WRITE (77, '(a)') '# single particle energies, units Ry'
        WRITE (77, '(a)') '# time[fs] spe_1 spe_2 spe_3 ...'
        CLOSE (77)
        OPEN (76, STATUS='unknown', FORM='FORMATTED', FILE='pspvariances.'//outnam)
        WRITE (76, '(a)') '# single particle energy energy variances, units Ry'
        WRITE (76, '(a)') '# time[fs] var_1 var_2 var_3 ...'
        CLOSE (76)
#if(extended)
        OPEN (75, STATUS='unknown', FORM='FORMATTED', FILE='pspvariancesp.'//outnam)
        WRITE (75, '(a)') '# s.p. energy variances (1ph-projected)'
        CLOSE (75)
#endif
      END IF
#if(extended)
      IF (jstboostinv /= 0) THEN
        OPEN (91, STATUS='unknown', FORM='FORMATTED', FILE='pspenergybi.'//outnam)
        WRITE (91, '(a)') '# Single particle energies -- boost invariant'
        CLOSE (91)
        OPEN (92, STATUS='unknown', FORM='FORMATTED', FILE='pspvariancesbi.'//outnam)
        WRITE (92, '(a)') '# s.p. energy variances -- boost invariant'
        CLOSE (92)
      END IF
#endif

      IF (jcharges /= 0) THEN
        OPEN (323, STATUS='unknown', FILE='pcharges.'//outnam)
        WRITE (323, '(a,2i6)') '# protocol of nr. of electrons in radial shells'
        WRITE (323, '(a,2i6)') '# time [fs] nr.el.radius_1 nr.el.radius_2 ...'
        DO i = 1, INT(nzsh*dz/drcharges)
          WRITE (323, '(a,i6,f8.3)') '#Column ', i + 2, i*drcharges
        END DO
        CLOSE (323)
      END IF

      IF (jangabso /= 0) THEN
        OPEN (47, STATUS='unknown', FORM='FORMATTED', FILE='pangabso.'//outnam)
        WRITE (47, '(a)') '# protocol of 2D angular distributions in blocks for each time'
        CLOSE (47)
        OPEN (29, STATUS='unknown', FORM='FORMATTED', FILE='pangabso1D.'//outnam)
        WRITE (29, '(a)') '# protocol of phi-integrated angular distributions in blocks for each time'
        CLOSE (29)
      END IF

      IF (ionmdtyp > 0) THEN

        IF (jpos /= 0) THEN
! Positions of cluster ions
          IF (nion > 0 .AND. jgeomion > 0) THEN
            OPEN (21, STATUS='unknown', FORM='FORMATTED', FILE='pposion.'//outnam)
            WRITE (21, '(a)') '# protocol of ionic positions (in blocks)'
            WRITE (21, '(a)') '# time[fs] x y z radius'
            CLOSE (21)
            OPEN (621, STATUS='unknown', FILE='pgeomion.'//outnam)
            WRITE (621, '(a)') '# global geometry parameters for ions'
            WRITE (621, '(a)') '# time[fs] <x> <y> <z> Cartesian quadrupole moments'
            CLOSE (621)
          END IF

#if(raregas)
! Positions of GSM cores, clouds and cations
          IF (isurf /= 0) THEN
            OPEN (24, STATUS='unknown', FORM='FORMATTED', FILE='pposcore.'//outnam)
            WRITE (24, '(a)') ' & '
            CLOSE (24)
            OPEN (25, STATUS='unknown', FORM='FORMATTED', FILE='pposdip.'//outnam)
            WRITE (25, '(a)') ' & '
            CLOSE (25)
            OPEN (125, STATUS='unknown', FORM='FORMATTED', FILE='pposkat.'//outnam)
            WRITE (125, '(a)') ' & '
            CLOSE (125)
          END IF
          OPEN (424, STATUS='unknown', FORM='FORMATTED', FILE='penerSurf.'//outnam)
          WRITE (424, '(a)') ' & '
          WRITE (424, *) '# tfs,ekinC,ekinE,ekinKf, ekinSurf'
          CLOSE (424)
#endif

        END IF

! Positions of center of mass
        IF (jposcm /= 0) THEN
          OPEN (281, STATUS='unknown', FORM='FORMATTED', FILE='pposCM.'//outnam)
          WRITE (281, '(a)') '# protocol of ionic c.m. positions (in blocks)'
          WRITE (281, '(a)') '# time[fs] x y z '
          CLOSE (281)
        END IF

        WRITE(*,*) ' test JVEL,NION=',jvel,nion
        IF (jvel /= 0) THEN
          IF (nion > 0) THEN
! Velocities of molecules ions
            OPEN (22, STATUS='unknown', FORM='FORMATTED', FILE='pvelion.'//outnam)
            WRITE (22, '(a)') '# protocol  of ionic velocities (in blocks)'
            WRITE (22, '(a)') '# time[fs]  v_x  v_y  v_z  kinetic_energy'
            CLOSE (22)
            OPEN (148, STATUS='unknown', FORM='FORMATTED', &
                       FILE='pkinenion.'//outnam)
            WRITE (148, '(a)') '# protcol of ionic kinetic energies'
            WRITE (148, '(a)') '# time Ekin,ion_x  EKin,ion_y  EKin,ion_z   Ekin,ion  T_ion[Ry] T_ion[K]'
            CLOSE (148)
          END IF
#if(raregas)
          IF (nion > 0 .AND. jvel /= 0) THEN
            OPEN (151, STATUS='unknown', FORM='FORMATTED', &
                       FILE='ptempion.'//outnam)
            WRITE (151, '(a)') '# protocol of ionic kinetic energy as temperature'
            WRITE (151, '(a)') '# time[fs]   T[K]    T[Ry] '
            CLOSE (151)
          END IF
          IF (isurf /= 0) THEN
! Velocities of substrate cores
            OPEN (26, STATUS='unknown', FORM='FORMATTED', FILE='pvelcore.'//outnam)
            WRITE (26, '(a)') ' & '
            CLOSE (26)
            OPEN (126, STATUS='unknown', FORM='FORMATTED', FILE='pvelkat.'//outnam)
            WRITE (126, '(a)') ' & '
            CLOSE (126)
            OPEN (148, STATUS='unknown', FORM='FORMATTED', FILE='pkinencore.'//outnam)
            WRITE (148, '(a)') ' & '
            CLOSE (148)
            OPEN (148, STATUS='unknown', FORM='FORMATTED', FILE='pkinenkat.'//outnam)
            WRITE (148, '(a)') ' & '
            CLOSE (148)
            OPEN (148, STATUS='unknown', FORM='FORMATTED', FILE='ptempsurf.'//outnam)
            WRITE (148, '(a)') ' & '
            CLOSE (148)
! Velocities of Argon clouds (if available)
            IF (ifadiadip /= 1) THEN
              OPEN (27, STATUS='unknown', FORM='FORMATTED', FILE='pveldip.'//outnam)
              WRITE (27, '(a)') ' & '
              CLOSE (27)
            END IF
          END IF
#endif
        END IF

        IF (jforce /= 0) THEN
          IF (projcharge /= 0D0) THEN
            DO ion = 1, nion
              WRITE (ext, '(i1)') ion
              OPEN (8886, STATUS='unknown', FILE='projforce.'//ext//'.'//outnam)
              WRITE (8886, *) '" protocol of forces from projectile'
              WRITE (8886, *) '# time projectile force in x, y, z'
              FLUSH  (255)
            END DO
          END IF
          IF (jlaser > 0) THEN
            DO ion = 1, nion
              WRITE (ext, '(i1)') ion
              OPEN (255, STATUS='unknown', FILE='plforce.'//ext//'.'//outnam)
              WRITE (255, *) '" protocol of forces from projectile'
              WRITE (255, *) '# time photon force in x, y, z'
              FLUSH  (255)
            END DO
          END IF
          DO ion = 1, nion
            WRITE (ext, '(i1)') ion
            OPEN (24, STATUS='unknown', FILE='pforce.'//ext//'.'//outnam)
            WRITE (24, *) '# protocol of force on ion "i"'
            WRITE (24, *) '# time[fs] force in x, y, z'
            FLUSH  (24)
          END DO
        END IF

#if(extended)
! IF(jforce /= 0) THEN

! IF(nion > 0)THEN
! ! Forces on cluster ions
! OPEN(30,STATUS='unknown',FORM='FORMATTED', FILE='pforce.3.'//outnam)
! WRITE(30,'(a)') ' & '
! CLOSE(30)

! IF(e0 /= 0D0) THEN
! ! Forces from laser
! OPEN(31,STATUS='unknown',FORM='FORMATTED', FILE='plforce.3.'//outnam)
! WRITE(31,'(a)') ' & '
! CLOSE(31)
! END IF

! END IF
! #if(raregas)
! IF(isurf /= 0) THEN
! IF(nc+NE+nk > 0)THEN
! ! Forces on Argon cores
! OPEN(28,STATUS='unknown',FORM='FORMATTED', FILE='pforce.1.'//outnam)
! WRITE(28,'(a)') ' & '
! CLOSE(28)

! IF(e0 /= 0D0) THEN
! ! Forces from laser
! OPEN(32,STATUS='unknown',FORM='FORMATTED', FILE='plforce.1.'//outnam)
! WRITE(32,'(a)') ' & '
! CLOSE(32)
! END IF

! IF(nelect > 0)THEN
! ! Forces on Argon clouds
! ! (ELSE, there have no forces)
! OPEN(29,STATUS='unknown',FORM='FORMATTED', FILE='pforce.2.'//outnam)
! WRITE(29,'(a)') ' & '
! CLOSE(29)

! IF(e0 /= 0D0) THEN
! ! Forces from laser
! OPEN(33,STATUS='unknown',FORM='FORMATTED', FILE='plforce.2.'//outnam)
! WRITE(33,'(a)') ' & '
! CLOSE(33)
! END IF

! END IF
! END IF
! END IF
! #endif
          ! END IF
#endif
#if(raregas)
        IF (jenergion > 0) THEN

! IF (nion > 0) THEN
!! Energies of cluster
! OPEN (34, STATUS='unknown', FORM='FORMATTED', FILE='penergion.'//outnam)
! WRITE (34, '(a)') '# ionic energies:'
! WRITE (34, '(a)') '# time e_kin,x e_kin,y e_kin,z e_kin,total'
! CLOSE (34)
! END IF
          IF (nc + NE + nk > 0) THEN
!Energies of the matrix
            OPEN (35, STATUS='unknown', FORM='FORMATTED', FILE='penermat.'//outnam)
            WRITE (35, '(a)') ' & '
            CLOSE (35)
          END IF
        END IF
#endif

      END IF

    END IF ! after jump over FILE headers

    IF (nelect > 0) THEN
      IF (myn == 0) THEN
        WRITE (6, *)
        WRITE (7, *)
        WRITE (6, *) ' **** dynamical calculation - TDLDA(+MD) ****'
        WRITE (7, *) ' **** dynamical calculation - TDLDA(+MD) ****'
        WRITE (6, *)
        WRITE (7, *)
      END IF
    ELSE
      IF (myn == 0) THEN
        WRITE (6, *)
        WRITE (7, *)
        WRITE (6, *) ' **** purely ionic dynamics ****'
        WRITE (7, *) ' **** purely ionic dynamics ****'
        WRITE (6, *)
        WRITE (7, *)

#if(raregas)
        IF (isurf /= 0) THEN
          CALL adjustdip(rho)
          DO i = 1, NE
            xeinit(i) = xe(i)
            yeinit(i) = ye(i)
            zeinit(i) = ze(i)
          END DO
          CALL info(psi, rho, aloc, 0, .true.)

          IF (myn == 0) THEN
            OPEN (308, POSITION='append', FILE='energies.res')
            dist = ABS(cz(1) - maxval(zc(1:nc)))
            WRITE (308, '(1f20.10,2e25.14)') dist, energy - enerinfty, energy
            CLOSE (308)
          END IF
        END IF
#endif

      END IF
    END IF

    RETURN
  END SUBROUTINE init_dynprotocol

!-----print_densdiff--------------------------------------------------

  SUBROUTINE print_densdiff(rho, it)

! Evaluates and prints difference of actual density 'rho'
! with initial density (retrieved from file 'densdiff'.
!
! Input:
! rho = actual local density
! it = time step in calling routine

    USE params
    USE util, ONLY: pm3dcut, printfield, inttostring, tmodbuf
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: rho(2*kdfull2)
    INTEGER, INTENT(IN) :: it

    REAL(DP), DIMENSION(:), ALLOCATABLE :: w1

    INTEGER :: i
!---------------------------------------------------------------------
    IF (tmodbuf(it,jdensitydiff)) THEN

      ALLOCATE (w1(kdfull2))
      OPEN (590, STATUS='unknown', FILE='densdiff')
      DO i = 1, kdfull2
        READ (590, *) w1(i)
      END DO
      CLOSE (590)

      SELECT CASE (it)
      CASE (0:999999999)
        OPEN (689, STATUS='unknown', FILE='pdensdiff.'//trim(adjustl(inttostring(it)))//'.'//outnam)
      CASE DEFAULT
        STOP '::too many time steps::'
      END SELECT

      w1 = w1 - rho
      CALL printfield(689, w1, 'x')
      CLOSE (689)

      DEALLOCATE (w1)

    END IF

!----------------------------------------
    IF (tmodbuf(it,jdensitydiff2d)) THEN

      ALLOCATE (w1(kdfull2))
      OPEN (590, STATUS='unknown', FILE='densdiff')
      DO i = 1, kdfull2*2
        READ (590, *) w1(i)
      END DO
      CLOSE (590)

      SELECT CASE (it)
      CASE (0:999999999)
        w1 = w1 - rho
        OPEN (689, STATUS='unknown', FILE='pdensdiff2Dxy.'//trim(adjustl(inttostring(it)))//'.'//outnam)
        CALL pm3dcut(689, 1, 2, 0D0, w1)
        CLOSE (689)
        OPEN (689, STATUS='unknown', FILE='pdensdiff2Dxz.'//trim(adjustl(inttostring(it)))//'.'//outnam)
        CALL pm3dcut(689, 1, 3, 0D0, w1)
        CLOSE (689)
        OPEN (689, STATUS='unknown', FILE='pdensdiff2Dyz.'//trim(adjustl(inttostring(it)))//'.'//outnam)
        CALL pm3dcut(689, 2, 3, 0D0, w1)
        CLOSE (689)
      CASE DEFAULT
        STOP '::too many time steps::'
      END SELECT

      DEALLOCATE (w1)

    END IF

!----------------------------------------
    IF (tmodbuf(it,jdensity2d)) THEN

      SELECT CASE (it)
      CASE (0:999999999)
        OPEN (689, STATUS='unknown', FILE='pdens2Dxy.'//trim(adjustl(inttostring(it)))//'.'//outnam)
        CALL pm3dcut(689, 1, 2, 0D0, rho)
        CLOSE (689)
        OPEN (689, STATUS='unknown', FILE='pdens2Dxz.'//trim(adjustl(inttostring(it)))//'.'//outnam)
        CALL pm3dcut(689, 1, 3, 0D0, rho)
        CLOSE (689)
        OPEN (689, STATUS='unknown', FILE='pdens2Dyz.'//trim(adjustl(inttostring(it)))//'.'//outnam)
        CALL pm3dcut(689, 2, 3, 0D0, rho)
        CLOSE (689)
      CASE DEFAULT
        STOP '::too many time steps::'
      END SELECT

    END IF

    RETURN
  END SUBROUTINE print_densdiff

!-----open_protok_el----------------------------------------------

  SUBROUTINE open_protok_el(it)

! OPEN protocol files for electronic properties.
!
! Input:
! it = time step in calling routine

    USE params
    USE util, ONLY: safeopen
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: it

!----------------------------------------------------------------

    IF (jlaser > 0) CALL safeopen(38, 0, 1, 'plaser')
    CALL safeopen(163, it, jenergy, 'penergies')
    CALL safeopen(323, it, jcharges, 'pcharges')

    IF (nelect > 0 .AND. myn == 0) THEN

      CALL safeopen(8, it, jdip, 'pdip')
#if(extended)
      IF (myn == 0 .AND. jdip /= 0) CALL safeopen(88, it, jdip, 'pprojdip')
#endif
      CALL safeopen(810, it, jdiporb, 'pdiporb.x')
      CALL safeopen(811, it, jdiporb, 'pdiporb.y')
      CALL safeopen(812, it, jdiporb, 'pdiporb.z')
      CALL safeopen(9, it, jquad, 'pquad')
      CALL safeopen(17, it, jinfo, 'infosp')
      CALL safeopen(47, 0, jangabso, 'pangabso')
      CALL safeopen(29, 0, jangabso, 'pangabso1D')
      CALL safeopen(68, it, jang, 'pangmo')
      CALL safeopen(78, it, jspdp, 'pspdip')
      CALL safeopen(608, it, jgeomel, 'pgeomel')
      CALL safeopen(806, it, jnorms, 'pescOrb')
      CALL safeopen(808, it, jnorms, 'pproba')

    END IF

    RETURN
  END SUBROUTINE open_protok_el

!-----analyze_ions-----------------------------------------------

  SUBROUTINE analyze_ions(it)

! OPEN protocol files and PRINT properties of ions.
!
! Input:
! it = time step in calling routine

    USE params
!#if(extended)
!    USE util, ONLY: gettemperature, getcm, safeopen, view3d, getclustergeometry
!#else
    USE util, ONLY: gettemperature, getcm, safeopen, getclustergeometry, tmodbuf
!#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: it

    INTEGER :: ion
    REAL(DP) :: amfac, ekx, eky, ekz, r2ion, r2iona, sumx, sumy, sumz, rkt
#if(raregas)
    INTEGER :: i
    REAL(DP) :: dist2, sumnc
    REAL(DP), EXTERNAL :: getdistance2
#endif

!----------------------------------------------------------------

    tfs = it*dt1*0.0484

    CALL getcm(1, 0, 0) ! c.m. now on 'rvectmp(1:3)'
    IF (jposcm > 0 .AND. MOD(it, MAX(jposcm,1)) == 0) THEN
      CALL safeopen(281, it, jposcm, 'pposCM')
      WRITE (281, '(1f15.6,3e17.8)') tfs, rvectmp(1), rvectmp(2), rvectmp(3)
      CLOSE (281)
    END IF

    IF (((jpos > 0 .AND. MOD(it, MAX(jpos,1)) == 0) &
         .OR. (jvel > 0 .AND. MOD(it, MAX(jvel,1)) == 0) &
         .OR. (jforce /= 0 .AND. MOD(it, MAX(jforce,1)) == 0))) THEN
      CALL safeopen(21, it, jpos, 'pposion')
      CALL safeopen(621, it, jgeomion, 'pgeomion')
      CALL safeopen(22, it, jvel, 'pvelion')
      CALL safeopen(149, it, jvel, 'pkinenion')
      CALL safeopen(151, it, jvel, 'ptempion')

      r2ion = SQRT(SUM((cx(1:nion) - rvectmp(1))**2 + (cy(1:nion) - rvectmp(2))**2 &
                       + (cz(1:nion) - rvectmp(3))**2)/nion)
      sumx = 0D0
      sumy = 0D0
      sumz = 0D0
      DO ion = 1, nion
        r2iona = SQRT(cx(ion)**2 + cy(ion)**2 + cz(ion)**2)
        IF (tmodbuf(it, jpos)) THEN
          WRITE (21, '(1f13.5,3e17.8,1pg13.5)') &
            tfs, cx(ion), cy(ion), cz(ion), r2iona
          FLUSH  (21)
        END IF
        IF (tmodbuf(it, jvel)) THEN
          WRITE (22, '(1f13.5,3e17.8,1pg13.5)') &
            tfs, cpx(ion), cpy(ion), cpz(ion), ekion
          FLUSH  (22)
        END IF
        sumx = sumx + (cpx(ion)**2)/amu(np(ion))/1836D0
        sumy = sumy + cpy(ion)**2/amu(np(ion))/1836D0
        sumz = sumz + cpz(ion)**2/amu(np(ion))/1836D0
      END DO
      IF (tmodbuf(it, jvel)) THEN
        CALL gettemperature(4,rkt)
        WRITE (149, '(1f13.5,6e17.8)') &
             tfs, sumx, sumy, sumz, sumx + sumy + sumz, rkt, rkt/(6.507D-6)
        FLUSH  (149)
      END IF

      IF (tmodbuf(it,jgeomion)) THEN
        CALL getclustergeometry
        WRITE (621, '(12e15.5)') tfs, comx, comy, comz, rmsion, &
          qtion(1, 1), qtion(2, 2), qtion(3, 3), &
          qtion(1, 2), qtion(1, 3), qtion(2, 3), dmdistion
        FLUSH  (621)
      END IF

      FLUSH  (21)
      FLUSH  (22)
#if(extended)
      IF (tmodbuf(it, jpos)) CALL view3d()
#endif
    END IF


#if(raregas)
! Kinetic energy in the three directions for cluster and matrix
! and total energy (both in penergion and penermat)

!!! This piece disabled. It is a doubling of 'pkinenion.*' above.
! IF (jenergion > 0 .AND. MOD(it, jenergion) == 0 .AND. (nion) > 0) THEN
! ekx = 0D0
! eky = 0D0
! ekz = 0D0
! amfac = amu(np(nion))*1836D0*ame*2D0
! DO ion = 1, nion
! ekx = ekx + cpx(ion)*cpx(ion)
! eky = eky + cpy(ion)*cpy(ion)
! ekz = ekz + cpz(ion)*cpz(ion)
! END DO
! ekx = ekx/amfac
! eky = eky/amfac
! ekz = ekz/amfac
! CALL safeopen(34, it, jenergion, 'penergion')
! WRITE (34, '(5f13.5)') tfs, ekx, eky, ekz, etot
! FLUSH  (34)
! END IF

! ?? not clear what that is doing
    IF (jpos > 0 .AND. MOD(it, MAX(jpos,1)) == 0 .AND. nelect == 0) THEN
      sumnc = 0D0
      DO i = 1, nc
        dist2 = getdistance2(i, i + nc)
        sumnc = sumnc + 0.5D0*cspr*dist2
      END DO
      WRITE (773, '(1f12.5,1e17.7)') tfs, sumnc
    END IF
#endif

    RETURN
  END SUBROUTINE analyze_ions

!-----analyze_elect-----------------------------------------------

  SUBROUTINE analyze_elect(psi, rho, aloc, it)

! Analysis and print of electronic properties during dynamics.
!
! Input:
! psi = set of COMPLEX s.p. wavefunctions
! rho = electron density
! aloc = local mean-field potential
! it = time step in calling routine

    USE params
    USE util, ONLY: safeopen, probab, stateoverl, &
                    calcchargdist, getelectrongeometry
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
    REAL(DP), INTENT(IN) :: rho(2*kdfull2)
    REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
    INTEGER, INTENT(IN) :: it

    LOGICAL, PARAMETER :: ttest = .FALSE.
    INTEGER :: nbe

!----------------------------------------------------------------

#if(extended)
    IF (ABS(phangle) > small) CALL phoverl(psi)
#endif

    IF (jstateoverlap > 0) THEN; IF (MOD(it, jstateoverlap) == 0) THEN
      CALL stateoverl(psi, psisavex)
    END IF; END IF

! ***** dynamic density plot *****

#if(extended)
    IF (idenspl > 0 .AND. MOD(it, idenspl) == 0) CALL fmtv_fld(psi, rho, it)
#endif
! NUMBER of escaped electrons

    IF (myn == 0 .AND. jesc > 0) THEN; IF (MOD(it, jesc) == 0) THEN
      CALL nescape(it, rho)
    END IF; END IF

    IF (jcharges /= 0) THEN; IF (MOD(it, jcharges) == 0) THEN
      CALL calcchargdist(rho)
    END IF; END IF

! Time-dependent Electron Localization function (TD-ELF)

    IF (jelf > 0) THEN; IF (MOD(it, jelf) == 0) THEN
#if(!mpi)
      CALL localize(rho, psi, it)
#else
      STOP ' LOCALIZE (switch JELF) should not be invoked in parallele code'
#endif
    END IF; END IF

! everything about single particles

    IF ((jstinf > 0 .AND. MOD(it, MAX(jstinf,1)) == 0) &
        .OR. (jinfo > 0 .AND. MOD(it, MAX(jinfo,1)) == 0) &
        .OR. (jenergy > 0 .AND. MOD(it, MAX(jenergy,1)) == 0) &
        .OR. (jesc > 0 .AND. jnorms > 0 .AND. MOD(it, MAX(jnorms,1)) == 0) &
        .OR. (jdiporb > 0 .AND. MOD(it, MAX(jdiporb,1)) == 0)) THEN
#if(mpi)
      IF (ttest) WRITE (*, *) ' ANALYZE before INFO: myn=', myn
#endif
#if(nompi)
      IF (ttest) WRITE (6, '(a,i4)') ' INFO from ANALYZE. it=', it
#endif
      CALL info(psi, rho, aloc, it, .true.)
#if(mpi)
      IF (ttest) WRITE (*, *) ' ANALYZE after INFO: myn=', myn
#endif
    END IF

! escaped electrons analyzed orbital by orbital

    IF (jnorms > 0) THEN; IF (MOD(it, jnorms) == 0) THEN
#if(!mpi)
      CALL safeopen(806, it, jnorms, 'pescOrb')
      WRITE (806, '(500(1pg13.5))') tfs, 1D0 - spnorm(1:nstate), &
        SUM(1D0 - spnorm(1:nstate))/nstate

      FLUSH  (806)
#endif

      CALL probab(psi)
    END IF; END IF

    IF (jang > 0) CALL instit(psi)

    IF (myn == 0) THEN
    IF (MOD(it, 100) == 0) THEN
      WRITE (6, '(a,i6,a,f12.6,a,f16.5)') &
        'iter=', it, ' tfs=', tfs, ' total energy=', etot
    ELSE IF (nelect > 0) THEN
      WRITE (6, '(a,i6,a,f12.6,a,f16.10)') &
        'iter=', it, ' tfs=', tfs, ' total energy=', etot
    END IF

      IF (jdip > 0) THEN; IF (MOD(it, jdip) == 0) THEN
        WRITE (8, '(f10.5,3e17.8)') tfs, qe(2), qe(3), qe(4)
        FLUSH  (8)
      END IF; END IF

#if(extended)
      IF (jdip > 0 .AND. MOD(it, jdip) == 0 .AND. eproj /= 0) THEN
        WRITE (88, '(f10.5,6e17.8)') tfs, qeproj(2), qeproj(3), qeproj(4), &
          qetarget(2), qetarget(3), qetarget(4)
        FLUSH  (88)
      END IF
#endif

      IF (jdiporb > 0) THEN; IF (MOD(it, jdiporb) == 0) THEN
        WRITE (810, '(f10.5,1000e17.8)') &
          tfs, (qeorb_all(nbe, 3), nbe=1, nstate_all)
        WRITE (811, '(f10.5,1000e17.8)') &
          tfs, (qeorb_all(nbe, 4), nbe=1, nstate_all)
        WRITE (812, '(f10.5,1000e17.8)') &
          tfs, (qeorb_all(nbe, 5), nbe=1, nstate_all)
        FLUSH  (810)
        FLUSH  (811)
        FLUSH  (812)
      END IF; END IF

      IF (jquad > 0) THEN; IF (MOD(it, jquad) == 0) THEN
        WRITE (9, '(f9.4,6f11.4)') tfs, qe(5), qe(6), qe(7), qe(8), qe(9), qe(10)
        FLUSH  (9)
      END IF; END IF
    END IF

    IF (jmp > 0 .AND. MOD(it, MAX(jmp,1)) == 0) CALL evalmp(803, psi)

    IF (myn == 0 .AND. it > 0 .AND. jangabso > 0 .AND. &
        MOD(it, MAX(jangabso,1)) == 0) CALL angular_distribution()

    IF (jgeomel > 0) THEN; IF (MOD(it, jgeomel) == 0) THEN
      CALL getelectrongeometry(psi, 0)
      WRITE (608, '(11e15.5)') tfs, codx, cody, codz, rmsel, &
        qtel(1, 1), qtel(2, 2), qtel(3, 3), qtel(1, 2), qtel(1, 3), qtel(2, 3)
    END IF; END IF

    IF (jang > 0) THEN; IF (MOD(it, jang) == 0) THEN
      WRITE (68, '(f12.5,3f11.5)') tfs, ajx, ajy, ajz
      FLUSH  (68)
    END IF; END IF

    IF (myn == 0) THEN

      IF (jspdp > 0) THEN; IF (MOD(it, jspdp) == 0) THEN
        WRITE (78, '(f10.5,3f13.5)') tfs, se(1), se(2), se(3)
        FLUSH  (78)
      END IF; END IF
  
      WRITE (7, '(a,f8.4,a,f7.2,a,3f9.2)') &
        't=', tfs, ' moments: monop.=', qe(1), ' dip.: ', qe(2), qe(3), qe(4)
      WRITE (7, '(a,i5,a/)') '** END of iter= ', it, '**'

      WRITE (6, '(a,f8.4,a,f7.2,a,3f9.2)') &
        't=', tfs, ' moments: monop.=', qe(1), ' dip.: ', qe(2), qe(3), qe(4)

    END IF

    RETURN
  END SUBROUTINE analyze_elect

!-----tstep_exp-------------------------------------------------------------

  SUBROUTINE tstep_exp(q0, aloc, rho, it, qwork, timagtime)

! One electronic time step by exponential evolution, consisting
! of a half step to produce an intermediate mean field followed
! by a full with using this mean field.
!
! q0 = s.p. wavefunctions to be propagated
! aloc = array for local mean field
! rho = array for local density
! it = nr. of time step
! qwork = work space for s.p. wavefunctions
! timagtime = LOGICAL variable, switch to imaginary time step

    USE params

#if(fsic)
    USE twost, ONLY: tnearest, eval_unitrot
#endif

    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN OUT) :: q0(kdfull2, kstate)
    REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
    REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
    INTEGER, INTENT(IN) :: it
    COMPLEX(DP), INTENT(OUT) :: qwork(kdfull2, kstate)
    LOGICAL, INTENT(IN) :: timagtime

    INTEGER :: nb, nterms, ithr
    COMPLEX(DP), ALLOCATABLE :: q1(:, :)
    COMPLEX(DP) :: cdtact
#if(extended)
!COMPLEX(DP),DIMENSION(:,:,:),POINTER :: q3D
#endif

#if(fsic)
! The parameter 'tnorotate' de-activates the subtraction of the
! Lagrangian matrix in the SIC step. The version of exponential
! evolution with subtraction of the Lagrangian matrix is found
! in 'exp_evolp'. The strategy needs yet testing. It was not
! really beneficial so far. (01/2013)
!
! what is the current STATUS of this issue ? F.L. (03/2017)
    LOGICAL, PARAMETER :: tnorotate = .true.
#else
    LOGICAL, PARAMETER :: tnorotate = .true. ! switch for full SIC (unsued here)
#endif

    LOGICAL, PARAMETER :: ttimestop = .TRUE.
    REAL(DP) :: t1, t2
    INTEGER :: it1, it2, it3, it4, it5
! acctime is accumulator for time usage
! acctime(1) = half step
! acctime(2) = mean field half
! acctime(3) = full step
! acctime(4) = mean field full
! acctime(6) = TSTEP_EXP
! acctime(7) = nr. of calls to TSTEP
    REAL(DP), SAVE :: acctime(1:10) = &
      (/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/)

    EXTERNAL :: exp_evol

!----------------------------------------------------------------------

#if(mpi)
    STOP 'exponential evolution not yet MPI parallelized'
#else
    myn = 0
#endif


#if(raregas)
    IF (nc + NE + nk > 0) STOP 'TSTEP_EXP not appropriate for rare gas'
#endif

    ALLOCATE (q1(kdfull2, 0:nthr))
#if(extended)
!#if(fftw_cpu)
! q3D(1:kxbox,1:kybox,1:kzbox) => q0(1:kdfull2,nb)
!#endif
#endif

! one half time step to define new mean field
! use exponential evolution to second order

    IF (ifsicp == 5) psisavex = q0

    IF (ttimestop) THEN
      acctime(1:5) = 0D0
      CALL system_clock(it1)
    END IF

    IF (.NOT. timagtime) THEN
      cdtact = CMPLX(dt1/2D0, 0D0, DP)
      IF (tnorotate .OR. ifsicp < 8) THEN
#if(dynomp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nb,ithr) SCHEDULE(STATIC)
#endif
        DO nb = 1, nstate
          qwork(:, nb) = q0(:, nb)
#if(dynomp)
          ithr = OMP_GET_THREAD_NUM()
#else
          ithr = 0
#endif
          CALL exp_evol(qwork(:, nb), aloc, nb, 4, cdtact, q1(:, ithr))
        END DO
#if(dynomp)
!$OMP END PARALLEL DO
#endif
#if(fsic)
      ELSE
        qwork = q0
        CALL exp_evolp(qwork, aloc, 4, cdtact, q1, q0)
#endif
      END IF
      IF (ttimestop) CALL system_clock(it2)

#if(fsic)
      IF (tnearest .AND. ifsicp .GE. 8) CALL eval_unitrot(qwork, q0)
#endif

! compute mean field at half time step
      CALL dyn_mfield(rho, aloc, qwork, dt1*0.5D0, it)

    END IF
    IF (ttimestop) CALL system_clock(it3)

! full time step to next wavefunctions
! use exponential evolution to fourth order

    IF (timagtime) THEN
      cdtact = CMPLX(0D0, -dt1, DP)
    ELSE
      cdtact = CMPLX(dt1, 0D0, DP)
    END IF

    nterms = 8 ! 4

    IF (tnorotate .OR. ifsicp < 8) THEN
#if(dynomp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(nb,ithr) SCHEDULE(STATIC)
#endif
      DO nb = 1, nstate
#if(dynomp)
        ithr = OMP_GET_THREAD_NUM()
#else
        ithr = 0
#endif
        CALL exp_evol(q0(:, nb), aloc, nb, nterms, cdtact, q1(:, ithr))
      END DO
#if(dynomp)
!$OMP END PARALLEL DO
#endif
#if(fsic)
    ELSE
      qwork = q0
      CALL exp_evolp(q0, aloc, nterms, cdtact, q1, qwork)
#endif
    END IF
    IF (ttimestop) CALL system_clock(it4)

#if(fsic)
    IF (tnearest .AND. ifsicp .GE. 8 .AND. .NOT. timagtime) CALL eval_unitrot(q0, qwork)
#endif

#if(extended)
! switch off absorption after times step 'ntref'
    IF (ntref > 0 .AND. it > ntref) nabsorb = 0
#endif

    DEALLOCATE (q1)

! compute mean field at new time

    CALL dyn_mfield(rho, aloc, q0, dt1, it)
    IF (ttimestop) THEN
      CALL system_clock(it5)
      acctime(1) = acctime(1) + (it2 - it1)*systime_factor
      acctime(2) = acctime(2) + (it3 - it2)*systime_factor
      acctime(3) = acctime(3) + (it4 - it3)*systime_factor
      acctime(4) = acctime(4) + (it5 - it4)*systime_factor
      acctime(6) = acctime(6) + SUM(acctime(1:4))
      acctime(7) = acctime(7) + 1D0
      WRITE (6, '(a,4(1pg13.5))') 'TIMES sub-steps=', acctime(1:4)
      WRITE (6, '(a,1pg13.5)') 'CPU in TSTEP_EXP=', acctime(6)/acctime(7)
    END IF

    RETURN
  END SUBROUTINE tstep_exp

!-----exp_evol-------------------------------------------------------------

  SUBROUTINE exp_evol(qact, aloc, nbe, norder, dtact, qwork)

! Propagation of the wavefunction of state 'nbe' by Taylor
! expanded evolution:
! qact = wavefunction on which H acts and resulting w.f.
! aloc = local potential for the actual spin component
! ak = kinetic energies in momentum space
! nbe = NUMBER of state
! norder = order of expansion (4 recommended for full step))
! dtact = time step

! Note: The propagation uses the action of the Hamiltonian
! where the diagonal element (s.p.energy) is subtracted.
! That diagonal element is evaluated in the first order
! CALL 'nterm=1'.

    USE params
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN OUT) :: qact(kdfull2)
    REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
    INTEGER, INTENT(IN) :: nbe
    INTEGER, INTENT(IN) :: norder
    COMPLEX(DP), INTENT(IN) :: dtact
    COMPLEX(DP), INTENT(OUT) :: qwork(kdfull2)

    COMPLEX(DP) :: dti, cfac
    INTEGER :: i, ilocbas, isig, nterm

!----------------------------------------------------------------------

    IF (ispin(nrel2abs(nbe)) == 1) THEN
      ilocbas = 1
    ELSE IF (ispin(nrel2abs(nbe)) == 2) THEN
      ilocbas = nxyz + 1
    ELSE
      STOP " EXPEVOL: spin index must be 1 or 2"
    END IF
    IF (ABS(AIMAG(dtact)) > 1D-10) THEN
      isig = -1
    ELSE
      isig = 1
    END IF

    dti = dtact*CMPLX(0D0, 1D0, DP)
    cfac = CMPLX(1D0, 0D0, DP)
    DO i = 1, nxyz
      qwork(i) = qact(i)
    END DO
    DO nterm = 1, norder
      CALL hpsi(qwork, aloc(ilocbas), nbe, isig*nterm)
      cfac = -dti/nterm*cfac
      DO i = 1, nxyz
        qact(i) = qact(i) + cfac*qwork(i)
      END DO
    END DO
    RETURN
  END SUBROUTINE exp_evol

#if(fsic)
!-----exp_evolp-------------------------------------------------------------

  SUBROUTINE exp_evolp(qact, aloc, norder, dtact, qwork, psi)

! Propagation of all s.p. wavefuntions by Taylor expanded
! exponential evolution (version needed for SIC):
! qact = wavefunction on which H acts and resulting w.f.
! aloc = local potential for the actual spin component
! ak = kinetic energies in momentum space
! norder = order of epxansion (4 recommended for full step))
! dtact = time step
! qwork = work space for COMPLEX wavefunction
! psi = set of wavefunctions before the step

! Note: The propagation uses the action of the Hamiltonian
! WHERE the diagonal element (s.p.energy) is subtracted.
! That diagonal element is evalutaed in the first order
! CALL 'nterm=1'.

    USE params
    USE util, ONLY: wfovlp
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN OUT) :: qact(kdfull2, kstate)
    REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
    INTEGER, INTENT(IN) :: norder
    COMPLEX(DP), INTENT(IN) :: dtact
    COMPLEX(DP), INTENT(OUT) :: qwork(kdfull2)
    COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)

    COMPLEX(DP), ALLOCATABLE :: chmatrix(:, :)

    COMPLEX(DP) :: dti, cfac, cacc(kstate)
    INTEGER :: i, ilocbas
    INTEGER :: na, nbe, ncs, nterm

!----------------------------------------------------------------------

    ALLOCATE (chmatrix(kstate, kstate))

    dti = dtact*CMPLX(0D0, 1D0, DP)

! compute H-matrix, store h*psi wavefunctions

    DO nbe = 1, nstate
      ilocbas = 1 + (ispin(nrel2abs(nbe)) - 1)*nxyz
      CALL hpsi(qact(1, nbe), aloc(ilocbas), nbe, 1)
      DO ncs = 1, nstate
        IF (ispin(nrel2abs(nbe)) == ispin(nrel2abs(ncs))) THEN
          chmatrix(ncs, nbe) = wfovlp(psi(:, ncs), qact(:, nbe))
        ELSE
          chmatrix(ncs, nbe) = CMPLX(0D0, 0D0, DP)
        END IF
      END DO
    END DO
! symmetrize H-matrix
    DO nbe = 1, nstate; DO ncs = 1, nbe - 1
        chmatrix(ncs, nbe) = (chmatrix(ncs, nbe) + CONJG(chmatrix(nbe, ncs)))/2D0
        chmatrix(nbe, ncs) = CONJG(chmatrix(ncs, nbe))
      END DO; END DO

! now the Taylor expansion (recycle stored h*psi in first step)

    DO nbe = 1, nstate
      ilocbas = 1 + (ispin(nrel2abs(nbe)) - 1)*nxyz
      cfac = CMPLX(1D0, 0D0, DP)
      DO nterm = 1, norder
        ! H*psi and prepare subtraction factors
        IF (nterm == 1) THEN
          qwork(:) = qact(:, nbe)
          qact(:, nbe) = psi(:, nbe) ! restore original wavefunctions
          cacc(:) = chmatrix(:, nbe)
        ELSE
          DO ncs = 1, nstate
            IF (ispin(nrel2abs(nbe)) == ispin(nrel2abs(ncs))) THEN
              cacc(ncs) = CMPLX(0D0, 0D0, DP)
              DO na = 1, nstate
                IF (ispin(nrel2abs(na)) == ispin(nrel2abs(ncs))) THEN
                  cacc(ncs) = cacc(ncs) &
                     + chmatrix(ncs, na)*wfovlp(psi(:, na), qwork)
                END IF
              END DO
            END IF
          END DO
          CALL hpsi(qwork, aloc(ilocbas), nbe, 2)
        END IF
        ! project H-matrix
        DO ncs = 1, nstate
          IF (ispin(nrel2abs(nbe)) == ispin(nrel2abs(ncs))) THEN
            qwork(:) = qwork(:) - psi(:, ncs)*cacc(ncs)
          END IF
        END DO
        ! accumulate to Taylor series
        cfac = -dti/nterm*cfac
        qact(:, nbe) = qact(:, nbe) + cfac*qwork(:)
      END DO
    END DO

    DEALLOCATE (chmatrix)

    RETURN
  END SUBROUTINE exp_evolp
#endif

!-----hpsi -------------------------------------------------------------

  SUBROUTINE hpsi(qact, aloc, nbe, itpri)

! ACTION of mean-field Hamiltonian on one s.p. wavefunction:
! qact = wavefunction on which H acts and resulting w.f.
! aloc = local potential for the actual spin component
! ak = kinetic energies in momentum space
! nbe = number of state
! itpri = switch for computing s.p. energies (for ABS(itpri)=1)
!         <0 switches to subtract mean-value of s.p. energy

    USE params
    USE util, ONLY: wfovlp
    USE kinetic
#if(fsic)
    USE twost
#endif

    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN OUT) :: qact(kdfull2)
    REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
!REAL(DP), INTENT(IN) :: akv(kdfull2)
    INTEGER, INTENT(IN) :: nbe
    INTEGER, INTENT(IN) :: itpri
    COMPLEX(DP), ALLOCATABLE :: qex(:)

! workspaces
    COMPLEX(DP), ALLOCATABLE :: q1(:), q2(:), q1fine(:), q2fine(:), qactfine(:)
    COMPLEX(DP), ALLOCATABLE :: qarray(:, :, :), qarrayfine(:, :, :)
    REAL(DP) :: wnorm
    LOGICAL :: tpri
    LOGICAL, PARAMETER :: tsubmean = .TRUE.
    LOGICAL, PARAMETER :: ttest = .FALSE.
    INTEGER :: i, is, na
    COMPLEX(DP) :: cf

!----------------------------------------------------------------------

    tpri = ABS(itpri) == 1

    ALLOCATE (q1(kdfull2), q2(kdfull2))
    q1 = (0D0, 0D0)
    q2 = (0D0, 0D0)
! ACTION of kinetic energy

#if(fftw_cpu)
    CALL kinetic_energy(qact, q2)
#endif
#if(netlib_fft)
    CALL fftf(qact, q1)
    DO i = 1, nxyz
      q1(i) = akv(i)*q1(i)
    END DO
    CALL fftback(q1, q2)
#endif
#if(findiff|numerov)
    CALL ckin3d(qact, q1)
#endif

! action of potential and non-local PsP (optionally)

    IF (ipsptyp == 1) THEN
      CALL nonlocalc(qact, q1, 0)

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
      CALL exchg(qact, qex, nbe)
      q1 = q1 + qex
      DEALLOCATE (qex)
    END IF

#if(fsic)
! subtract SIC potential for state NBE
    IF (ifsicp .GE. 8) THEN
      is = ispin(nrel2abs(nbe))
      DO na = 1, nstate
        IF (ispin(nrel2abs(na)) == is) THEN
          cf = wfovlp(psiut(:, na), qact)
          DO i = 1, nxyz
            q1(i) = q1(i) - qnewut(i, na)*cf
          END DO
        END IF
      END DO
    END IF
#endif

    IF (tpri) THEN
      epotsp(nbe) = wfovlp(qact, q1)
      amoy(nbe) = ekinsp(nbe) + epotsp(nbe)
      q2 = q1 + q2
      spvariance(nbe) = SQRT(MAX(REAL(wfovlp(q2, q2), DP) &
                                 - ABS(wfovlp(qact, q2))**2, 1D-99))
      is = ispin(nrel2abs(nbe))
      IF (ttest) WRITE (*, '(a,2i4,5(1pg13.5))') &
        ' HPSI: nbe,is,esp,var=', nbe, is, amoy(nbe), spvariance(nbe), &
        ekinsp(nbe), epotsp(nbe), REAL(wfovlp(q2, q2), DP)
      FLUSH  (6)
    ELSE
      q2 = q1 + q2
    END IF

    IF (itpri < 0) THEN
      qact = q2 - amoy(nbe)*qact
    ELSE
      qact = q2
    END IF

    DEALLOCATE (q1, q2)

    RETURN
  END SUBROUTINE hpsi

