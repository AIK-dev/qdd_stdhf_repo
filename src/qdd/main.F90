PROGRAM tdlda_m

! Main PROGRAM of code package for "electronic dynamics in clusters"
! according to TDLDA-MD with various options for a SIC.

! *******************************************
! declarations
! *******************************************

  USE params
  USE kinetic
  USE util
  USE coulsolv
  USE orthmat
  USE RTA_module

#if(fsic)
  USE twostr
  USE twost
  USE twost_util
#endif

#ifdef QDD_INFO
  USE QDDInfo, ONLY: writeQDDVersionInfo
#endif /* QDD_INFO */

  IMPLICIT NONE

! psir = REAL wavefunctions in coord space
! psi = COMPLEX wavefunctions in coord space (for dynamics)
! rho = electronic density in coord space
! aloc = mean-field-potential (to be initialized before dynamics)
! chpcoul = Coulomb-potential of electrons
! c$s/cp$s = auxiliary field to store ionic coords and momenta in between
! the trial ionic propagation
! rhojel = jellium density

#if(mpi)
  INCLUDE 'mpif.h'
  INTEGER :: is(mpi_status_size)
#endif

  REAL(DP), ALLOCATABLE :: aloc(:), rho(:)
  REAL(DP), ALLOCATABLE :: psir(:, :)
  COMPLEX(DP), ALLOCATABLE, TARGET :: psi(:, :)
  COMPLEX(DP), ALLOCATABLE :: psiw(:, :)

  INTEGER :: ion, it, nspup
  INTEGER :: i
#if(raregas)
  REAL(DP):: dt
#endif
  REAL(DP):: time_absfin

  LOGICAL, PARAMETER :: imaginary_time = .true. ! activate static afterburn

  INTERFACE
    SUBROUTINE dyn_propag(psi, rho, aloc)
      USE params
      COMPLEX(DP), INTENT(IN OUT), TARGET :: psi(kdfull2, kstate)
      REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
      REAL(DP), INTENT(IN OUT) :: rho(2*kdfull2)
    END SUBROUTINE
  END INTERFACE

! *******************************************
!
! initializations
!
! *******************************************

#ifdef QDD_INFO
  CALL writeQDDVersionInfo()
#endif

  CALL cpu_time(time_absinit)

  CALL init_parallele()

  CALL initnamelists ! READ all input parameters

! Update Dynamic thread adjustment and global no. of threads IF set in input FILE
#if(omp)
!  CALL OMP_SET_NUM_THREADS(numthr) ! to also set the system wide NUMBER of OMP threads
!  CALL OMP_SET_DYNAMIC(setdyn)
#if(omp_debug)
  WRITE (*, *) "setdyn = ", setdyn
  WRITE (*, *) "numthr = ", numthr
  WRITE (*, *) "OMP_GET_MAX_THREADS() = ", OMP_GET_MAX_THREADS()
  WRITE (*, *) "OMP_GET_NUM_PROCS() = ", OMP_GET_NUM_PROCS()
  WRITE (*, *) "OMP_GET_DYNAMIC() = ", OMP_GET_DYNAMIC()
  WRITE (*, *) "OMP_GET_NESTED() = ", OMP_GET_NESTED()
!$OMP PARALLEL
  WRITE (*, *) "OMP_GET_NUM_THREADS() = ", OMP_GET_NUM_THREADS()
  WRITE (*, *) "OMP_GET_THREAD_NUM() = ", OMP_GET_THREAD_NUM()
!$OMP END PARALLEL
#endif
  nthr = numthr - 1
#endif

  IF (icooltyp == 3) CALL init_simann() ! initialize simulated annealing

  CALL init_baseparams() !init grid size, number of states ...

#if(raregas)
  CALL check_isrtyp ! check short range part of substrate
#endif

  CALL iperio ! initializing the 'periodic table'

  CALL changeperio ! pseudo-potential parameters

  CALL iparams() ! check dynamic parameters

  CALL checkoptions() !check coherence of preprocessor option

  CALL init_grid() ! init coulomb solver, kinetic energy, grid properties

  CALL init_fields() ! allocate basic arrays in module params

  ALLOCATE (psir(kdfull2, kstate))
  psir = 0D0
  ALLOCATE (aloc(2*kdfull2), rho(2*kdfull2))
  aloc = 0D0
  rho = 0D0

  IF (myn == 0) CALL ocoption(7) ! output compiled options
  IF (myn == 0) CALL ocoption(8) ! output compiled options

  CALL init_output() ! headers for basic output files

! Choice of ionic background :
  WRITE (*, *) ' nion2=', nion2
  SELECT CASE (nion2)
  CASE (0)
    CALL init_jellium() ! initialize jellium background
  CASE (1)
    WRITE (*, *) ' ions switch'
    CALL initions() ! reading ionic positions and inits (for005ion.<NAME>)
  CASE (2)
    WRITE (*, *) ' EXTERNAL background potential '
    CALL pseudo_external() ! read exernal ionic pseudopotential from a file
  CASE DEFAULT
    STOP 'nion2 has incorrect value. nion2 must be 0, 1 or 2'
  END SELECT

  CALL initwf(psir) ! init wf, jellium, static parameters

#if(raregas)
  ! initialize substrate
  IF (isurf == 1) THEN
    CALL initfunctions
    CALL initsurface
  END IF
#endif

#if(mpi)
  CALL init_boxpara()
  WRITE (*, *) 'lengnod:', lengnod
#endif

#if(extended)
!IF(ihome == 1) CALL init_homfield() ! OPTIONAL homogeneous exernal E field
#endif

  CALL timer(1) ! set timer

#if(fsic)
  IF (nelect > 0 .AND. ifsicp > 7) THEN
    CALL init_fsicr() ! initialize parameters for full SIC
  END IF
#endif

! *******************************************
!
! static
!
! *******************************************

  IF (nelect > 0 .AND. irest == 0 .AND. ismax > 0) THEN
    CALL statit(psir, rho, aloc)
  END IF

! switch to Monte-Carlo optimization of ionic configurations

  IF (icooltyp == 3) THEN
    CALL simann(psir, rho, aloc)
    STOP ' Monte Carlo completed'
  END IF

  DEALLOCATE (psir)

#if(extended)
! imaginary-time iteration (to improve statics)
  IF (imaginary_time .AND. isitmax > 0 .AND. ismax > 0) THEN
    ifexpevol = .TRUE.
    CALL afterburn(psir, rho, aloc)
  END IF
#endif

! *******************************************
!
! dynamic
!
! *******************************************

  CALL init_dynamic

  CALL dyn_propag(psi, rho, aloc)

#if(raregas)
  IF (iuselast == -3) THEN
    OPEN (89, STATUS='unknown', FILE='for005surf.cool')
    WRITE (89, '(a,i6,a,i6)') '# nc = ', nc, ' nk= ', nk
    DO i = 1, nc
      WRITE (89, '(6e20.9,2i6)') xc(i), yc(i), zc(i), xe(i), ye(i), ze(i), imobc(i), &
        imobe(i)
    END DO
    DO i = 1, nk
      WRITE (89, '(3e20.9,i6)') xk(i), yk(i), zk(i), imobk(i)
    END DO
    CLOSE (89)
  END IF
#endif

! epilog

  IF (jcharges /= 0) CLOSE (323)

#if(raregas)
  IF (isurf == 1) THEN
    CLOSE (121)
    CLOSE (122)
    CLOSE (123)
    CLOSE (124)
    CLOSE (126)
    CLOSE (24)
    CLOSE (128)
    CLOSE (130)
    CLOSE (157)
  END IF
#endif

#if(extended)
  IF (jattach > 0) THEN
    DEALLOCATE (psi_target)
    DEALLOCATE (ispin_target)
    DEALLOCATE (occ_target)
    DEALLOCATE (spe_target)
    DEALLOCATE (match)
  END IF
  IF (nproj_states > 0) DEALLOCATE (proj_states)
#endif

  CLOSE (163)
  CLOSE (68)
  CLOSE (8)
  CLOSE (9)
  CLOSE (78)
  IF (jnorms > 0) THEN
    CLOSE (806)
    CLOSE (808)
  END IF

! ******************** END of main PROGRAM ****************************

  DEALLOCATE (psi)
!  DEALLOCATE (psitophi)

#if(mpi)
  CALL mpi_finalize(mpi_ierror)
#endif


CONTAINS

  SUBROUTINE init_dynamic()
!------------------------------------------------------------

    IMPLICIT NONE

    REAL(DP), EXTERNAL:: energ_ions ! declared in ion_md
    INTEGER, PARAMETER :: itest = 0
    INTEGER :: i, ion, nspup

! optionally initialize wavefunction and work arrays
    ALLOCATE (psi(kdfull2, kstate))
    psi = CMPLX(0D0, 0D0, DP)

#if(extended)
    IF (nabsorb > 0 .AND. jescmaskorb /= 0) &
         ALLOCATE (rhoabsoorb(kdfull2, kstate))
#endif

    IF (ifexpevol) ALLOCATE (psiw(kdfull2, kstate))

! initialize protocol files
    IF (nelect > 0 .AND. nabsorb > 0) CALL init_absbc(rho)
    IF (nelect > 0 .AND. jmp > 0) CALL initmeasurepoints
    IF (myn == 0 .OR. knode == 1) THEN
WRITE(*,*) 'test MYN=',myn
#if(raregas)
      CALL init_dynprotocol(rho, aloc, psi)
#else
      CALL init_dynprotocol()
#endif
      IF (nelect == 0) CALL getforces(rho, psi, -1, 0) ! initial call if MD only
    END IF

#if(raregas)
    IF (surftemp > 0) CALL init_surftemp()
#endif

    IF (nelect > 0) THEN

#if(fsic)
      IF (ifsicp >= 7) THEN
        CALL init_fsic()
      END IF
#endif

      IF (nabsorb > 0) CALL init_absbc(rho)
      IF (jmp > 0) CALL initmeasurepoints

! *** HOW TO USE THE RESTART ***
! 1. to start a dynamic run using saved REAL wave functions:
! set irest=0 and tstat=.TRUE. in the dynamic part of for005.outnam
! 2. to resume a dynamic run, set irest=1 in dynamic part
! of for005.outnam
! 3. save data every n iterations by setting isave=n
! ******************************

! initialize or recover wavefunctions
      IF (irest == 0) THEN
        CALL init_dynwf(psi)
        IF (nabsorb > 0) CALL init_abs_accum()
        WRITE (*, *) ' INIT_DYNWF branch' !test
      ELSE
#if(extended)
        IF (iscatterelectron /= 0) CALL init_scattel(psi)
#endif
        CALL restart2(psi, outnam, .false.)
        WRITE (7, '(a,i3)') 'restart irest=', irest
      END IF
      IF(nelect>3) THEN
        WRITE (*, *) 'INIT_DYN: psi', SUM(ABS(psi(1:kdfull2, 1))**2) !test
        WRITE (*, *) 'INIT_DYN: psi', SUM(ABS(psi(1:kdfull2, 2))**2) !test
        WRITE (*, *) 'INIT_DYN: psi', SUM(ABS(psi(1:kdfull2, 3))**2) !test
        WRITE (*, *) 'INIT_DYN: psi', SUM(ABS(psi(1:kdfull2, 4))**2) !test
      END IF
  WRITE(*,'(a,200(1pg13.5))') 'OCCUP after init_dynwf:',occup(1:nstate)

      ! order states in contingent blocks of spin
!      ALLOCATE (psitophi(nstate, nstate))
!      psitophi = 0D0
!      DO i = 1, nstate
!        psitophi(i, i) = 1D0
!      END DO
      IF (irest == 0) THEN !only for first pass. if irest <>0 no temperature
        CALL calcpseudo()
        CALL dyn_mfield(rho, aloc, psi, 0D0, 0)
        CALL info(psi, rho, aloc, 0, .true.) ! update occupation numbers and amoy
!
! in case of RTA, set an initial distribution with temperature 'rtatempinit'
! to give the RTA procedure space for energies below the initial energy 
!
        IF (rtatempinit > 1D-3) THEN
#if(extended)
          IF(nhstate==0) THEN
            CALL fermi_init(amoy, rtatempinit, occup, 1)
            CALL fermi_init(amoy, rtatempinit, occup, 2)
            WRITE(*,*) 'running Fermi - extended branch',nhstate
          END IF
          WRITE(*,*) 'after Fermi - extended branch'
#else
          CALL fermi_init(amoy, rtatempinit, occup, 1)
          CALL fermi_init(amoy, rtatempinit, occup, 2)
#endif
          CALL dyn_mfield(rho, aloc, psi, 0D0, 0)
        END IF
        CALL info(psi, rho, aloc, 0, .true.) !update occupation numbers and amoy
      END IF
  WRITE(*,'(a,200(1pg13.5))') 'OCCUP after RTA init:',occup(1:nstate)
      nspup = 0
      DO i = 1, nstate
        IF (ispin(i) == 1) nspup = nspup + 1 !count nspup
      END DO
      DO i = 1, nstate
        IF (i <= nspup) THEN
          ispin(i) = 1
        ELSE
          ispin(i) = 2
        END IF
      END DO
      eqstnspup = nspup
      eqstnspdw = nstate - eqstnspup


! optinally refresh (pseudo)potentials
      IF (irest > 0 .OR. tstat) THEN
        CALL calcpseudo()
        CALL calclocal(rho, aloc) ! ??
        IF (ifsicp > 0) CALL calc_sic(rho, aloc, psi)
        IF (ipsptyp == 1) THEN
        DO ion = 1, nion
#if(extended)
        IF (tinterpol) THEN
          CALL calc_projFine(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
        ELSE
          CALL calc_proj(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
        END IF
#else
        CALL calc_proj(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
#endif
        END DO
        END IF
      END IF
      CALL dyn_mfield(rho, aloc, psi, 0D0, 0)

! intialize protocols
      IF (irest == 0) CALL info(psi, rho, aloc, 0, .true.)
      IF (jang > 0) CALL instit(psi) ! notes
      IF (jdensity1d > 0) THEN
        CALL rhointxy(rho, 0)
        CALL rhointxz(rho, 0)
        CALL rhointyz(rho, 0)
      END IF

      time = irest*dt1*0.0484D0/(2D0*ame)
      IF (myn == 0 .OR. knode == 1) THEN
        WRITE (7, '(a,i6,a,f12.6,a,f16.5)') &
          'iter=', irest, ' tfs=', time, ' total energy=', etot
        WRITE (7, '(a,f8.4,a,f7.2,a,3f9.2)') 't=', time,&
          ' moments: monop.=', qe(1), ' dip.: ', qe(2), qe(3), qe(4)
        WRITE (6, '(a,i6,a,f12.6,a,f16.5)') &
          'iter=', irest, ' tfs=', time, ' total energy=', etot
        WRITE (6, '(a,f8.4,a,f7.2,a,3f9.2)') 't=',time,&
          ' moments: monop.=', qe(1), ' dip.: ', qe(2), qe(3), qe(4)
      END IF

      IF (myn == 0 .OR. knode == 1) CALL open_protok_el(0)

    ELSE ! case of ionic dynamics only
      ecorr = energ_ions()
    END IF 
    WRITE (*, *) '5.:cpx,y,z:', cpx(1:nion), cpy(1:nion), cpz(1:nion)
    IF (nelect > 0 .AND. itest == 0 .AND. irest == 0) &
       CALL analyze_elect(psi, rho, aloc, 0)

! intialize ionic dynamics
    IF (ionmdtyp == 1 .AND. irest == 0) THEN
! leap frog first step : propagation of momenta by half a time step
      CALL lffirststep(rho, psi)
    END IF
    ekionold = 0D0

    IF (irest == 0) THEN
      totintegprob = 0.D0
      reference_energy = etot
    END IF

  END SUBROUTINE init_dynamic

END PROGRAM tdlda_m

