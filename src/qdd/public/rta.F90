 MODULE rta_module

  USE params, ONLY: DP

! parameters for RTA calculation

  INTEGER :: jrtaint = 0 ! frequency for evaluating RTA step
  INTEGER :: itmaxDCMF = 2000 !maximum number DCMF iterations
!  INTEGER :: rtaforcetemperature = 0 ! force same temperature for both spins
  REAL(DP) :: rtamu = -1D0 ! mu from augmented Lagrangian for rho
  REAL(DP) :: rtamuj = -1D0 ! mu from augmented Lagrangian for j
  REAL(DP) :: rtaeps = 0D0 ! step size for damped gradient algorithm
  REAL(DP) :: rtae0dmp = 0D0 ! E0dmp for damped gradient algorithm
  REAL(DP) :: rtatempinit = 0D0 ! initial temperature for RTA step
  REAL(DP) :: rtasigee = -1D0 !cross section, default for Na is 6.5
  REAL(DP) :: rtars = -1D0 !Wigner-Seitz radius, default for Na 3.7
  REAL(DP) :: rtaferminit = 0.5D0 ! maximum value of T for first fermi call
  REAL(DP) :: rtadT1 = 0.001D0 ! maximum temperature step in loop fermi calls
  REAL(DP) :: rtasumvar2max = 1D-4 ! variance of s.p. energies, unit [Ry]
  REAL(DP) :: rtadiffenergmax = 0D0 ! maximum allowed energy mismatch [Ry]
  REAL(DP) :: rtaExitErr = 1D-6 ! upper limit on density deviation
  REAL(DP) :: rtaExitErrj = 0D0 ! upper limit on current deviation
                                ! if 0D0 --> set to rtaExitErr/10
  REAL(DP) :: rtaerr1 = 0D0 ! density dev. below which fermi routine called
                            ! if 0D0 --> set to  rtaExitErr*20
  REAL(DP) :: rtatempact ! minimal temperature limit, fraction of 'rtatempinit'
  REAL(DP), PARAMETER :: Tfrac=1D1 ! fraction of 'rtatempinit' for 'rtatempact'

! variables for RTA
  INTEGER :: eqstnspup, eqstnspdw

! variables and parameter for timer
  INTEGER :: itimerta(11),iraterta
  LOGICAL,PARAMETER :: ttimerta=.TRUE.     ! detailed timing information

CONTAINS
!
!_____________________________rta_________________________________________
  SUBROUTINE rta(psi, aloc, rho, iterat)
!
! Driver routine for the RTA procedure. Takes given state and applies
! one RTA step to it.
! The state is described by the set of s.p. wavefunctions 'psi'
! communicated as list parameter together with its occupation 'occup'
! communicated via module 'params'.
! The change of state is accompanied by a change of density 'rho' and
! local Kohn-Sham potential 'aloc' which both are communicated as
! list parameters.

    USE params, ONLY: DP, kdfull2, kstate, nstate, dvol, ispin, ifsicp, outnam, &
                      nyf, nxyf, centfx, nx2, ny, nz, enonlo, amoy, energy, occup, &
                      hbar, dt1, ekinsp, tfs
    USE kinetic
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate)
    REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2), rho(2*kdfull2)
    INTEGER, INTENT(IN) :: iterat

    COMPLEX(DP),ALLOCATABLE :: psiorth(:,:), psieq(:,:)
    COMPLEX(DP),ALLOCATABLE :: rhopsieq(:,:)
!    COMPLEX(DP),ALLOCATABLE :: scal(:,:)
    REAL(DP),ALLOCATABLE :: occuppsi(:), occuporth(:)
    REAL(DP) :: Estar, trace, sigee, rs, oneontaurelax, Nel0, Nel
    REAL(DP) :: Occspin(2), Occref(2)
    REAL(DP) :: Occloss(2)=(/0D0,0D0/), Occsav
    REAL(DP) :: Eperspinref(2), Eperspintarget(2)
    REAL(DP) :: mu, T, Eref, eta, EspT0
!    REAL(DP) :: entropy             ! external function
    INTEGER :: i, j, ii, nspup, nspdw
    LOGICAL, SAVE :: tinitoutput = .TRUE.

    ALLOCATE(psiorth(kdfull2, kstate))
    ALLOCATE(rhopsieq(nstate, nstate))
!    ALLOCATE(scal(nstate, nstate))
    ALLOCATE(occuppsi(kstate), occuporth(kstate))

!   set default parameters if not entered explicitely
    IF(tinitoutput) THEN
      IF(rtaerr1==0D0) rtaerr1 = rtaExitErr*100D0
      IF(rtaExitErrj==0D0) rtaExitErrj = rtaExitErr*200D0
      IF(rtadiffenergmax == 0D0) rtadiffenergmax = rtasumvar2max
      rtasumvar2max = rtasumvar2max**2     ! checking variance**2
      rtatempact = rtatempinit/Tfrac
      WRITE(*,'(a,5(/a,1pg13.5)/a,i6)') &
        'Running RTA with termination parameters:', &
        ' rtadiffenergmax=',rtadiffenergmax, &
        ' rtasumvar2max=  ',SQRT(rtasumvar2max), &
        ' rtaExitErr=     ',rtaExitErr, &
        ' rtaExitErrj=    ',rtaExitErrj, &
        ' rtaerr1=        ',rtaerr1, &
        ' itmaxDCMF=      ',itmaxDCMF
    END IF

!   copy given parameters to internal parameters
    sigee = rtasigee    !for Na 6.5D0 (rta2)
    rs = rtars          !for Na 3.7D0 (rta2)
    nspup = eqstnspup
    nspdw = eqstnspdw
    WRITE (*, *) 'RTA entered: time step nr.',iterat

!   initialize time
    IF(ttimerta) THEN
      CALL system_clock(itimerta(1), iraterta)
      CALL system_clock(itimerta(1))
    END IF

!
!   compute mean field and basic observables from scratch
!
    CALL dyn_mfield(rho, aloc, psi, 0D0, iterat)
    CALL info(psi, rho, aloc, iterat, .false.)
    occuppsi = occup
    CALL PartNumPerSpin('start rta', Occref)
    CALL calc_Eref(occup,ispin,amoy,Eperspinref)

    IF(ttimerta) CALL system_clock(itimerta(11))

!
!   transformation to natural orbitals
!
    CALL srhomat(psi, psiorth, occuporth)
    occup = occuporth       ! because dyn_mfield communicates via 'occup'

    IF(ttimerta) CALL system_clock(itimerta(2))

!
!   re-compute mean field for new basis
!
    Nel0 = sum(occup(1:nstate))         ! total number of electrons
    CALL dyn_mfield(rho, aloc, psiorth, 1D-10, iterat) ! to check after transfo

    IF(ttimerta) CALL system_clock(itimerta(10))

    WRITE(*,*) 'RTA: properties after transformation to natural orbitals'
    CALL info(psiorth, rho, aloc, iterat, .false.)
    Eref = energy
    WRITE (*, '(a,200f10.6)') 'occup ', (occup(i), i=1, nstate)
    WRITE (*, '(a,200i5)') 'ispin', (ispin(i), i=1, nstate)
    CALL PartNumPerSpin('start rta', Occref)    !memorize for final correction

    IF(ttimerta) CALL system_clock(itimerta(3))

    psi(:,1:nstate) = psiorth(:,1:nstate)
    DEALLOCATE(psiorth)

    WRITE(*,*) 'RTAinit: mfield,srhomat,mfield,info:', &
        DBLE((itimerta(11)-itimerta(1))/iraterta),&
        DBLE((itimerta(2)-itimerta(11))/iraterta),&
        DBLE((itimerta(10)-itimerta(2))/iraterta),&
        DBLE((itimerta(3)-itimerta(10))/iraterta)

!
!   compute the local-instantenous equivalent thermal state by DCMF
!
    ALLOCATE(psieq(kdfull2, kstate))
    CALL eqstate(psi, aloc, rho, psieq, occuporth, iterat)

    IF(ttimerta) CALL system_clock(itimerta(4))


!
!   several observables for the thermal equilibrium state
!
    WRITE (*, *) '_____________________ after eqstate_________________________'
    CALL PartNumPerSpin('END eqstate', Occspin)
    CALL occupT0(occup(1:nstate), amoy(1:nstate), Estar, EspT0)  ! exc. energy 'Estar'
    CALL dyn_mfield(rho, aloc, psieq, 0D0, iterat)
    CALL info(psieq, rho, aloc, iterat, .false.)
    Nel = SUM(occup(1:nstate))       ! new number of electrons

!
!   mix the required part of the thermal equilibirum state into the total state
!
    oneontaurelax = 0.4d0*sigee/rs**2*Estar/Nel     ! relaxation rate
    eta = dt1*oneontaurelax*jrtaint                 ! mixing coefficient
    Occsav = SUM(occup(1:nspup))
    CALL calcrhoeq(psi(:, 1:nspup), psieq(:, 1:nspup), &
                   psi(:, 1:nspup), occuporth(1:nspup), &
                   occup(1:nspup), nspup)           ! mixing spin-up
    Occloss(1) = Occloss(1) + Occsav-SUM(occup(1:nspup))
    Occsav = SUM(occup(nspup + 1:nstate))
    CALL calcrhoeq(psi(:, nspup + 1:nstate), psieq(:, nspup + 1:nstate), &
                   psi(:, nspup + 1:nstate), occuporth(nspup + 1:nstate), &
                   occup(nspup + 1:nstate), nspdw)  ! mixing spin-down
    Occloss(2) = Occloss(2) + Occsav-SUM(occup(nspup + 1:nstate))
    CALL PartNumPerSpin('END calcrhoeq', Occspin)
!    trace = 0D0
!    DO i = 1, nstate
!      trace = trace + REAL(scal(i, i))
!    END DO
    WRITE (*,'(a,5(1pg13.5))') &
       'RTA after mix: Nel,oneontaurelax,eta,Estar=',&
          Nel, oneontaurelax, eta, Estar
    WRITE (*,*) &
       '               Loss of particles by mixing=',Occloss

    IF(ttimerta) CALL system_clock(itimerta(5))


!
!   compute and correct total energy in two subsequent strokes
!
    WRITE(*,*) 'Observables before correction steps',i
    CALL dyn_mfield(rho, aloc, psi, 0D0, iterat)    !update field to new density
    CALL info(psi, rho, aloc, iterat, .false.)   !to have energy
    CALL PartNumPerSpin('in corrective step', Occspin)
    CALL calc_Eref(occup, ispin, amoy, Eperspinref)
    DO i=1,2
      WRITE(*,*) 'CorrectEnergy2: Energy mismatch=',Eref-energy
      IF(ABS(Eref-energy) < rtasumvar2max) EXIT
      Eperspintarget = Eperspinref + (Eref-energy)*Eperspinref/sum(Eperspinref)
      CALL CorrectEnergy2(occref(1), Eperspintarget(1),occup(1:nspup), &
                          amoy(1:nspup), occup(1:nspup), nspup)
      CALL CorrectEnergy2(occref(2), Eperspintarget(2),occup(nspup+1:nstate),&
                          amoy(nspup+1:nstate),occup(nspup+1:nstate), nspdw)
      WRITE(*,*) 'Observables after correction step nr.',i
      CALL dyn_mfield(rho, aloc, psi, 0D0, iterat)  !update field to new density
      CALL info(psi, rho, aloc, iterat, .false.)   !to have energy
      CALL PartNumPerSpin('in corrective step', Occspin)
      CALL calc_Eref(occup, ispin, amoy, Eperspinref)
    END DO
      WRITE(*,*) 'after CorrectEnergy2: Energy mismatch=',Eref-energy

    CALL temperature(mu, T)      ! effective temperaur for 'occup'


    IF(ttimerta) THEN
      WRITE(*,*) 'RTAfin final prints'
      CALL system_clock(itimerta(6))
      WRITE(*,'(10(a,f9.4/))') &
           'RTAfin prologue 1:',DBLE(itimerta(2)-itimerta(1))/iraterta,&
           'RTAfin prologue 2:',DBLE(itimerta(3)-itimerta(2))/iraterta,&
           'RTAfin core DCMF: ',DBLE(itimerta(4)-itimerta(3))/iraterta,&
           'RTAfin mix step:  ',DBLE(itimerta(5)-itimerta(4))/iraterta,&
           'RTAfin epilogue:  ',DBLE(itimerta(6)-itimerta(5))/iraterta
    END IF

    IF (tinitoutput) THEN
      OPEN (1002, STATUS='unknown', FILE='prta.'//outnam)
      WRITE (1002, '(a)') '# tfs entropy mu T'
      CLOSE (1002)
      tinitoutput = .FALSE.
    END IF
    OPEN (1002, POSITION='append', FILE='prta.'//outnam)
    WRITE (1002, '(11f14.7)') tfs, entropy(occup), mu, T
    CLOSE (1002)

    DEALLOCATE(psieq)
    DEALLOCATE(rhopsieq)
    DEALLOCATE(occuppsi, occuporth)

    RETURN

  CONTAINS
!_________________________SUBROUTINE calcrhoeq()_____________________________
    SUBROUTINE calcrhoeq(psiorthloc, psieqloc, psiloc, occuporthloc, &
                         occuploc, nstateloc)
!
! Mix actual dynamic state with state of local equilibrium to
! new dynamic state, handling states in terms of s.p. wavefunctions
! and its occupation probabilities.
! This routine is done for one sort of spin at a time.
! I/O parameters are:
!  psiorthloc   = set of actual s.p. wavefunctions (before DCMF)
!  occuporthloc = occupations of the actual wavefunctions
!  psieqloc     = s.p. wavefunctions of local eequilibrium state
!  pisloc       = new set of s.p. wavefunctions after mixing
!  occuploc     = at entry occupations for local equilibrium state
!                 at exit occupations for the new set of wavefunctions
!  nstateloc    = nr. of s.p. states

      USE util, ONLY: wfovlp
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: nstateloc
      COMPLEX(DP),INTENT(IN) :: psiorthloc(kdfull2, nstateloc)
      COMPLEX(DP),INTENT(IN) :: psieqloc(kdfull2, nstateloc)
      COMPLEX(DP),INTENT(OUT) :: psiloc(kdfull2, nstateloc)
      REAL(DP),INTENT(IN) :: occuporthloc(nstateloc)
      REAL(DP),INTENT(IN OUT) :: occuploc(nstateloc)

      REAL(DP) :: coeffoccup
      COMPLEX(DP),ALLOCATABLE :: scal0(:,:), scal1(:,:), tab(:,:)
!      COMPLEX(DP),ALLOCATABLE :: taborth(:,:)
      COMPLEX(DP),ALLOCATABLE :: occloc(:,:)
      COMPLEX(DP),ALLOCATABLE :: uloc(:,:), vloc(:,:), vploc(:,:)
      COMPLEX(DP),ALLOCATABLE :: wloc(:,:), D1onrholoc(:,:)
      REAL(DP),ALLOCATABLE :: eigenloc(:)
      INTEGER indx(2*nstateloc), i, j, ncall, k, l

      WRITE (*, *) 'CALCRHOEQ: nstateloc=', nstateloc

      !allocate matrices
      ALLOCATE(scal0(2*nstateloc, 2*nstateloc))
      ALLOCATE(scal1(2*nstateloc, 2*nstateloc))
      ALLOCATE(tab(kdfull2, 2*nstateloc))
!      ALLOCATE(taborth(kdfull2, 2*nstateloc))
      ALLOCATE(occloc(2*nstateloc, 2*nstateloc))
      ALLOCATE(uloc(2*nstateloc, 2*nstateloc))
      ALLOCATE(vloc(2*nstateloc, 2*nstateloc))
      ALLOCATE(vploc(2*nstateloc, 2*nstateloc))
      ALLOCATE(wloc(2*nstateloc, 2*nstateloc))
      ALLOCATE(D1onrholoc(2*nstateloc, 2*nstateloc))
      ALLOCATE(eigenloc(2*nstateloc))



      ! tab is the double set of s.p. wave functions
      tab(:, 1:nstateloc) = psiorthloc(:, 1:nstateloc)
      tab(:, nstateloc + 1:2*nstateloc) = psieqloc(:, 1:nstateloc)

      ! matrix of overlaps for the double set
      CALL OverlapMatrix2(tab, tab, scal1, 'double set', 2*nstateloc)

      ! initialize occloc as diagonal matrix of occupation numbers
      occloc = 0D0
      DO i = 1, nstateloc
        occloc(i, i) = (1 - eta)*occuporthloc(i)
        occloc(nstateloc + i, nstateloc + i) = eta*occuploc(i)
      END DO
      ! multiply by oveerlap matrix from booth sides
      occloc = matmul(transpose(conjg(scal1)), occloc)
      occloc = matmul(occloc, scal1)

      ! diagonalize overlap matrix
      CALL cdiagmat(scal1, eigenloc, vloc, 2*nstateloc)
      ! compute diagonal matrices D1on_rho
      D1onrholoc = 0D0
      DO i = 1, 2*nstateloc
        D1onrholoc(i, i) = 1D0/sqrt(eigenloc(i))
      END DO
      ! scale eigenvectors 'v' into 'vp' to have a unity matrix
      vploc = matmul(vloc, D1onrholoc)
      ! use 'vp' to compute Vp^T*occloc*Vp
      scal0 = matmul(transpose(conjg(Vploc)), occloc)
      scal0 = matmul(scal0, vploc)

      ! diagonalize that
      CALL cdiagmat(scal0, eigenloc, uloc, 2*nstateloc)
      ! compose to transformation matrix for change of basis
      wloc = matmul(vploc, uloc)

      ! new occupation numbers
      coeffoccup = sum(eigenloc(1:2*nstateloc))
      ! ordering of the occupation numbers
      ncall = 2*nstateloc
      CALL indexx(ncall, eigenloc, indx)
      WRITE (*, *) 'IN calcrhoeq after indxx, new occupations sorted:'
      WRITE (*, '(10f14.10)') (eigenloc(indx(i)), i=2*nstateloc, 1, -1)

!     new occupations
      j = 0
      DO i = 2*nstateloc, nstateloc + 1, -1
        j = j + 1
        occuploc(j) = eigenloc(indx(i))
      END DO
      l = 0

      ! compute the new basis functions phi_j in 'taborth'
!      taborth = 0D0
      psiloc = 0D0
#if(omp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k,l) SCHEDULE(STATIC)
#endif
      DO i = 2*nstateloc, nstateloc + 1, -1
        j = 2*nstateloc-i+1 ! j + 1
        k = indx(i)
        DO l = 1, 2*nstateloc
          psiloc(:, j) = psiloc(:, j) + tab(:, l)*wloc(l, k)
        END DO
      END DO
#if(omp)
!$OMP END PARALLEL DO
#endif


!#if(omp)
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,i) SCHEDULE(STATIC)
!#endif
!      DO j = 1, 2*nstateloc
!        DO i = 1, 2*nstateloc
!          taborth(:, j) = taborth(:, j) + wloc(i, j)*tab(:, i)
!        END DO
!      END DO
!#if(omp)
!!$OMP END PARALLEL DO
!#endif
!      ! store the nstateloc states with highest occupation numbers
!#if(omp)
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,i) SCHEDULE(STATIC)
!#endif
!      DO i = 2*nstateloc, nstateloc + 1, -1
!        j = 2*nstateloc-i+1 ! j + 1
!        psiloc(:, j) = taborth(:, indx(i))
!      END DO
!#if(omp)
!!$OMP END PARALLEL DO
!#endif


      DEALLOCATE(scal0)
      DEALLOCATE(scal1)
      DEALLOCATE(tab)
!      DEALLOCATE(taborth)
      DEALLOCATE(occloc)
      DEALLOCATE(uloc)
      DEALLOCATE(vloc)
      DEALLOCATE(vploc)
      DEALLOCATE(wloc)
      DEALLOCATE(D1onrholoc)
      DEALLOCATE(eigenloc)


    END SUBROUTINE calcrhoeq

  END SUBROUTINE rta

!________________________________________cdiagmat_________________________
  SUBROUTINE cdiagmat(mat, eigen, vect, N)
!mat : complex matrix to be diagonalised.
!eigen :vector of eigenvalues
!vect :eigenvectors written in columns
!N :dimension, assuming that all states are occupied
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(1D0) ! PRECISION setting
    INTEGER, INTENT(IN) :: N
    COMPLEX(DP), INTENT(IN) :: mat(N, N)
    REAL(DP), INTENT(OUT) :: eigen(N)
    COMPLEX(DP), INTENT(OUT) :: Vect(N, N)

    COMPLEX(DP) :: mat0(N, N)
    INTEGER :: i, j
    vect = 0D0
    eigen = 0D0
    mat0 = mat
    CALL HEigensystem(N, mat0, N, eigen, vect, N, 0)
    vect = transpose(conjg(vect))

  END SUBROUTINE cdiagmat


  SUBROUTINE calc_Eref(occup, ispin, Ei, Eref)
!
! sum of s.p. energies accumuklated per spin

    USE params, ONLY: DP, nstate
    IMPLICIT NONE
    REAL(DP), INTENT(IN):: occup(nstate), Ei(nstate)
    REAL(DP), INTENT(OUT):: Eref(2)
    INTEGER, INTENT(IN) ::ispin(nstate)

    INTEGER:: i

    Eref = 0D0
    DO i = 1, nstate
      Eref(ispin(i)) = Eref(ispin(i)) + occup(i)*Ei(i)
    END DO

  END SUBROUTINE calc_Eref
!_________________________________eqstate_______________________________________

  SUBROUTINE eqstate(psi, aloc, rho, psi1, occuporth, iterat)
!update psi into psi1, and occuporth(entrance) into occup(IN params)
    USE params, ONLY: DP, kdfull2, kstate, nstate, dvol, ispin, ifsicp, outnam, &
                      nyf, nxyf, centfx, nx2, ny, nz, enonlo, &
                      amoy, energy, occup, nxsh, nysh, nzsh, h2m
    USE kinetic
    USE util, ONLY: prifld, wfnorm
    IMPLICIT NONE
    COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2, kstate), psi1(kdfull2, kstate)
    REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2), rho(2*kdfull2), occuporth(kstate)
    INTEGER, INTENT(IN) :: iterat

    COMPLEX(DP) :: dpsi(kdfull2)

    REAL(DP) :: rhotot0(kdfull2), rhotot1(kdfull2), ERR, ma0, ma1, err1
    REAL(DP) :: lambda(kdfull2), mu, time0, time1
    REAL(DP) :: curr0(kdfull2, 3), curr1(kdfull2, 3), lambdaj(kdfull2, 3)
    REAL(DP) :: errj, muj !parameters for current
    REAL(DP) :: sumvar2, eal, fac, EspPerMod(nstate), EspTotRef, temp(2)
    REAL(DP) :: EspPerSpinRef(2), EspPerSpinAchieved(2), EspPerSpinTarget(2)
    REAL(DP) :: mut, EspTotAchieved, dummy, ekincoll0, ekincoll1, errcoll1
    REAL(DP) :: errf1, errf2, dT1, dT2 ! local parameters, to fermi
    REAL(DP) :: rho2aver
    INTEGER :: i, j, ishift, ii, nspup, nspdw, ind
    LOGICAL :: topen
    INTEGER, EXTERNAL :: conv3to1
    LOGICAL, PARAMETER :: tdetails = .TRUE. ! switch to more control output
    LOGICAL, PARAMETER :: toccupT0 = .FALSE. ! switch to T=0D0 initial occup


    WRITE (*, *) "EQSTATE: iterat=", iterat
    IF(ttimerta) CALL system_clock(itimerta(4))
    nspup = eqstnspup
    nspdw = eqstnspdw
    EspPerMod(1:nstate) = amoy(1:nstate)

    ! initialize psi, rhotot,lambda,j lamdaj
    sumvar2 = 1D3          ! a big value to start with
    ERR = 1D3
    ERRj = 1D3
    psi1 = psi             ! copy set of wf's to work space
    lambda = 0D0
    lambdaj = 0D0
    mu = rtamu
    muj = rtamuj

    ! internal convergence parameters.
    errf1 = rtaerr1
    errf2 = errf1/2D0
    dT1 = rtadT1
    dT2 = dT1*2.5D0

!    CALL calcrhotot(rhotot0, psi)
    rhotot0 = rho(1:kdfull2)
    rhotot1 = rhotot0
    CALL calc_current(curr0, psi)
    ekincoll0 = 0D0
    DO ind=1,kdfull2
      ekincoll0 = ekincoll0 + &
                  (curr0(ind,1)**2+curr0(ind,2)**2+curr0(ind,3)**2)/&
                  MAX(rho(ind),1D-20)
    END DO
    ekincoll0 = ekincoll0*dvol*h2m
    curr1 = curr0
    ma0 = dvol*SUM(rhotot0)
    rho2aver = SUM(rhotot0**2)*dvol    ! reference to regularize square current
    WRITE(*,*) 'j ma0:',j,ma0
    WRITE (*, *) ' j centfx ERR errj mu muj  ma0 ma1'
    j = 0
! ALOC remains unchanged throughout the DCMF iterations
!    CALL coul_mfield(rho)

    CALL calc_Eref(occuporth,ispin,amoy,EspPerSpinRef) ! reference total sp en.
    EspTotRef = EspPerSpinRef(1) + EspPerSpinRef(2)
    WRITE (*, *) 'before Fermi, j,EspTotRef=',j,EspTotRef
     IF(toccupT0) THEN
       CALL OccT1(SUM(occup(1:nspup)), amoy(1:nspup), &
                  EspPerSpinAchieved(1), mut, dummy, nspup, &
                  rtatempact, occup(1:nspup))
       CALL OccT1(SUM(occup(nspup+1:nstate)), amoy(nspup+1:nstate), &
                  EspPerSpinAchieved(2), mut, dummy, nspdw, &
                  rtatempact, occup(nspup+1:nstate))
     ELSE
       CALL fermi_total(amoy, SUM(EspPerSpinRef), occup, &
                        rtatempact, rtaferminit, temp(1), mut)
     END IF
    !order against sp Energy
    CALL OrderPsiOccEsp(nspup, psi1(:, 1:nspup), occup(1:nspup), amoy(1:nspup))
    CALL OrderPsiOccEsp(nspdw, psi1(:, nspup + 1:nstate), &
                        occup(nspup + 1:nstate), amoy(nspup + 1:nstate))

    WRITE (*,'(a,i5,4(1pg13.5))') 'Fermi finished: j,entropy,sum s.p.en.=',&
        j,entropy(occup),EspPerSpinRef
    EspPerSpinTarget = EspPerSpinRef  !initial value of EspPerSpinTarget
    EspTotAchieved = 0D0              !to avoid EXIT at first pass

    INQUIRE (UNIT=177, OPENED=topen)
    IF (.NOT. topen) THEN
      IF (tdetails) THEN
        OPEN (177, FILE='convergence_RTA.'//outnam)
        OPEN (1000, STATUS='unknown', FORM='FORMATTED', FILE='pspeed.'//outnam)
        WRITE (1000, '(a)') &
          '# final densities and currents after DCMF ot one reference point'
        WRITE (1000, '(a)') 'nx-grid rho-up_goal rho-down_goal rho-up_obt. '// &
          ' rho-down_obt. curr-up_goal curr-up_obt.'
        CLOSE (1000)
      END IF
      OPEN (1001, STATUS='unknown', FILE='peqstate.'//outnam)
      WRITE (1001, '(a)') '# final precisions after DCMF step', &
        '# nr.step nr.iterat s.p.variance rho_error j_error mue muej energy'
      CLOSE (1001)
    END IF
    IF(ttimerta) THEN
      CALL system_clock(itimerta(5))
      WRITE(*,'(a,f9.4)') 'EQSTATE: prologue:', &
                  DBLE(itimerta(5)-itimerta(4))/iraterta
    END IF
    WRITE (177, '(/a,i5/a)') '#RTA protocol at iterat=', iterat, &
      ' j ERR errj diffrho sumj0 sumj1 diffE sumvar2 entropy'//&
      ' j-error**2  error-ekin   ekin'

    ! now the main loop for the RTA iteration
    DO WHILE ( &
              ( sumvar2 > rtasumvar2max &
                .OR. ABS(EspTotRef - EspTotAchieved) >  rtadiffenergmax &
                .OR. ERR > rtaExitErr  &
                .OR. errj > rtaExitErrj &
                .OR. j<6                    ) &
                .AND. &
              ( j .LE. itmaxDCMF  &
                .OR. &
                (ABS(EspTotRef - EspTotAchieved) >  rtadiffenergmax &
                 .AND. j .LE. 2*itmaxDCMF )  ) &
             )
      IF(ttimerta) CALL system_clock(itimerta(4))
      mu = rtamu

!
! at high level of convergence fully update occupation numbers at every step
!
      IF (j>10 .AND. &
          ((sumvar2 < rtasumvar2max*1D3 .AND. ERR < errf2) &
           .OR. ERR < errf2/4D0)   ) THEN
        mu = rtamu*2
        EspPerSpinAchieved = 0D0
        DO ii = 1, nstate
          EspPerSpinAchieved(ispin(ii)) = &
               EspPerSpinAchieved(ispin(ii)) + occup(ii)*EspPerMod(ii)
        END DO
        WRITE (*,'(a,i4,20f14.7)') 'iter,EspPerSpinAchieved-EspPerSpinRef', &
            j,EspPerSpinAchieved(1) - EspPerSpinRef(1), &
            EspPerSpinAchieved(2) - EspPerSpinRef(2)
!          CALL fermi_total(EspPerMod, SUM(EspPerSpinTarget), occup, &
          CALL fermi_total(EspPerMod,  EspTotRef, occup, &
                           max(temp(1) - dT2, rtatempact), &
                           temp(1) + dT2, temp(1), mut)
        EspTotAchieved = SUM(EspPerMod(1:nstate)*occup(1:nstate))
        IF(tdetails) THEN
          WRITE(*,'(a,i5,4(1pg13.5))') 'Fermi 1: j,entropy,T,E-diffs=',&
              j,entropy(occup),temp(1),&
              EspTotAchieved- EspTotRef,&
              eal- EspTotRef
        END IF

!
! at preliminary level of convergence evaluate new occupations
! only in small window of temperatures and only every 10. step
!
      ELSE IF((sumvar2 < rtasumvar2max*1D5 .AND. ERR < errf1) &
              .AND. MOD(j, 10) == 0) THEN
          CALL fermi_total(EspPerMod, EspTotRef, occup, &
                             max(temp(1) - dT1, rtatempact), &
                             temp(1) + dT1, temp(1), mut)
        EspTotAchieved = SUM(EspPerMod(1:nstate)*occup(1:nstate))
        WRITE(*,*) 'Fermi 2: j,entropy,T,E-diff=',j,entropy(occup),&
            temp(1),EspTotAchieved-EspTotRef

      END IF

      j = j + 1

      ! iterate set of s.p. wavefunctions
      IF(ttimerta) CALL system_clock(itimerta(5))
      CALL calc_psi1(psi1, aloc, rhotot0, rhotot1, curr0, curr1, j, &
                     lambda, mu, lambdaj, muj, sumvar2, eal, EspPerMod)
      IF(ttimerta) CALL system_clock(itimerta(6))

      ! update Lagrangian fields according to augmented Lagrangian method
      lambda = lambda + (rhotot1 - rhotot0)*2D0*mu
      lambdaj = lambdaj + (curr1 - curr0)*2D0*muj

      ii = conv3to1(nxsh, nysh, nzsh)
      WRITE(*,'(a,i4,5(1pg13.5))') 'lambda observables:',j,&
         SUM(lambda)*dvol,SUM(ABS(lambda))*dvol,MAXVAL(lambda),&
         MINVAL(lambda),rhotot1(ii)-rhotot0(ii)

!     scale to relative errors
      ERR = SUM(ABS(rhotot1 - rhotot0))*dvol/ma0
!old      errj = SUM(ABS(curr1-curr0)**2)/MAX(SUM(ABS(curr0)**2),ma0*1D-4/dvol)
      errj = SUM(ABS(curr1 - curr0)**2)/MAX(SUM(ABS(curr0)**2),1D-10)

      ma1 = dvol*SUM(rhotot1)
      IF (mod(j, 1) .eq. 0) THEN
        IF (tdetails) WRITE (*, '(a,i6,20f14.7)') ' DCMF-iteration:', &
          j, centfx, ERR, errj, mu, muj, ma0, ma1, EspTotRef, EspTotAchieved &
          , EspTotRef - EspTotAchieved, SQRT(sumvar2), &
          EspPerSpinAchieved(1) - EspPerSpinAchieved(2)
        ekincoll1 = 0D0
        DO ind=1,kdfull2
          ekincoll1 = ekincoll1 + &
                     (curr1(ind,1)**2+curr1(ind,2)**2+curr1(ind,3)**2)/&
                     MAX(rhotot1(ind),1D-20)
        END DO
        ekincoll1 = ekincoll1*dvol*h2m
        errcoll1 = 0D0
        DO ind=1,kdfull2
          errcoll1 = errcoll1 + &
                     ((curr1(ind,1)-curr0(ind,1))**2+&
                      (curr1(ind,2)-curr0(ind,2))**2+&
                      (curr1(ind,3)-curr0(ind,3))**2)/&
                     MAX(rhotot1(ind),1D-20)
        END DO
        errcoll1 = errcoll1*dvol*h2m
        WRITE(*,'(a,i5,6(1pg13.5))') '                j,ekins:',&
          j,errcoll1,ekincoll0,ekincoll1,&
          ABS(ekincoll0-ekincoll1)/ekincoll0,errcoll1/ekincoll0,errj
        WRITE(*,'(a,3(1pg13.5))') 'norms,occup:',SUM(occup(1:nstate)) &
          ,SUM(occup(1:nspup)),SUM(occup(nspup+1:nstate))
        WRITE(*,'(400(1pg13.5))') occup(1:nstate)
        FLUSH (6)
        WRITE (177, '(a,i6,20(1pg13.5))') ' DCMF-iteration:', &
          j, ERR, errj, ma0 - ma1, SUM(abs(curr0))*dvol/rho2aver, &
          SUM(abs(curr1))*dvol/rho2aver,&
          EspTotRef - EspTotAchieved, SQRT(sumvar2),&
          entropy(occup),errcoll1/ma0,ABS(ekincoll0-ekincoll1)/ma0,ekincoll0
        FLUSH (177)
      END IF

      IF(ttimerta) THEN
        CALL system_clock(itimerta(7))
        WRITE(*,'(9(a,f9.4/))') &
                'LOOP prologue: ',DBLE(itimerta(5)-itimerta(4))/iraterta,&
                'LOOP calc_psi1:',  DBLE(itimerta(6)-itimerta(5))/iraterta,&
                'LOOP epilogue: ',  DBLE(itimerta(7)-itimerta(6))/iraterta
      END IF

    END DO   ! end of main DCMF loop

     amoy(1:nstate) = EspPerMod(1:nstate)
     dummy = ABS(EspTotRef - EspTotAchieved)
     IF(dummy >  rtadiffenergmax) &
       WRITE(*,'(a,1pg13.5))') 'DCMF state: unresolved energy mismatch=',dummy

    IF (tdetails) THEN
      OPEN (1000, STATUS='unknown', POSITION='append', FORM='FORMATTED', FILE='pspeed.'//outnam)
      WRITE (1000, '(a,i8)') '# at dynamic step nr.:', iterat
      IF(dummy >  rtadiffenergmax) &
       WRITE(1000,'(a,1pg13.5))') 'DCMF state: unresolved energy mismatch=',dummy
      WRITE (1000, '(i6,6(1pg13.5))') (ii, rhotot0(ny*(nyf + nxyf) + ii), &
                                       rhotot0(ny*(nyf + nxyf) + ii), &
                                       rhotot1(ny*(nyf + nxyf) + ii), &
                                       rhotot1(ny*(nyf + nxyf) + ii), &
                                       curr0(ny*(nyf + nxyf) + ii,1), &
                                       curr1(ny*(nyf + nxyf) + ii,1), &
                                       ii=1, nx2)
      WRITE (1000, *)
      CLOSE (1000)
    END IF
    CALL cpu_time(time1)
    WRITE (*, *) 'j,time', j, (time1 - time0)
    CALL dyn_mfield(rho, aloc, psi1, 0D0, iterat)
    WRITE (*, *) '_____________in eqstate 2: iterat=', iterat
    CALL info(psi1, rho, aloc, iterat, .false.)
    WRITE (*, '(2i4,11(1pg13.5))') iterat, j, sumvar2, ERR, &
        errj, mu, muj, energy
    OPEN (1001, POSITION='append', FILE='peqstate.'//outnam)
    WRITE (1001, '(2i6,11(1pg13.5))') iterat, j, SQRT(sumvar2), ERR, &
        errj, mu, muj, energy
      IF(dummy >  rtadiffenergmax) &
       WRITE(1001,'(a,1pg13.5))') 'unresolved energy mismatch=',dummy
    CLOSE (1001)

    RETURN

  END SUBROUTINE eqstate


  SUBROUTINE OccupT0(occloc, esploc, Estar, EspT0)
!
! This routine computes occupation numbers at temperature=0 for each spin
! in order of s.p. energies, putting occupation numbers to 1, until
! particle number is reached.
! Then computes excitation energy Estar as output of the routine.
! The calling occupation numbers 'occloc' are not changed.
!
    USE params, ONLY: DP, ispin, nstate
    IMPLICIT NONE

    REAL(DP),INTENT(IN) :: occloc(nstate), esploc(nstate)
    REAL(DP),INTENT(OUT) :: Estar, EspT0

    INTEGER :: isp, i, n, indx(nstate), order(nstate)
    REAL(DP) :: Etot, OccTot, EspSpin(nstate), OccSpin(nstate), occup
    REAL(DP) :: EspSpinT0

    LOGICAL, PARAMETER :: ttestprint=.FALSE.
! first compute tables of energies per spin
! here without assumptionn that spins are in blocks
    EspSpinT0 = 0D0
    Etot = dot_product(Occloc(1:nstate), Esploc(1:nstate))
    DO isp = 1, 2
      n = 0
      DO i = 1, nstate
        IF (ispin(i) .eq. isp) THEN
          n = n + 1
          OccSpin(n) = occloc(i)
          EspSpin(n) = EspLoc(i)
          indx(N) = i
        END IF
      END DO
      OccTot = SUM(Occspin(1:N))
      CALL indexx(n, EspSpin, order)

      ! now we have energies per spin in energy order
      ! do a second loop to populate the occupation numbers
      DO i = 1, n
        IF (OccTot >= 1D0) THEN
          occup = 1D0
          OccTot = OccTot - 1D0
        ELSEIF (Occtot < 1D0 .and. OccTot >= 0) THEN
          occup = OccTot
          OccTot = OccTot - 1D0
        ELSE
          occup = 0D0
          OccTot = OccTot - 1D0
        END IF
        EspSpinT0 = EspSpinT0 + occup*EspLoc(indx(order(i)))

        IF(ttestprint) WRITE (*, '(a,4f14.8)') 'occup, Esp', &
           occup, occloc(indx(order(i))), EspLoc(indx(order(i)))
      END DO
    END DO
    Estar = Etot - EspSpinT0
    EspT0 =  EspSpinT0
    WRITE (*, '(a,3f14.7)') 'OCCUPT0: estar,Etot,EspSpinT0=', &
          estar, Etot, EspSpinT0

  END SUBROUTINE OccupT0

!_______________________________Calc_psi1____________________________________

! one DCMF step for the s.p. wavefunctions, including update of mean field
!
  SUBROUTINE Calc_psi1(psi1, aloc, rhotot0, rhototloc, curr0, curr1, j, &
                       lambda, mu, lambdaj, muj, sumvar2, eal, EspMod)
#if(dynomp)
    USE params, ONLY: DP, kdfull2, kstate, nstate, dvol, ispin, eye, ipsptyp, &
                      h2m, hbar, ame, pi, nx2, ny2, nz2, occup, enonlo, nxyz, &
                      omp_get_thread_num
#else
    USE params, ONLY: DP, kdfull2, kstate, nstate, dvol, ispin, eye, ipsptyp, &
                      h2m, hbar, ame, pi, nx2, ny2, nz2, occup, enonlo, nxyz
#endif
    USE kinetic
    USE util, ONLY: wfovlp,c_wfovlp
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: rhotot0(kdfull2)
    REAL(DP), INTENT(IN) :: mu, muj
    INTEGER, INTENT(IN):: j
    COMPLEX(DP), INTENT(IN OUT) :: psi1(kdfull2, kstate)
    REAL(DP), INTENT(IN OUT) :: EspMod(nstate)      ! s.p. energies
    REAL(DP), INTENT(IN OUT) :: aloc(2*kdfull2)
    REAL(DP), INTENT(IN OUT) :: rhototloc(kdfull2)
    REAL(DP), INTENT(IN OUT) :: curr0(kdfull2, 3), curr1(kdfull2, 3)
    REAL(DP), INTENT(IN OUT) :: lambda(kdfull2)
    REAL(DP), INTENT(IN OUT) :: lambdaj(kdfull2, 3)
    REAL(DP), INTENT(OUT) :: eal      ! sum of Lagr. s.p. energies
    REAL(DP), INTENT(OUT) :: sumvar2

    COMPLEX(DP),ALLOCATABLE :: psik(:,:)     ! wf in k space
    COMPLEX(DP),ALLOCATABLE :: ALpsi(:,:)    ! DCMF potential acting on wf
    COMPLEX(DP),ALLOCATABLE :: vpsi(:,:)     ! KS potential acting on wf
    COMPLEX(DP),ALLOCATABLE :: Apsi(:,:)     ! workspace
!    COMPLEX(DP),ALLOCATABLE :: hpsip(:,:)    ! save H*psi in k space
    COMPLEX(DP),ALLOCATABLE :: Apot(:,:)   ! Lagrangrian for density and current

    REAL(DP) :: epswf, e0dmp, deltalambda, EspLag(nstate), epot(nstate)
    REAL(DP) :: sumvar2off, dkvol, eee, ekin
    REAL(DP):: hbarm
    INTEGER :: N, ii, ishift, ii1, ithr, nacthr, ntm, nta
    INTEGER,SAVE :: iprintcount=0    ! counter for printing HPSI information

    INTEGER,PARAMETER :: mfreq_var=1      ! frequency of variance evaluation
    LOGICAL,PARAMETER :: ttestprint=.FALSE.

    ! to compute new s.p. energies after DCMF step
    LOGICAL,PARAMETER :: tactual_energy=.FALSE.

#if(omp)
    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS,OMP_GET_NUM_THREADS
#endif

    ithr = 0
    hbarm = hbar/ame
    dkvol = ((2*pi)**3)*dvol
    epswf = rtaeps
    e0dmp = rtae0dmp
    eal = 0D0
    WRITE(*,*) 'CALC_PSI1:',SUM(rhotot0)*dvol,SUM(rhototloc)*dvol,&
       SUM(aloc)*dvol

!    itimerta(7) = 0
!    itimerta(8) = 0
!    itimerta(10) = itimerta(6)
    nacthr = numthr-1

    ALLOCATE(Apot(kdfull2,0:3))
    ALLOCATE(Apsi(kdfull2,0:nacthr))
    ALLOCATE(ALpsi(kdfull2,0:nacthr))
!    ALLOCATE(vpsi(kdfull2,0:nacthr))
!    ALLOCATE(hpsip(kdfull2,0:nacthr))
!    ALLOCATE(psik(kdfull2,0:nacthr))

!    ALLOCATE(hpsi(kdfull2,kstate))
!    ALLOCATE(scal(nstate,nstate), psi1new(kdfull2,kstate), HpsiN(kdfull2))

    IF(ttimerta) CALL system_clock(itimerta(6))

!   prepare effective constraint fields
    Apot(:,0) = (lambda(:)+2D0*mu*(rhototloc(:)-rhotot0(:)))   ! density constr.
    Apot(:,1) = eye*(2*muj*(curr1(:, 1) - curr0(:, 1)) + lambdaj(:, 1))
    Apot(:,2) = eye*(2*muj*(curr1(:, 2) - curr0(:, 2)) + lambdaj(:, 2))
    Apot(:,3) = eye*(2*muj*(curr1(:, 3) - curr0(:, 3)) + lambdaj(:, 3))
!   divergence of j-constraints, use 'vpsi' and 'Apsi' as workspaces
!    WRITE(*,'(a,i3,4(1pg13.5))') 'Apot etc:',j, &
!      (REAL(wfovlp(CONJG(Apot(:,ii)),Apot(:,ii))),ii=0,3)

     IF(ttimerta) CALL system_clock(itimerta(11))

!    IF(MOD(j-1,mfreq_var)==0) THEN
      sumvar2 = 0D0
      sumvar2off = 0D0
!    END IF
!    WRITE(*,*) 'XGRADIENT start:'
    FLUSH(6)
#if(dynomp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ithr) SCHEDULE(STATIC) REDUCTION(+: eal,sumvar2,sumvar2off)
#endif
    DO N = 1, nstate
#if(dynomp)
      ithr = OMP_GET_THREAD_NUM()
      IF(N==1) THEN
        ntm = OMP_GET_MAX_THREADS()
        nta = OMP_GET_NUM_THREADS()
        IF(ntm.NE.nta) WRITE(*,*) &
        'RTA PSI1: threads not exhausted: N,ithr,maxthreads,actthreads=',&
         N,ithr,ntm,nta
      END IF
      IF(ttestprint) THEN
        WRITE(*,*) 'PSI1: N,ithr,norm,numthr=',&
            N,ithr,wfovlp(psi1(:,N),psi1(:,N)),&
            OMP_GET_MAX_THREADS(),OMP_GET_NUM_THREADS()
      END IF
#endif

      IF(N==1 .AND. ttimerta) CALL system_clock(itimerta(10))

!     action of the constrained mean-field Hamiltonian on the s.p. wf
!     returns action solf Kohn-Sham Hamiltonian on woorkspace 'Apsi'
!     returns action of constraint on workspace 'ALpsi'
      CALL calc_hamiltonien(N)

      IF(N==1 .AND. ttimerta) CALL system_clock(itimerta(7))

      ! transformations into k space
!      CALL fftf(psi1(:, N), psik(:,ithr))         ! s.p. wavefunction
      CALL fftf(psi1(:, N), psi1(:, N))         ! s.p. wavefunction
      CALL fftf(ALpsi(:,ithr), ALpsi(:,ithr))     ! full DCMF potential on wf
      CALL fftf(Apsi(:,ithr), Apsi(:,ithr))       ! potential acting on wf

      ! m.f. Ham. in k-space
!      Apsi(:,ithr) = Apsi(:,ithr)+h2m*akv*psik(:,ithr)
      Apsi(:,ithr) = Apsi(:,ithr)+h2m*akv*psi1(:, N)
      ! new s.p. energy
      EspMod(N) = REAL(dot_product(Apsi(:,ithr), psi1(:, N)))*dkvol

      ! action of full DCMF Hamilt.
      Apsi(:,ithr) = Apsi(:,ithr) + ALpsi(:,ithr)
      ! s.p. energy with constraint
      EspLag(N) = REAL(dot_product(Apsi(:,ithr), psi1(:, N)))*dkvol
      ! subtract expectation value
      Apsi(:,ithr) = Apsi(:,ithr) - EspLag(N)*psi1(:, N)
      ! add up to sum of s.p. energies
      eal = eal + EspLag(N)*occup(N)
      ! add up to variance
      sumvar2 = sumvar2 + occup(N)* &
             REAL(dot_product(Apsi(:,ithr),Apsi(:,ithr)))*dkvol

!      WRITE(*,*) 'EVARSP-P:',N,&
!           REAL(dot_product(Apsi(:,ithr),Apsi(:,ithr)))*dkvol

      IF(N==1 .AND. ttimerta) CALL system_clock(itimerta(8))

!      IF(MOD(j-1,mfreq_var)==0) THEN
!        hpsip(:,ithr) = Apsi(:,ithr)
!        CALL calc_varsp(hpsip(:,ithr), psi1, sumvar2, sumvar2off, N) ! variance
!      END IF

      !gradient step in k space and backtransformation to real space
      psi1(:, N) = psi1(:, N) - epswf/(e0dmp + akv)*(Apsi(:,ithr))
      CALL fftback(psi1(:, N), psi1(1, N))

      IF(N==1 .AND. ttimerta .AND. j<3) THEN
        CALL system_clock(itimerta(9))
        WRITE(*,'(9(a,f9.4/))') &
                 'HPSI 1:    ',DBLE(itimerta(7)-itimerta(10))/iraterta,&
                 'HPSI 2:',DBLE(itimerta(8)-itimerta(7))/iraterta,&
                 'HPSI 3: ',DBLE(itimerta(9)-itimerta(8))/iraterta
!        itimerta(10) = itimerta(9)
      END IF

    END DO    ! end loop over states N
#if(dynomp)
!$OMP END PARALLEL DO
#endif
    DEALLOCATE(ALpsi)
    DEALLOCATE(Apot)

!    WRITE(*,*) 'XGRADIENT end:'
    FLUSH(6)

    IF(MOD(j-1,mfreq_var)==0) THEN
      sumvar2 = sumvar2/SUM(occup(1:nstate)) ! normalize per particle
      sumvar2off = sumvar2off/SUM(occup(1:nstate)) ! normalize per particle
    END IF

    IF(ttimerta) CALL system_clock(itimerta(7))

    CALL cschmidt(psi1)   ! ortho-normalization

    IF(ttimerta) CALL system_clock(itimerta(8))

    IF(MOD(j-1,mfreq_var)==0) THEN
!      CALL calc_var(hpsip, psi1, sumvar2, sumvar2off)  ! variance
      WRITE (*, '(a,2(1pg13.5))') 's.p. variances:',SQRT(sumvar2),SQRT(sumvar2off)
      WRITE (*, '(a,(21f10.6))') 'occup', (occup(N), N=1, nstate)
      WRITE (*, '(a,i5,(200f10.6))') 'j,EspLag (h+Lagr.)',j,(EspLag(N), N=1, nstate)
      WRITE (*, '(a,i5,(200f10.6))') 'j,EspMod (h only) ',j,(EspMod(N), N=1, nstate)
    END IF


    ! new mean field and expectation values
    CALL calcrhotot(rhototloc, psi1)
    CALL calc_current(curr1, psi1)

  IF(tactual_energy) THEN
!
! new s.p. energy, use 'Apsi' as work space for action of m.f. potential
!
#if(dynomp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ishift,ii,ithr) SCHEDULE(STATIC)
#endif
    DO N = 1, nstate
#if(dynomp)
      ithr = OMP_GET_THREAD_NUM()
      IF(N==1) THEN
        ntm = OMP_GET_MAX_THREADS()
        nta = OMP_GET_NUM_THREADS()
        IF(ntm.NE.nta) WRITE(*,*) &
        'RTA SPE: threads not exhausted: N,ithr,maxthreads,actthreads=',&
         N,ithr,ntm,nta
      END IF
      IF(ttestprint) WRITE(*,*) 'new s.p. energies: ithr=',ithr
#else
      ithr = 0
#endif
      CALL calc_ekin(psi1(:, N), EspMod(N))
      ishift = (ispin(N) - 1)*nxyz
      Apsi(:,ithr) = 0D0
      IF (ipsptyp == 1) THEN
        CALL nonlocalc(psi1(:, N), Apsi(:,ithr), 0)  ! Apsi is the term Vnonloc *psi
      END IF

      ! add local potential (here the DFT potential for goal-density)
      DO ii = 1, kdfull2
        Apsi(ii,ithr) = Apsi(ii,ithr) + aloc(ishift + ii)*psi1(ii, N)
      END DO
!      epot(N) = 0D0
!      DO ii = 1, kdfull2
!        epot(N) = epot(N) + conjg(psi1(ii, N))*Apsi(ii,ithr)*dvol
!      END DO
      epot(N) = SUM(conjg(psi1(:, N))*Apsi(:,ithr))*dvol
      EspMod(N) = EspMod(N) + epot(N)    ! EspMod is sp energy,
    END DO
#if(dynomp)
!$OMP END PARALLEL DO
#endif
    FLUSH(6)
    IF ((mod(j, 1) .eq. 0)) THEN
      WRITE (*, '(a,i5,(200f10.6))') 'j,EspMod (redone) ',j,(EspMod(N), N=1, nstate)
    END IF
  END IF
!    DEALLOCATE(hpsip)
!    DEALLOCATE(vpsi)
    DEALLOCATE(Apsi)
!    DEALLOCATE(psik)
!    DEALLOCATE(hpsi, HpsiN)
!    DEALLOCATE(scal, psi1new)
    IF(ttimerta) THEN
      CALL system_clock(itimerta(10))
      WRITE(*,'(9(a,f9.4/))') &
                'PSI1 Apot:    ',DBLE(itimerta(11)-itimerta(6))/iraterta,&
                'PSI1 Loop 1: ',DBLE(itimerta(7)-itimerta(11))/iraterta,&
                'PSI1 cschmid: ',DBLE(itimerta(8)-itimerta(7))/iraterta,&
                'PSI1 epilogue:',DBLE(itimerta(10)-itimerta(8))/iraterta
    END IF
  CONTAINS

!
! Computes the action of the DCMF potential on the N-th. s.p  wavefunction.
! Output is
!   Apsi  = mean-field potential acting on psi
!           (before final composition, 'Apsi' is used as workspace)
!   ALpsi = full DCMF potential acting on psi
! all in coordinate space.
!
    SUBROUTINE calc_hamiltonien(Nact)
      INTEGER,INTENT(IN) :: Nact
!      COMPLEX(DP), ALLOCATABLE :: A(:), grad(:), gradApsi(:), jgradpsi(:)
      INTEGER :: ishift2,ii3,ithr2
      COMPLEX(DP) :: testdata(9)

      ! auxiliary arrays
!      ALLOCATE (A(kdfull2),grad(kdfull2),gradApsi(kdfull2),jgradpsi(kdfull2))
!      ALLOCATE (gradApsi(kdfull2))
      ishift2 = (ispin(Nact) - 1)*kdfull2
#if(dynomp)
      ithr2 = OMP_GET_THREAD_NUM()
      IF(ttestprint) WRITE(*,*) 'in Hamiltonien: ithr2=',ithr2,&
                      OMP_GET_MAX_THREADS(),OMP_GET_NUM_THREADS()
#else
      ithr2 = 0
#endif
!      WRITE(*,*) 'inside Hamiltonien:',ithr2,ishift,Nact,ispin(Nact),kdfull2
!      FLUSH(6)

      ! compute lambda^(eff)*psi (in real space)
!      ALpsi(:,ithr2) = Apot(:, 0) * psi1(:, Nact)

!      testdata(1) = wfovlp(psi1(:,Nact),vpsi(:,ithr2))
!      testdata(2) = wfovlp(psi1(:,Nact),ALpsi(:,ithr2))

      IF(Nact==1 .AND. ttimerta) CALL system_clock(itimerta(8))

!
!     the momentum dependent constraint
!     computes (muj (j-j0)+lambda) nabla phi+nabla((muj (j-j0)+lambda) phi)
!     uses 'Apsi' as workspace
      ! compute (muj (j-j0)+lambda) nabla phi+nabla((muj (j-j0)+lambda) phi)
      CALL xgradient_rspace(psi1(1, Nact), Apsi(:,ithr2))
      ALpsi(:,ithr2) = Apot(:,1)* Apsi(:,ithr2)
      Apsi(:,ithr2) = Apot(:,1)*psi1(:, Nact)
      CALL xgradient_rspace(Apsi(:,ithr2), Apsi(:,ithr2))
      ALpsi(:,ithr2) = ALpsi(:,ithr2) + Apsi(:,ithr2)

      CALL ygradient_rspace(psi1(1, Nact),  Apsi(:,ithr2))
      ALpsi(:,ithr2) = ALpsi(:,ithr2) + Apot(:,2)* Apsi(:,ithr2)
      Apsi(:,ithr2) = Apot(:,2)*psi1(:, Nact)
      CALL ygradient_rspace(Apsi(:,ithr2), Apsi(:,ithr2))
      ALpsi(:,ithr2) = ALpsi(:,ithr2) + Apsi(:,ithr2)

      CALL zgradient_rspace(psi1(1, Nact),  Apsi(:,ithr2))
      ALpsi(:,ithr2) = ALpsi(:,ithr2) + Apot(:,3)* Apsi(:,ithr2)
      Apsi(:,ithr2) = Apot(:,3)*psi1(:, Nact)
      CALL zgradient_rspace(Apsi(:,ithr2), Apsi(:,ithr2))
      ALpsi(:,ithr2) = ALpsi(:,ithr2) + Apsi(:,ithr2)
      ALpsi(:,ithr2) = Apot(:, 0) * psi1(:, Nact) - hbarm*ALpsi(:,ithr2)
!      testdata(3) = wfovlp(psi1(:,Nact),ALpsi(:,ithr2))
!      testdata(4) = wfovlp(ALpsi(:,ithr2),ALpsi(:,ithr2))



      IF(Nact==1 .AND. ttimerta) CALL system_clock(itimerta(7))

!
!     compute action of Kohn-Sham potential on wf 'psi'
!
      IF (ipsptyp == 1) THEN
        CALL nonlocalc(psi1(1, Nact), Apsi(:,ithr2), 0)
        DO ii3 = 1, kdfull2
          Apsi(ii3,ithr2) = Apsi(ii3,ithr2)+aloc(ishift2 + ii3)*psi1(ii3, Nact)
        END DO
      ELSE
        DO ii3 = 1, kdfull2
          Apsi(ii3,ithr2) = aloc(ishift2 + ii3)*psi1(ii3, Nact)
        END DO
      END IF


!      WRITE(*,*) &
!        'Hamilt: N,tests:',Nact,REAL(testdata(3:4),DP)
!      WRITE(*,'(a,3i4,6(1pg13.5))') &
!        'Hamilt: j,N,ithr,tests:',j,Nact,ithr2,testdata(1:3)
      IF(Nact==1 .AND. ttimerta .AND. iprintcount<20) THEN
        CALL system_clock(itimerta(9))
        WRITE(*,'(9(a,f9.4/))') &
                 'Ham 1:',DBLE(itimerta(7)-itimerta(8))/iraterta,&
                 'Ham 2:',DBLE(itimerta(9)-itimerta(7))/iraterta
        iprintcount = iprintcount + 1
      END IF

!      DEALLOCATE(gradApsi)

    END SUBROUTINE calc_hamiltonien
  END SUBROUTINE Calc_psi1

!____________________________________calcrhospin______________________________
  SUBROUTINE calcrhotot(rho, q0)
! total density 'rho' for COMPLEX wavefunctions 'q0'

    USE params
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: q0(kdfull2, kstate) ! wavefunctions
    REAL(DP), INTENT(OUT) :: rho(kdfull2)          ! total density

    INTEGER nb, ishift, ind

    rho = 0D0
    DO nb = 1, nstate
      DO ind = 1, nxyz
        rho(ind) = rho(ind) + occup(nb)* &
                   REAL((CONJG(q0(ind, nb)))*q0(ind, nb),DP)
      END DO
    END DO
    RETURN
  END SUBROUTINE calcrhotot




  SUBROUTINE fermi_total(ekmod, ereftot, occup, T0i, T1i, t2, chempot)
! Modifies occupations to reproduce the total sp energy 'eref' and
! the total particle number. This is here done for each spin separately.
    USE params, ONLY: DP, kstate, nstate
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ekmod(kstate)   ! s.p. energies
    REAL(DP), INTENT(IN) :: Ereftot         ! reference energy, to be reached
    REAL(DP), INTENT(IN) :: T0i,T1i         ! upper/lower bounds of temperature
    REAL(DP), INTENT(IN OUT) :: occup(kstate) ! final Fermi occupation
    REAL(DP), INTENT(IN OUT) :: chempot            ! Fermi energy
    REAL(DP), INTENT(OUT) :: t2               ! final temperature

    REAL(DP) :: et(2), T0, T1
    REAL(DP) :: occref(2)       ! reference particle number per spin
    REAL(DP) :: occtot(2)       ! actual particle number per spin
!    REAL(DP) :: T(2)            !
    INTEGER :: nspup, nspdw
    LOGICAL,PARAMETER :: ttest=.FALSE.

    nspup = eqstnspup
    nspdw = eqstnspdw
    occref(1) = SUM(occup(1:nspup))!total occupation spin1
    occref(2) = SUM(occup(nspup + 1:nstate))

    IF(ttest) WRITE (*, *) 'Fermi_Total entered:', nspup, nspdw, &
                            occref(1:2), T0i, T1i

    ! T0, T1 are lower and upper boundaries of temp
    T0 = T0i
    T1 = T1i

    !matching temperature by bracketing
    DO WHILE (abs(T1 - T0) > 1D-10)
      T2 = (T1 + T0)/2D0
      CALL occT1(occref(1), ekmod(1:nspup), et(1), chempot, occtot(1), nspup, &
                 T2, occup(1:nspup))
      CALL occT1(occref(2), ekmod(nspup + 1:nstate), et(2), chempot, occtot(2), &
                 nspdw, T2, occup(nspup + 1:nstate))
      IF (SUM(et) > ereftot) THEN ! set new bounds in T to embrace solution
        T1 = T2
      ELSE
        T0 = T2
      END IF
      IF(ttest) WRITE (*,'(a,8(1pg13.5))') 'T loop:', &
          t0, t1, ereftot, SUM(et), occref(1:2), occtot(1:2)
    END DO
    IF(ttest) WRITE (*, *) 'FERMI_TOTAL: T,chempot,part.nr.=', T2, chempot, SUM(occtot(:))
  END SUBROUTINE fermi_total

!_________________________________OccT1_______________________________________
  SUBROUTINE OccT1(occrefloc, enerloc, Etotloc, chempot, occtotloc, n, T, &
                   occuploc)
! Determines chemical potential 'chempot' and Fermi occpuation at given
! temperature to reproduce wanted particle number.
! This version assumes the same spin for all states.
!    occrefloc: total particle number to be reached
!    enerloc:   sp energy per mode
!    Etotloc:   resulting total s.p. energy
!    chempot:     Fermi energy
!    occtotloc: target total particle number
!    T:         temperature
!    occuploc:  return the resulting occupation numbers

    USE params, ONLY: DP, kstate, nstate
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: n
    REAL(DP),INTENT(IN) :: occrefloc, enerloc(n), T
    REAL(DP),INTENT(OUT) :: Etotloc, chempot, occuploc(n), occtotloc

    INTEGER, PARAMETER :: maxit = 50
    REAL(DP) :: chempot0, chempot1, fac, occ0, occ1
    INTEGER :: i, iter
    INTEGER :: orderloc(n)
    LOGICAL,PARAMETER :: ttest=.FALSE.

    chempot0 = minval(enerloc)-2D0*T
    chempot1 = maxval(enerloc)

    IF(ttest) WRITE(*,'(a,200(1pg13.5))') 'OCCT1 energies:',enerloc(1:n)
    occ0 = SUM(1D0/(1D0 + exp((enerloc(1:n) - chempot0)/T)))
    occ1 = SUM(1D0/(1D0 + exp((enerloc(1:n) - chempot1)/T)))
    IF(ttest) WRITE(*,*) occ0,occ1,chempot0,chempot1
    IF(occ0>occrefloc) STOP 'OCCT1 start: lower Fermi energy too high'
    IF(occ1<occrefloc) STOP 'OCCT1 start: upper Fermi energy too low'

    DO iter = 1, maxit
      IF (abs(chempot0 - chempot1) < 1D-10) EXIT
      chempot = (chempot0 + chempot1)/2D0
      occtotloc = 0D0
      Etotloc = 0D0
      DO i = 1, n
        occuploc(i) = 1D0/(1D0 + exp((enerloc(i) - chempot)/T))
        occtotloc = occtotloc + occuploc(i)
        Etotloc = Etotloc + occuploc(i)*enerloc(i)
      END DO
      IF(ttest) WRITE(*,'(a,i4,6(1pg12.5))') 'OCCT1:',&
          iter,occ0,occ1,occtotloc,chempot0,chempot1,chempot
      IF (occtotloc < occrefloc) THEN
        chempot0 = chempot
        occ0 = occtotloc
      ELSE
        chempot1 = chempot
        occ1 = occtotloc
      END IF
      IF((occ0-occrefloc)*(occ1-occrefloc)>0D0) &
           STOP 'OCCT1: particle number not embraced'
    END DO

  END SUBROUTINE occT1

  SUBROUTINE fermi_init(ekmod, T, occup, ispinact)
!
! Initializes chemical potential for a Fermi distribution of
! given temperature and spin.
! This routine is similar to 'OccT1', but can be used for arrays
! with mixed spin.
! I/O parameters:
!  ekmod    = s.p. energies (spin field 'ispin' communicated via 'params')
!  T        = given temperature
!  occup    = initially given occupation and output distribution
!  ispinact = actual spin

    USE params, ONLY: DP, kdfull2, kstate, nstate, ispin
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ekmod(nstate), T
    REAL(DP), INTENT(IN OUT) :: occup(kstate)
    INTEGER, INTENT(IN) :: ispinact
    REAL(DP)::mu1, mu2, mu3, o1, o2, o3, Oref !,Occ
    INTEGER::i

!   compute total particle number to be reached
    Oref = 0D0
    DO i = 1, nstate
      IF (ispin(i) .eq. ispinact) Oref = Oref + occup(i)
    END DO

!   find chemical checmical potential and final distrobution by bisection
    mu1 = minval(ekmod)-2D0*T
    mu2 = maxval(ekmod)
    o1 = occ(mu1, ekmod, T, occup, Oref, ispinact)
    o2 = occ(mu2, ekmod, T, occup, Oref, ispinact)
    DO WHILE (abs(mu1 - mu2) > 1D-5)
      mu3 = (mu1 + mu2)/2D0
      o3 = occ(mu3, ekmod, T, occup, Oref, ispinact)
      IF (o3 < 0D0) THEN
        mu1 = mu3
        o1 = o3
      ELSE
        mu2 = mu3
        o2 = o3
      END IF
    END DO
    DO i = 1, nstate
      IF (ispin(i) == ispinact) occup(i) = 1D0/(1D0 + exp((ekmod(i) - mu3)/T))
    END DO
  CONTAINS
    REAL(DP) FUNCTION occ(mu, ekmod, T, occup, O, ispinact)
      USE params, ONLY: DP, nstate, ispin, kstate
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: ekmod(nstate), O
      REAL(DP), INTENT(IN) :: occup(kstate)
      REAL(DP)::mu, T
      INTEGER, INTENT(IN)::ispinact
      REAL(DP)::mu1, mu2, mu3, o1, o2, o3
      INTEGER::i
      Occ = 0D0
      DO i = 1, nstate
        IF (ispin(i) == ispinact) Occ = Occ + 1D0/(1D0 + exp((ekmod(i) - mu)/T))
      END DO
      Occ = Occ - O
    END FUNCTION Occ
  END SUBROUTINE fermi_init


  SUBROUTINE srhomat(psi, psiorth, occuporth)
!
! Transform given state to natural orbitals.
!  psi      = given set of s.p. states (occupations come via 'params')
!  psiorth  = emerging natural orbitals
!  occuorth = corresponding new occcupations

    USE params, ONLY: DP, kdfull2, kstate, nstate, ispin, nrel2abs, nxyz, occup, apnum
    IMPLICIT NONE

!    REAL(DP), INTENT(IN) :: aloc(2*kdfull2)
    COMPLEX(DP), INTENT(IN) :: psi(kdfull2, kstate)
    COMPLEX(DP), INTENT(OUT) :: psiorth(kdfull2, kstate)
    REAL(DP), INTENT(OUT) :: occuporth(kstate)

    INTEGER:: i, j
    COMPLEX(DP) :: scal(nstate, nstate), scal0(nstate, nstate)
    COMPLEX(DP) :: D1on_rho(nstate, nstate), rhomat(nstate, nstate)
    REAL(DP) :: eigen(nstate), eigen0(nstate) ! eigenvalues
    COMPLEX(8) :: v(nstate, nstate), vp(nstate, nstate)
    COMPLEX(8) :: u(nstate, nstate), w(nstate, nstate)
    COMPLEX(8) :: winv(nstate, nstate)  !inverse of w
    REAL(DP) :: occup1(nstate)
    REAL(DP) :: devorthnorm

!   overlap matrix <psi_i |psi_j>
    CALL OverlapMatrix1(psi, psi, scal, ispin, 'matrix <psi|psi>')

!   check ortho-normality
    devorthnorm = 0D0
    DO i=1,nstate
      devorthnorm = devorthnorm + &
                    SUM(scal(:,i)**2)-scal(i,i)**2+(1D0-scal(i,i))**2
    END DO
    WRITE(*,*) 'SRHOMAT: deviation from ortho-normality=',devorthnorm

!   compute the density matrix 'rhomat' in the basis psi
    rhomat = 0D0
    DO i = 1, nstate
      rhomat(i, i) = occup(i)
    END DO
    scal0 = matmul(transpose(conjg(scal)), rhomat)
    rhomat = matmul(scal0, scal)
    WRITE (*, '(a,200(1pg13.5))') &
         'SRHOMAT: occup= ', (occup(i), i=1, nstate)

!   diagonalize overlap matrix
    CALL cdiagmat(scal, eigen0, v, nstate)
    WRITE (*, '(a,200(1pg13.5))') &
         'SRHOMAT: eigen0=', (eigen0(i), i=1, nstate)

!   diagonal matrices of inverse sqrt of eigenvalues
    D1on_rho = 0D0
    DO i = 1, nstate
      D1on_rho(i, i) = 1/sqrt(eigen0(i))
    END DO
!   scale eigenvector 'v' into 'vp' to transform to a unity matrix
    vp = matmul(v, D1on_rho)

!   compute transformed density matrix Vp^T*rhomat*Vp
    scal0 = matmul(transpose(conjg(Vp)), rhomat)
    scal0 = matmul(scal0, vp)
!   diagonalize Vp^T*rhomat*Vp
    CALL cdiagmat(scal0, eigen, u, nstate)
!   compute the composed transformation matrix
    w = matmul(vp, u)
    winv = matmul(transpose(conjg(u)), transpose(conjg(v)))

!   evaluate new basis states
    psiorth = 0D0
    DO j = 1, nstate
      DO i = 1, nstate
        psiorth(:, j) = psiorth(:, j) + w(i, j)*psi(:, i)
      END DO
    END DO
!   update the change of basis matrix
!    psitophi = matmul(winv, psitophi)

    WRITE (*, '(a,200(1pg13.5))') &
        'SRHOMAT: new occups=', (eigen(i), i=1, nstate)
    occuporth(1:nstate) = eigen

  END SUBROUTINE srhomat


  SUBROUTINE OverlapMatrix1(tab1, tab2, scal, ispin, mess)
!
! Overlap matrix between two (possibly different) wavefunction arrays.
! This version checks spin of wavefunctions and runs of the basic
! number of states.

    USE params, ONLY: DP, kdfull2, nstate, kstate, nrel2abs

    USE util, ONLY: wfovlp

    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: tab1(kdfull2, kstate) ! first set of wfs
    COMPLEX(DP), INTENT(IN) :: tab2(kdfull2, kstate) ! second set of wfs
    COMPLEX(DP), INTENT(OUT) :: scal(nstate, nstate) ! overlap matrix
    INTEGER, INTENT(IN) :: ispin(nstate)  ! spin assignment
    CHARACTER(*),INTENT(IN) :: mess       ! message for printing


    INTEGER:: i, j

#if(omp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,i) SCHEDULE(STATIC)
#endif
    DO j = 1, nstate
      DO i = 1, nstate
        IF (ispin(nrel2abs(i)) == ispin(nrel2abs(j))) THEN
          scal(i, j) = wfovlp(tab1(:, i), tab2(:, j))
        ELSE
          scal(i, j) = 0D0
        END IF
      END DO
    END DO
#if(omp)
!$OMP END PARALLEL DO
#endif

    WRITE (*, *) 'overlap matrix of '//mess//':'
    DO i = 1, nstate
      WRITE (*, '(200(1pg13.5)))') (scal(i, j), j=1, nstate)
    END DO

  END SUBROUTINE OverlapMatrix1


  SUBROUTINE OverlapMatrix2(tab1, tab2, scal, mess, nact)
!
! Overlap matrix between two (possibly different) wavefunction arrays.
! This version assuems that both sets of wavefunction have all the
! same spin.

    USE params, ONLY: DP, kdfull2, nstate, kstate

    USE util, ONLY: wfovlp

    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: tab1(kdfull2, kstate) ! first set of wfs
    COMPLEX(DP), INTENT(IN) :: tab2(kdfull2, kstate) ! second set of wfs
    COMPLEX(DP), INTENT(OUT) :: scal(nstate, nstate) ! overlap matrix
    INTEGER, INTENT(IN) :: nact  ! number of states in the wf-fields
    CHARACTER(*),INTENT(IN) :: mess       ! message for printing


    INTEGER:: i, j

#if(omp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,i) SCHEDULE(STATIC)
#endif
    DO j = 1, nact
      DO i = 1, nact
          scal(i, j) = wfovlp(tab1(:, i), tab2(:, j))
      END DO
    END DO
#if(omp)
!$OMP END PARALLEL DO
#endif

    WRITE (*, *) 'Overlap matrix of '//mess//':'
    DO i = 1, nact
      WRITE (*, '(200(1pg13.5)))') (scal(i, j), j=1, nact)
    END DO

  END SUBROUTINE OverlapMatrix2

!_______________________________cdiagspin_____________________________________
  SUBROUTINE cdiagspin(mat, eigen, vect, N)
!
! Driver routine to diagonalize the complex symmettric matrix in
! contiguous block of same spin.
!  mat   = complex symmetry matrix to be diagonalised.
!  eigen = real array of  eigenvalues
!  vect  = complex array of eigenvectors
!  N     = dimension of matrices and vectors

    USE params, ONLY: DP
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    COMPLEX(DP), INTENT(IN) :: mat(N, N)
    REAL(DP), INTENT(OUT) :: eigen(N)
    COMPLEX(DP), INTENT(OUT) :: Vect(N, N)

    INTEGER :: nspup, nspdw, i, j

    nspup = eqstnspup   ! nr. spin up states
    nspdw = eqstnspdw   ! nr. spin down states
    WRITE (*, *) 'CDIAGSPIN: nspup,nspdw=', nspup, nspdw

    vect = 0D0
    eigen = 0D0
    CALL CDIAG(mat(1:nspup, 1:nspup), eigen(1:nspup), &
               vect(1:nspup, 1:nspup), nspup, nspup)
    CALL CDIAG(mat(nspup+1:N, nspup+1:N), eigen(nspup+1:N), &
               vect(nspup+1:N, nspup+1:N), nspdw, nspdw)

  END SUBROUTINE cdiagspin


  SUBROUTINE indexx(n, arrin, indx)

! Sorting of an array 'arrin' by an index field 'indx' without
! changing the array 'arrin'.
! The output array 'indx' renders 'arrin(indx(j))' in ascending
! order for j = 1,2,...,n.
! This routine is adapted from Numerical Recipes.

    INTEGER :: i, n, indx(n), indxt, ir, j, l
    REAL(DP) :: arrin(n), q

    DO j = 1, n
      indx(j) = j
    END DO

    IF (n .eq. 1) RETURN

    l = n/2 + 1
    ir = n
10  CONTINUE
    IF (l .gt. 1) THEN
      l = l - 1
      indxt = indx(l)
      q = arrin(indxt)
    ELSE
      indxt = indx(ir)
      q = arrin(indxt)
      indx(ir) = indx(1)
      ir = ir - 1
      IF (ir .eq. 1) THEN
        indx(1) = indxt
        RETURN
      END IF
    END IF
    i = l
    j = l + l
20  CONTINUE
    IF (j .le. ir) THEN
      IF (j .lt. ir) THEN
        IF (arrin(indx(j)) .lt. arrin(indx(j + 1))) j = j + 1
      END IF
      IF (q .lt. arrin(indx(j))) THEN
        indx(i) = indx(j)
        i = j
        j = j + j
      ELSE
        j = ir + 1
      END IF
      GO TO 20
    END IF
    indx(i) = indxt
    GO TO 10

  END SUBROUTINE indexx

  SUBROUTINE PartNumPerSpin(mess, PartNum)
!
! Total number of electrons for given spin:
!   mess     = message to be printed
!   PartNum(1:2) = particle numbers spin up/down

    USE params, ONLY: DP, nstate, ispin, occup
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: mess
    REAL(DP), INTENT(OUT) :: PartNum(2)

    INTEGER::i

    PartNum = 0D0
    DO i = 1, nstate
      PartNum(ispin(i)) = PartNum(ispin(i)) + occup(i)
    END DO

    WRITE (*, '(a,a,2f12.8)') 'PartNumPerSpin', mess, PartNum(:)

  END SUBROUTINE PartNumPerSpin


  SUBROUTINE CorrectEnergy2(Wref, Eref, w, E, Wout, nloc)
!
! Corrects inital occupation 'w' to improved repproduction of
! particle number and total s.p. energy.
! This is done by a perturbative treatment as if input and output
! were Fermi distributions.

    USE params, ONLY: DP
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nloc ! number of states
    REAL(DP),INTENT(IN) :: Wref   ! reference value for total particle number
    REAL(DP),INTENT(IN) :: Eref   ! reference value for total s.p.energy
    REAL(DP),INTENT(IN) :: w(nloc)  ! array of given s.p. occupations
    REAL(DP),INTENT(IN) :: E(nloc)  ! array of s.p. energies
    REAL(DP),INTENT(OUT) :: Wout(nloc) ! corrected s.p. occupations


    REAL(DP):: Wloc, Eloc, dE, dW, fi(Nloc), A0, A1, A2, D, lam, mu
    INTEGER i

    Wloc = SUM(W)            !  total particle number from input occupations
    Eloc = dot_product(W, E) ! total s.p. energy from input arrays
    dw = Wref - Wloc         ! initial deviation particle number
    dE = Eref - Eloc         ! initial deviation sum of s.p. energies

    A0 = 0D0
    A1 = 0D0
    A2 = 0D0
    fi = 0D0
    DO i = 1, nloc
      IF (w(i) == 0.5D0) THEN
        fi(i) = 1D20
      ELSE IF (w(i) > 1D-20 .and. w(i) < (1D0 - 1D-15)) THEN
        fi(i) = 1/(LOG(w(i)) - LOG(1 - w(i)))**2
      ELSE
        fi(i) = 0D0
      END IF
      A0 = A0 + fi(i)
      A1 = A1 + fi(i)*E(i)
      A2 = A2 + fi(i)*E(i)**2
    END DO
    D = A0*A2 - A1*A1
    lam = (+dE*A1 - dW*A2)/D
    mu = (dW*A1 - dE*A0)/D
    DO i = 1, nloc
      Wout(i) = W(i) - lam*fi(i) - mu*fi(i)*E(i)
    END DO

    WRITE (*, '(a,10(1pg15.7))') &
     'CorrectEnergy2: Wreference,Ereference,Wachieved,Eachieved=',&
      Wref, Eref, SUM(Wout), dot_product(Wout, E)
    WRITE (*, '(a,200(1pg13.5))') 'E array:', E

  END SUBROUTINE CorrectEnergy2

!______________________________Temperature________________________________
  SUBROUTINE temperature(mu, T)
!
! Compute sort of a temperature by fitting a Fermi distribution to
! the occcupation distrobution in 'occup'.

    USE params, ONLY: DP, occup, amoy, nstate
    IMPLICIT NONE

    REAL(DP),INTENT(OUT) :: mu, T

    REAL(DP) :: TOL
    INTEGER INFO, lwa, i
    REAL(DP), ALLOCATABLE :: x(:), fvec(:), wa(:)
    INTEGER, ALLOCATABLE :: iwa(:)
    INTEGER, PARAMETER :: n = 2  ! number of independent variables, here(mu,t)

    lwa = nstate*n + 5*n + nstate       ! length of workspaces
    ALLOCATE (iwa(nstate), wa(lwa), x(n), fvec(nstate))

    mu = (amoy(1) + amoy(nstate))/2.d0   ! first guess
    T = 0.1D0
    x(1) = mu
    x(2) = T
    TOL = 1.d-5

    WRITE (*, *) 'before LMDIF1:', nstate, x(1), x(2), mu, T
!   LMDIF1 minimizes SUM(ff**2), ff=array of differences to Fermi distr.
    CALL lmdif1(ff, nstate, 2, x, fvec, tol, info, iwa, wa, lwa)

    mu = x(1)
    T = x(2)
    WRITE(*,*) 'after TEMPERATURE: s.p. energies, occcup_input, occcup_fit'
    DO i = 1, nstate
      WRITE (*,'(i4,3(1pg13.5))') i,amoy(i), occup(i), fvec(i) + occup(i)
    END DO

    DEALLOCATE (iwa, wa, x, fvec)

    RETURN
!________________________________ff__________________________________________
  CONTAINS
    SUBROUTINE ff(m, n, X, FVEC, IFLAG)

!   Difference between Fermi distribution for (mu,T) and given occupation.
!   As subroutine to be usable for the fitting routine LMDIF1.

      USE params, ONLY: DP, occup, amoy

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m, n, iflag
      REAL(DP), INTENT(IN) :: X(n)
      REAL(DP),INTENT(OUT) :: Fvec(m)

      INTEGER :: i

!      mu <--> X(1)
!      T  <--> X(2)

!      WRITE (*, *) 'in ff m,n,mu, T', m, n, mu, T
      DO i = 1, m
        fvec(i) = 1D0/(1D0 + exp((amoy(i) - X(1))/X(2))) - occup(i)
      END DO
    END SUBROUTINE ff

  END SUBROUTINE temperature

  REAL(DP) FUNCTION entropy(occupin)
!
! The entropy in the occcupation distribution 'occupin'.
! Number of states 'nstate' via module 'params'.

    USE params, ONLY: DP, nstate

    REAL(DP),INTENT(IN) :: occupin(1:nstate)

    REAL(DP) :: entropacc
    INTEGER :: i

    entropacc = 0D0
    DO i = 1, nstate
      IF ((occupin(i) > 1D-15) .and. (occupin(i) .lt. (1D0 - 1D-15))) THEN
        entropacc = entropacc - occupin(i)*LOG(occupin(i)) &
                              - (1D0-occupin(i))*LOG(1D0-occupin(i))
      END IF
    END DO

    entropy = entropacc

  END FUNCTION entropy


  SUBROUTINE OrderPsiOccEsp(n, psiin, occin, Esp)
!
! Physically reorder wavefunction and occcupation arrays in order
! of increasing s.p. energies.
!! Attention 1.: take care to order only contiguous arrays of same spin!
!! Attention 2.: coordination between states and s.p. energies is lost!

    USE params, ONLY: DP, kdfull2
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n          ! actual number of states to be sorted
    COMPLEX(DP), INTENT(IN OUT) :: psiin(kdfull2, n)  ! wfs to be sorted
    REAL(DP), INTENT(IN OUT) :: occin(n)              ! occups to be sorted
    REAL(DP), INTENT(IN OUT) :: Esp(n)                ! s.p. energies

    INTEGER:: indx(n), i, j
    COMPLEX(DP),ALLOCATABLE :: psiloc(:,:)
    REAL(DP),ALLOCATABLE :: occloc(:)
    LOGICAL :: tnew

    CALL indexx(n, Esp, indx)

! check need for reshuffling
    tnew=.FALSE.
    DO i = 1, n
      IF(indx(i).NE.i) THEN
        tnew = .TRUE.
        EXIT
      END IF
    END DO

    IF(tnew) THEN
      ALLOCATE(psiloc(kdfull2, n),occloc(n))

!
!   reorder wavefunctions
!
#if(omp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
#endif
      DO i = 1, n
        psiloc(:,i) = psiin(:, indx(i))
      END DO
#if(omp)
!$OMP END PARALLEL DO
#endif
#if(omp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(STATIC)
#endif
      DO i = 1, n
        psiin(:,i) = psiloc(:,i)
      END DO
#if(omp)
!$OMP END PARALLEL DO
#endif
!
!   reorder occupation numbers
!
      DO i = 1, n
        occloc(i) = occin(indx(i))
      END DO
      occin = occloc
!
!   reorder s.p. energies
!
      DO i = 1, n
        occloc(i) = Esp(indx(i))
      END DO
      Esp = occloc


      DEALLOCATE(psiloc,occloc)
    END IF

  END SUBROUTINE OrderPsiOccEsp





END MODULE rta_module
