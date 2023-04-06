
MODULE params
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: DP = KIND(1D0) ! PRECISION setting

!
! General settings for dimensions and basic arrays.
! Main communicator for dynamical variables, options, and various fields.
!

! number of nodes (=1 for serial version)
  INTEGER :: knode = 1
! max. nr. electron states per node
  INTEGER :: kstate = 20
! max. total nr. electron states
  INTEGER :: ksttot

! settings ad definitions for openmp parallel computing
#if(omp)
  INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS, OMP_GET_NUM_PROCS, OMP_NUM_THREADS
  INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  LOGICAL, EXTERNAL :: OMP_GET_DYNAMIC, OMP_GET_NESTED
!  EXTERNAL :: OMP_SET_NUM_THREADS, OMP_SET_DYNAMIC
  INTEGER :: numthr = 0 ! actual number of threads in openmp
  LOGICAL :: setdyn = .TRUE.
#else
  INTEGER, PARAMETER :: numthr = 1 ! actual number of threads in openmp
#endif
  INTEGER :: nthr ! iterator over num. of OMP threads. Starts counting at 0, so 1 less than the actual number of threads

! maximum number of ions
  INTEGER :: ng

  INTEGER, PARAMETER :: maxnang = 648 ! max. angular bins for PAD
  INTEGER, PARAMETER :: maxmps = 180 ! max. nr. analyzing points for PES

! physical units and constants (here Rydberg units)
  REAL(DP), PARAMETER :: e2 = 2.0D0
  REAL(DP), PARAMETER :: hbar = 1.0D0
  REAL(DP), PARAMETER :: ame = 0.5D0
  REAL(DP), PARAMETER :: h2m = hbar*hbar/2.0D0/ame ! h2m = hbar**2/2m

! frequently used mathematical constants
  REAL(DP), PARAMETER :: zero = 0.0D0
  REAL(DP), PARAMETER :: half = 0.5D0
  REAL(DP), PARAMETER :: one = 1.0D0
  REAL(DP), PARAMETER :: PI = 3.141592653589793D0
  REAL(DP), PARAMETER :: tPI = 6.28318530717958647692528676655900577D0
  REAL(DP), PARAMETER :: fourpi = 4.0D0*PI
  REAL(DP), PARAMETER :: small = 1.0D-20
  REAL(DP), PARAMETER :: sq2 = 1.4142135623730950D0
  REAL(DP), PARAMETER :: tsqrt2 = 2.82842712474619D0
  REAL(DP), PARAMETER :: sqrtpi = 1.772453850905516D0
  COMPLEX(DP), PARAMETER :: eye = (0D0, 1D0)

! maximum sizes of the box in x,y, and z
  INTEGER :: kxbox = 0, kybox = 0, kzbox = 0
! deduced grid parameters
  INTEGER :: kxmax, kymax, kzmax
  INTEGER :: nx2, ny2, nz2
  INTEGER :: nx2fine, ny2fine, nxyfine
  INTEGER :: nxyz, nyf, nxyf
  INTEGER :: kdfull2, kdfull2fine
  INTEGER :: nx, ny, nz
  INTEGER :: nxfine, nyfine, nzfine
  INTEGER :: knodem ! =knode-1
! bounds of loops
  INTEGER :: minx = 1, maxx
  INTEGER :: miny = 1, maxy
  INTEGER :: minz = 1, maxz
! bounds of esc. el.
  INTEGER :: nbnx = 2, nbxx
  INTEGER :: nbny = 2, nbxy
  INTEGER :: nbnz = 2, nbxz
! mid-point for x,y,z-values
  INTEGER :: nxsh = 0, nysh = 0, nzsh = 0

  REAL(DP) :: dx = 0D0, dy = 0D0, dz = 0D0, dvol ! mesh spacing, volume
  REAL(DP), ALLOCATABLE :: xval(:), yval(:), zval(:) ! grid coordinates
  REAL(DP), ALLOCATABLE :: xt2(:), yt2(:), zt2(:) ! coordinates**2

  REAL(DP), ALLOCATABLE :: enonlo(:) ! non-local s.p.energy
  REAL(DP), ALLOCATABLE :: amoy(:) ! single particle energies
  REAL(DP), ALLOCATABLE :: spvariance(:) ! s.p. energy variances
  REAL(DP), ALLOCATABLE :: spvariancebi(:) ! s.p. en.varian. boost-inv.
  REAL(DP), ALLOCATABLE :: spenergybi(:) ! s.p. energy boost-inv.
  REAL(DP), ALLOCATABLE :: spnorm(:) ! norm of s.p. wf
  REAL(DP), ALLOCATABLE :: occup(:) ! occupation weight
  INTEGER :: nstate = 1, nspdw = -9999 ! Nr. of states
  INTEGER :: nstate_all ! total Nr. of states (mpi)
  INTEGER, ALLOCATABLE :: nq(:, :) ! nodes of init. states
  INTEGER, ALLOCATABLE :: ipar(:, :) ! xyz-parities
  INTEGER, ALLOCATABLE :: ispin(:) ! spin of states
  INTEGER, ALLOCATABLE :: nrel2abs(:), nabs2rel(:) ! POINTER to wfs
  INTEGER, ALLOCATABLE :: nrel2abs_other(:, :) ! POINTER to wfs
  INTEGER, ALLOCATABLE :: nhome(:) ! home node of wf


! basic parameters

  LOGICAL :: directenergy = .false. ! to compute energy directly
  INTEGER :: numspin = 2 ! number of spin components
  REAL(DP) :: omeg, eferm ! initial Gaussians
  REAL(DP) :: xfac, yfac, zfac ! initial. auxiliary
  REAL(DP) :: epswf = 0.2D0, e0dmp = 2D0, epsorc = 1D-8 ! convergence
  REAL(DP) :: b2occ = 0.4D0, gamocc = 10D0, deocc = 0D0, osfac = 1D0 ! h.o.initalization
  REAL(DP) :: temp = 0D0, occmix = 0.5D0, epsoro = 1D-6 ! electron temperature
  REAL(DP) :: endcon = 1D-5, radjel = 4D0, surjel = 1D0 ! jellium params
  INTEGER :: itback = 200 ! jellium params
  REAL(DP) :: bbeta = 0D0, gamma = 0D0, beta4 = 0D0
  REAL(DP) :: falph, fbeta, fhexe ! jellium auxiliary
  REAL(DP) :: dpolx = 0D0, dpoly = 0D0, dpolz = 0D0 ! static dipole potential
  LOGICAL :: tdipolxyz ! switch to dipole fields
  INTEGER :: irotat = 0 ! rotational initialisation
  LOGICAL :: tspindip = .FALSE. ! spin dipole initialization
  REAL(DP) :: phirot = 0D0
  INTEGER :: nelect = -9999, nion = -9999, nion2 = 1
  REAL(DP) :: charge ! Nr. el. & ions
  REAL(DP) :: scaleion = 1D0, scaleionx = 1D0, scaleiony = 1D0
  REAL(DP) :: scaleionz = 1D0
  LOGICAL :: tmoms_rel_cm = .TRUE. ! .TRUE. = relative to c.m., .FALSE. = relative to box
  REAL(DP) :: shiftionx = 0D0, shiftiony = 0D0, shiftionz = 0D0
  REAL(DP) :: rotionx = 0D0, rotiony = 0D0, rotionz = 0D0
  LOGICAL :: tfixcmion = .FALSE.
  REAL(DP) :: shiftWFx = 0D0, shiftWFy = 0D0, shiftWFz = 0D0
  INTEGER :: ispinsep = 0
  REAL(DP) :: comx, comy, comz, dion(3), qtion(3, 3)
  REAL(DP) :: apnum, rmsion, dmdistion, rmsel, rhopss, rhomix
  REAL(DP) :: codx, cody, codz, del(3), qtel(3, 3)
  REAL(DP) :: time_absinit

  LOGICAL :: tshiftcmtoorigin = .FALSE., tfreezekspot = .FALSE.
  INTEGER :: jdensitydiff = 0, jdensity2d = 0, jdensitydiff2d = 0
  INTEGER :: nmptheta = 2, nmpphi = 1, nmps, jmp = 0, imps(maxmps)
  INTEGER :: jnorms = 0
  INTEGER :: jcharges = 0

  INTEGER, ALLOCATABLE :: ispin_target(:)
  REAL(DP), ALLOCATABLE :: occ_target(:)
  REAL(DP), ALLOCATABLE :: spe_target(:)
  COMPLEX(DP), ALLOCATABLE :: psi_target(:, :)
  INTEGER, ALLOCATABLE :: match(:, :)
  REAL(DP) :: totintegprob
  REAL(DP) :: aver_estar, emin_target, emax_target
  INTEGER :: nstate_target, nmatch

  LOGICAL :: tplotorbitals = .FALSE.
  INTEGER :: ievaluate = 0 ! ????
  REAL(DP) :: ekin0pp = 0D0, vxn0 = 0D0, vyn0 = 0D0, vzn0 = -1D0
  REAL(DP) :: rheatclust
  REAL(DP) :: reference_energy = 0D0
  REAL(DP) :: drcharges = 5D0

  LOGICAL :: iflocaliz = .FALSE. ! evaluate localization
  INTEGER :: ifile = 60 ! UNIT number for SAVE
  INTEGER :: myn = 0 ! nr. of actual node (zero in sequetiential code)
  INTEGER :: ismax = 1000 ! maximum number of static iterations
  INTEGER :: itmax = 1000 ! number of time steps for electronic propagation
  INTEGER :: idyniter = 0 ! number iterations to start dynamic E0DMP

! Printing informations

  CHARACTER(LEN=13) :: outnam ! output file name
  CHARACTER(LEN=13) :: outname ! output file name with node number in it.
  INTEGER :: istinf = 10 ! frequency for printing info from static iteration

  LOGICAL :: iffastpropag = .TRUE. ! a bit more efficient time step in TV splitting
  LOGICAL :: ifexpevol = .FALSE. ! switch to exponential evolution
  ! to be set at compile time

  INTEGER :: irest = 0 ! switch to restart/READ dynamics from file ’SAVE’
  INTEGER :: isaves = 0, isaved = 0 ! saves results after every ’isave..’ steps
  LOGICAL :: tstat = .FALSE. ! switch to restart/READ statics from file ’RSAVE’
  ! on file ’RSAVE’ in and after static iteration
  ! on file ’SAVE’ in dynamic propagation

  INTEGER :: jdensity1d = 0 ! print densities integrated over xy, xz, yz in 3 separate files

  INTEGER :: jpos = -9999, jvel = -9999, jesc = -9999, jforce = 0, jposcm = 0, jgeomion = 0
  INTEGER :: jinfo = 10, jdip = -9999, jdiporb = 0, jquad = 0, jang = 0, jspdp = 0, jenergy = 10
  INTEGER :: jgeomel = 0
  INTEGER :: jangabso = 0
  INTEGER :: jelf = 0
  INTEGER :: jstinf = 0
  INTEGER :: jstateoverlap = 0
  INTEGER :: nabsorb = 0, ifsicp = 2, ionmdtyp = 0, icooltyp = 0
  INTEGER :: ipsptyp = -999, idenfunc = 1, icooltimes = 0
  INTEGER :: jheatmod = 0 ! modulus for re-heating the system
  INTEGER :: itersicp6
  LOGICAL :: tstinf, init_ao = .FALSE.
  LOGICAL :: ifspemoms = .FALSE., iftransme = .FALSE., ifredmas = .FALSE.

  REAL(DP) :: rheattemp = 0D0 ! re-heat temperature
  REAL(DP) :: tempion = 0D0, dt1 = -1D99
  REAL(DP) :: centfx = 0D0, centfy = 0D0, centfz = 0D0
  REAL(DP) :: shiftinix = 0D0, shiftiniy = 0D0, shiftiniz = 0D0
  REAL(DP) :: systime_factor ! factor to convert results of 'systime'

  INTEGER :: ifhamdiag = 20 ! frequency for Hamiltonian diagonalization
  REAL(DP) :: variance_gain = 0.33333D0 ! variance criterion for diagonalization

! spatial fields as densities and potentials
  REAL(DP), ALLOCATABLE :: rhojel(:) ! jellium density
  REAL(DP), ALLOCATABLE :: potion(:) ! pseudopotentials
  REAL(DP), ALLOCATABLE :: potFixedIon(:) ! potential from frozen ions
  REAL(DP), ALLOCATABLE :: chpcoul(:) ! Coulomb potential

! fields and variables for analysis of electron emission

  REAL(DP), ALLOCATABLE :: rhoabso(:) ! storage for absorbed density
  REAL(DP), ALLOCATABLE :: spherMask(:) ! mask for spherical absorbing bounds
  REAL(DP), ALLOCATABLE :: spherloss(:) ! loss factor for spher. abs. bounds
  REAL(DP) :: bcrad, powabso = 0.0675D0 ! width & power of abs. bounds
  INTEGER :: ispherAbso = 1 ! switch to spherical abs. bounds
  INTEGER :: iangabso = 0, nangtheta = 1, nangphi = 1
  INTEGER :: indicesmp(maxnang*maxnang)
  REAL(DP) :: angthetah = PI, angthetal = 0D0
  REAL(DP) :: angphil = 0D0, angphih = 2*PI, delomega, xango, yango, zango
  LOGICAL, ALLOCATABLE :: tgridabso(:) ! array tagging absorbing points
  REAL(DP), ALLOCATABLE :: rhoabsoorb(:, :)


  INTEGER, PARAMETER :: kmom = 35
  INTEGER :: nrmom
  REAL(DP) :: qe(kmom), se(5), ajx, ajy, ajz
  REAL(DP), ALLOCATABLE :: qeorb_all(:, :)


! storage for base wavefunctions in case of dynamics with exact exchange
! or in CASE of computing overlaps with initial state
  COMPLEX(DP), ALLOCATABLE :: psisavex(:, :)

! the energy transmitted from calc-lda to info etc

  REAL(DP) :: enrear ! rearrangement energy
  REAL(DP) :: ecback, ecrho, ecorr, dt12, sgaus, ekion, energy
  REAL(DP) :: energ2, enerpw, encoulsp, entrop, epot, espnb, esh1
  REAL(DP) :: etot, ekionold, qold2, qold3, qold4
  REAL(DP) :: ekmat = 0D0, engg, enii, enig, ecrhoimage = 0D0
  REAL(DP), ALLOCATABLE :: ekinsp(:), evarsp(:), evarsp2(:), epotsp(:)
  INTEGER :: jekion, iquery4
  REAL(DP), ALLOCATABLE :: hmatrix(:, :)
  COMPLEX(DP), ALLOCATABLE :: expmatrix(:, :)
  REAL(DP) :: symcon
  INTEGER :: ndims(2)

! dynamic variables of ionic motion

  INTEGER, PARAMETER :: nxsg = 7, nysg = 7, nzsg = 7 ! size of subgrids
  INTEGER :: modionstep = 1 ! modulus for ion step
  INTEGER :: inewforce
  INTEGER, ALLOCATABLE :: nfixed(:) ! Nr. of fixed ions
  LOGICAL, ALLOCATABLE :: tblock(:)
  REAL(DP), ALLOCATABLE :: cx(:), cy(:), cz(:), cpx(:), cpy(:), cpz(:)
  REAL(DP), ALLOCATABLE :: dcx(:), dcy(:), dcz(:), dcpx(:), dcpy(:), dcpz(:)
  REAL(DP), ALLOCATABLE :: fx(:), fy(:), fz(:), flx(:), fly(:), flz(:)
  REAL(DP), ALLOCATABLE :: fprojx(:), fprojy(:), fprojz(:)

! book keeping for atomic orbital initialization
  REAL(DP), ALLOCATABLE :: radini(:)
  INTEGER, ALLOCATABLE :: initord(:, :), ipol(:)
  INTEGER, ALLOCATABLE :: nmaxst(:)

! parameters for simulated annealing
  REAL(DP), ALLOCATABLE :: eloc(:), enoloc(:), eion(:, :), eiinew(:)
  REAL(DP) :: cptemp, delpos, ERR, binerg, errtot, erfac1, trfac1, prfac1
  REAL(DP) :: trfac2, prfac2, errsim, eiontot, enoloctot, facann
  REAL(DP) :: eloctot, errks0, errks, sumvar, sumvar2
  INTEGER :: ionsin, nrun, nloop1, nloop2, loop1, nyes, ncon, ncsim, iknow
  INTEGER :: isize_seed
  INTEGER, ALLOCATABLE :: rand_seed(:)

! parameters for external excitations by laser or projectile
  INTEGER :: idenspl = 0 ! print densities integrated over x, y and z in 3 separated files
  ! (at the end in static, every 'idenspl' steps in dynamic)
  INTEGER :: i3dz = 0, i3dx = 0, i3dstate = 0
  INTEGER :: itft = 3
  REAL(DP) :: tnode = 0D0, deltat = 0D0, tpeak = 0D0, omega = 0D0, e0 = 0D0, time, tfs = 0D0
  REAL(DP) :: e1x = 1D0, e1y = 0D0, e1z = 1D0, phi = 0D0
  REAL(DP) :: e2x = 0D0, e2y = 0D0, e2z = 1D0, phase2 = 0D0, omega2 = 0D0, e0_2 = 0D0
  REAL(DP) :: tstart2 = 0D0, deltat2 = 0D0
  REAL(DP) :: fl(6), power
  INTEGER :: jlaser = 0, ilas = 0
  REAL(DP), TARGET :: datalaser(18)
  REAL(DP), POINTER :: ppower, elaser, fpulseinteg1, fpulseinteg2, fXUVinteg, &
                       ex1, ey1, ez1, ex2, ey2, ez2, exXUV, eyXUV, ezXUV, &
                       fpulse2integ1, fpulse2integ2, f2pulse2integ1, f2pulse2integ2
! datalaser summarizes ppower, elaser, fpulseinteg for 1st pulse, 2nd pulse,
! XUV attotrain ex, ey and ez of the field resulting from the 1st
! same for the 2nd (standard) pulse, same for the attopulse train
! fpulse2integ1, fpulse2integ2, f2pulse2integ1, f2pulse2integ2 --> for phase corrected PES
! fpulseinteg1, fpulseinteg2 ! integrated pulses for gauge transformation
  REAL(DP), TARGET :: dataold(7)
! dataold contains old values of time, foft1, foft2, foftXUV, acc1, acc2, accXUV
  REAL(DP), POINTER :: timeold, foft1old, foft2old, foftXUVold, &
                       acc1old, acc2old, accXUVold
  REAL(DP) :: projcharge = 0D0 ! projectile charge
  REAL(DP) :: projvelx = 0D0, projvely = 0D0, projvelz = 0D0 ! projectile velocity
  REAL(DP) :: projinix = 0D0, projiniy = 0D0, projiniz = 0D0 ! initial projectile POSITION
  ! impact PARAMETER = min(projinix,projiniy,projiniz)

! parameters for extended cases and observables

! workspace for communication
  REAL(DP) :: rvectmp2(3), rvectmp(3)
  INTEGER :: iindtmp(3)

! POINTER for energy-density functional
  PROCEDURE(), POINTER :: calc_lda

! these includes should be shifted to own modules

#include "pseudo.F90"

CONTAINS

  SUBROUTINE init_baseparams()

! Initializes those parameters which are derived from the basic
! parameters.

    nthr = numthr - 1
! deduced grid parameters
    kxmax = kxbox/2 + 1; kymax = kybox/2 + 1; kzmax = kzbox/2 + 1
    nx2 = kxbox; ny2 = kybox; nz2 = kzbox
    nxyz = nx2*ny2*nz2; nyf = nx2; nxyf = nx2*ny2
    kdfull2 = kxbox*kybox*kzbox
    nx = nx2/2; ny = ny2/2; nz = nz2/2

! deduce nr. of states per node, readjust total nr. of states
! note: the input variable 'kstate' means total nr. of states
! and is copied first to the correct variable 'KSTTOT'.
    ksttot = kstate
    kstate = ksttot/knode ! trial value
    IF ((knode*kstate) < ksttot) kstate = kstate + 1
    ksttot = knode*kstate ! final value

    knodem = knode - 1

! bounds of loops
    minx = 1; maxx = nx2
    miny = 1; maxy = ny2
    minz = 1; maxz = nz2
! bounds of esc. el. ???
    nbnx = 2; nbxx = nx2 - 1
    nbny = 2; nbxy = ny2 - 1
    nbnz = 2; nbxz = nz2 - 1
! mid-point for x,y,z-values
    nxsh = nx2/2; nysh = ny2/2; nzsh = nz2/2

! max. nr. of ions
    ng = nion

! max. nr. of ions
    ng = nion

! pointers to laser properties

    ppower => datalaser(1)
    elaser => datalaser(2)
    fpulseinteg1 => datalaser(3)
    fpulseinteg2 => datalaser(4)
    fXUVinteg => datalaser(5)
    ex1 => datalaser(6)
    ey1 => datalaser(7)
    ez1 => datalaser(8)
    ex2 => datalaser(9)
    ey2 => datalaser(10)
    ez2 => datalaser(11)
    exXUV => datalaser(12)
    eyXUV => datalaser(13)
    ezXUV => datalaser(14)
    fpulse2integ1 => datalaser(15)
    fpulse2integ2 => datalaser(16)
    f2pulse2integ1 => datalaser(17)
    f2pulse2integ2 => datalaser(18)

    timeold => dataold(1)
    foft1old => dataold(2)
    foft2old => dataold(3)
    foftXUVold => dataold(4)
    acc1old => dataold(5)
    acc2old => dataold(6)
    accXUVold => dataold(7)

  END SUBROUTINE init_baseparams

  SUBROUTINE init_fields()

! Allocates basic fields and initializes them if needed.
! Sets defaults for many input parameters.

    ALLOCATE (nq(3, ksttot)) ! nodes of init. states
    ALLOCATE (ipar(3, ksttot)) ! xyz-parities
    ALLOCATE (ispin(ksttot)) ! spin of states
    ALLOCATE (nrel2abs(ksttot), nabs2rel(ksttot)) ! pointer to wfs
    ALLOCATE (nrel2abs_other(kstate, 0:knode - 1)) ! pointer to wfs
    ALLOCATE (nhome(ksttot)) ! home node of wf

    ALLOCATE (amoy(kstate)) ! single particle energies
    amoy = 0D0
    ALLOCATE (enonlo(kstate)) ! single particle energies
    enonlo = 0D0
    ALLOCATE (spvariance(kstate)) ! s.p. energy variances
    ALLOCATE (spvariancebi(kstate)) ! s.p. energy variances
    ALLOCATE (qeorb_all(ksttot, 11)) ! s.p. dipole moments
    ALLOCATE (spenergybi(kstate))
    ALLOCATE (spnorm(kstate)) ! norm of s.p. wf
    ALLOCATE (occup(ksttot)) ! occupation weight
    ALLOCATE (rhojel(kdfull2)) ! jellium density
    ALLOCATE (potion(kdfull2)) ! pseudopotentials
    ALLOCATE (potFixedIon(kdfull2)) ! potential from frozen ions
    ALLOCATE (chpcoul(kdfull2)) ! Coulomb potential

    IF (nabsorb > 0) THEN
      ALLOCATE (rhoabso(kdfull2)) ! storage for absorbed density
      ALLOCATE (spherMask(kdfull2)) ! mask for spherical absorbing bounds
      ALLOCATE (spherloss(kdfull2)) ! loss factor for spher. abs. bounds
      ALLOCATE (tgridabso(kdfull2))
    END IF

    ALLOCATE (ekinsp(kstate), evarsp(kstate), evarsp2(kstate), epotsp(kstate))
    ekinsp = 0D0
    evarsp = 0D0
    evarsp2 = 0D0
    epotsp = 0D0

    IF (ifsicp > 7) ALLOCATE (hmatrix(kstate, kstate))

! fields for PsP projectors
    IF (ipsptyp /= 0) THEN
      ALLOCATE (ifin(0:ng), icount(knl, 0:ng))
      ALLOCATE (p0_1(knl, 0:ng), p0_2(knl, 0:ng), p1_1(knl, 0:ng), p1_1x(knl, 0:ng))
      ALLOCATE (p1_1y(knl, 0:ng), p1_1z(knl, 0:ng))
      ALLOCATE (p0_3(knl, 0:ng), p1_2(knl, 0:ng), p1_2x(knl, 0:ng))
      ALLOCATE (p1_2y(knl, 0:ng), p1_2z(knl, 0:ng))
      ALLOCATE (p2_1(knl, 0:ng), p2_xy(knl, 0:ng), p2_xz(knl, 0:ng))
      ALLOCATE (p2_yz(knl, 0:ng), p2_xy2(knl, 0:ng), p2_z2(knl, 0:ng))

      ALLOCATE (ifinfine(0:ng), icountfine(knl, 0:ng))
      ALLOCATE (icountfinesp(knl*ng))
      ALLOCATE (icountsp(knl*ng))
      ALLOCATE (p0_1fine(knl, 0:ng), p0_2fine(knl, 0:ng), p1_1fine(knl, 0:ng), p1_1xfine(knl, 0:ng))
      ALLOCATE (p1_1yfine(knl, 0:ng), p1_1zfine(knl, 0:ng))
      ALLOCATE (p0_3fine(knl, 0:ng), p1_2fine(knl, 0:ng), p1_2xfine(knl, 0:ng))
      ALLOCATE (p1_2yfine(knl, 0:ng), p1_2zfine(knl, 0:ng))
      ALLOCATE (p2_1fine(knl, 0:ng), p2_xyfine(knl, 0:ng), p2_xzfine(knl, 0:ng))
      ALLOCATE (p2_yzfine(knl, 0:ng), p2_xy2fine(knl, 0:ng), p2_z2fine(knl, 0:ng))
      ALLOCATE (tnonloc(0:ng))
    END IF

! fields for ionic variables
    ALLOCATE (tblock(0:ng))
    ALLOCATE (nfixed(ng), np(0:ng))
    ALLOCATE (cx(ng), cy(ng), cz(ng), cpx(ng), cpy(ng), cpz(ng))
    ALLOCATE (dcx(ng), dcy(ng), dcz(ng), dcpx(ng), dcpy(ng), dcpz(ng))
    ALLOCATE (fx(ng), fy(ng), fz(ng), flx(ng), fly(ng), flz(ng))
    ALLOCATE (fprojx(ng), fprojy(ng), fprojz(ng))

! book keeping for atpmic orbital initialization
    ALLOCATE (radini(1:ng), ipol(1:ng))
    ALLOCATE (initord(1:3, 1:ng), nmaxst(1:ng))

    IF (icooltyp == 3) ALLOCATE (eloc(0:ng), enoloc(0:ng), eion(ng, ng), eiinew(ng))


    RETURN

  END SUBROUTINE init_fields


END MODULE params
