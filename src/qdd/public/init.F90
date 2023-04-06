
!------------------------------------------------------------

SUBROUTINE initnamelists

! Sets defaults for input variables
! and THEN reads input variables through NAMELIST.
! Parameters are communicated via module 'params'.
! For meaning of the parameters see documentation.

!------------------------------------------------------------
  USE params
  USE kinetic
  USE coulsolv, ONLY: tcoulfalr
  USE RTA_module
  IMPLICIT NONE

  CHARACTER(LEN=80) :: title

!INTEGER,PARAMETER :: kparall=11
  INTEGER:: iu, nnx2, maxnum
  CHARACTER(LEN=3) :: num

  REAL(DP)::dx2

  NAMELIST /global/ nelect, nion, nspdw, numspin, nion2,&
#if(omp)
    numthr, &
#endif
    ipsptyp, &
    temp, b2occ, gamocc, deocc, osfac, &
    init_ao, kstate, kxbox, kybox, kzbox, dx, dy, dz, &
    endcon, itback, &
    dpolx, dpoly, dpolz, &
    tcoulfalr, idenfunc, ifsicp,&
    scaleion, scaleionx, scaleiony, scaleionz, &
    shiftionx, shiftiony, shiftionz, &
    rotionx, rotiony, rotionz, &
    tmoms_rel_cm, tshiftcmtoorigin, &
    shiftwfx, shiftwfy, shiftwfz, ispinsep

  NAMELIST /static/ ismax,  istinf, isaves, idyniter, ifhamdiag, &
    epswf, e0dmp, epsoro,  occmix,  variance_gain, tstat, &
    ifspemoms, iftransme,  iflocaliz, tplotorbitals

  NAMELIST /dynamic/ nabsorb, powabso, ispherabso,&
    itmax, isaved, dt1, irest, &
    centfx, centfy, centfz, shiftinix, shiftiniy, shiftiniz, &
    tspindip, &
    irotat, phirot, &
    iffastpropag, ifexpevol, &
    ionmdtyp, ifredmas, modionstep, icooltyp, &
    tempion, &
    jpos, jvel, jesc, jforce, jgeomion, &
    jdip, jdiporb, jquad, jang, jangabso, jspdp, jinfo, jenergy, &
    jposcm, jgeomel, jelf, jstinf, &
    jlaser, &
    itft, tnode, deltat, tpeak, omega, e0, &
    projcharge, projvelx, projvely, projvelz, &
    projinix, projiniy, projiniz, &
    e1x, e1y, e1z, e2x, e2y, e2z, phi, &
    phase2, omega2, e0_2, tstart2, deltat2, &
    iangabso, nangtheta, nangphi, &
    delomega, angthetal, angthetah, angphil, angphih, &
    tfreezekspot, tfixcmion, &
    nmptheta, nmpphi, jmp, &
    jnorms, jdensitydiff, jdensitydiff2d, &
    jdensity2d, jdensity1d, jcharges, drcharges, &
    jstateoverlap, &
    jrtaint

  NAMELIST /RTAparams/ rtamu, rtamuj, rtasumvar2max, rtaeps, rtae0dmp, &
    rtatempinit, rtasigee, rtars, &
    rtaferminit, rtaerr1, rtadT1, rtaExitErr, rtaExitErrj, itmaxDCMF
!   rtaforcetemperature



  WRITE (6, *) 'Reading for005.// ...'

! initialize the variables with DEFAULT values

  delomega = (angthetah - angthetal)/4D0/nangtheta
  delomega = MIN(delomega, (angphih - angphil)/4D0/nangphi)



! OPEN input files


  IF (myn == 0) THEN
    OPEN (UNIT=5, STATUS='old', FORM='FORMATTED', FILE='for005')
    iu = 7
    WRITE (*, *) ' enter title (=qualifier) for that run:'
    READ (5, *) title
    outnam = title(1:13)
    IF (outnam == ' ') STOP " no title given "
  outname = trim(outnam)

    OPEN (UNIT=7, STATUS='unknown', FILE='out_details.'//outname)
    WRITE (6, *) ' title is now: '//outnam
    WRITE (iu, *) ' title is now: '//outnam
    CLOSE (5)
    WRITE (6, *) outnam, "**************************************"
    OPEN (5, STATUS='old', FORM='FORMATTED', FILE='for005.'//outnam)

! READ input from namelists

    READ (5, global, END=99999)
    WRITE (*, *) ' GLOBAL READ'

    IF (dx == 0D0) THEN
      STOP ' DX must be given explicitely'
    ELSE
      IF(dy == 0D0) dy = dx
      IF(dz == 0D0) dz = dx
    END IF
    IF (dx < 0D0 .OR. dy < 0D0 .OR. dz < 0D0) &
      STOP 'any grid spacing must be positive'
    IF (kxbox == 0) THEN
      STOP ' KXBOX must be given explicitely'
    ELSE
      IF(kybox == 0) kybox = kxbox
      IF(kzbox == 0) kzbox = kxbox
    END IF

    READ (5, static, END=99999)
    WRITE (*, *) ' STATIC READ'

    READ (5, dynamic, END=99999)
    WRITE (*, *) ' DYNAMIC READ'

    IF (jrtaint > 0) READ (5, RTAparams, END=99999)
    WRITE (*, *) ' RTA READ'

99999 CLOSE (5)
    IF (itmax > 0) THEN
      IF (jdip < 0) STOP "you must specify JDIP in NAMELIST DYNAMIC"
      IF (jesc < 0) STOP "you must specify JESC in NAMELIST DYNAMIC"
    END IF

  END IF

  IF (ionmdtyp == 0) THEN
    jpos = 0
    jvel = 0
  ELSE
    IF (jpos < 0) STOP "you must specify JPOS in NAMELIST DYNAMIC"
    IF (jvel < 0) STOP "you must specify JVEL in NAMELIST DYNAMIC"
    IF (nion2 == 0) STOP "IONMDTYP==0 required for jellium"
  END IF
! adapt input parameters IF necessary
  IF (nion2 == 0) tmoms_rel_cm = .FALSE. ! relat. to center of box for jellium

  IF (e0 /= 0D0 .and. jlaser < 0) STOP 'you must specify JLASER in NAMELIST DYNAMIC'
  IF (e0 == 0D0) jlaser = 0

  tdipolxyz = dpolx*dpolx + dpoly*dpoly + dpolz*dpolz .GT. 0D0

  RETURN
END SUBROUTINE initnamelists

!-----changePerio--------------------------------------------

SUBROUTINE changeperio

! Reads pseudo potential parameters from NAMELIST thus over-riding
! DEFAULT settings.

  USE params
  IMPLICIT NONE

  NAMELIST /perio/  ch, amu, cc1, cc2, crloc, &
    h0_11g, h0_12g, h0_22g, h0_33g, h1_11g, h1_22g, h2_11g, &
    dr1, dr2, prho1, prho2, r0g, r1g, r2g, radiong, &
    nion2, radjel, surjel, bbeta, gamma, beta4

  OPEN (5, STATUS='old', FORM='FORMATTED', FILE='for005.'//outnam)
  READ (5, perio, END=99999)
  WRITE (*, *) ' PERIO READ'
99999 CONTINUE


END SUBROUTINE changeperio

!-----iparams--------------------------------------------------

SUBROUTINE iparams()

! Check consistency of input parameters, do some initializations

  USE params
  USE RTA_module
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE calc_lda_gunnar(rho, chpdft)
      USE params, ONLY: DP, kdfull2
      REAL(DP), INTENT(IN) :: rho(2*kdfull2)
      REAL(DP), INTENT(OUT) :: chpdft(2*kdfull2)
    END SUBROUTINE calc_lda_gunnar
  END INTERFACE
  INTERFACE
    SUBROUTINE calc_lda_pw92(rho, chpdft)
      USE params, ONLY: DP, kdfull2
      REAL(DP), INTENT(IN) :: rho(2*kdfull2)
      REAL(DP), INTENT(OUT) :: chpdft(2*kdfull2)
    END SUBROUTINE calc_lda_pw92
  END INTERFACE


  INTEGER :: itime, irate

!--------------------------------------------------------------

! check undefined input parameters
  IF(nelect < 0) STOP "NELECT must be given a value"
  IF(nspdw < 0) STOP "NSPDW must be given a value"
  IF(nion < 0) STOP "NION must be given a value"
  IF(nion2 == 1 .AND. ipsptyp < 0) STOP "IPSPTYP must be given a value"
  IF(dt1 < 0D0) STOP "DT1 must be given a value"

! some initializations
  phi = phi*pi/180D0 ! convert input 'phi' from degree

! determine factor to PRINT system time in seconds
  CALL system_clock(itime, irate)
  systime_factor = 1D0/irate


  qold2 = 0.01D0
  qold3 = 0.01D0
  qold4 = 0.01D0
  jekion = 0
  iquery4 = 0
  WRITE (6, *) 'jekion=', jekion

! check consistency of input options

  IF (numspin .NE. 2 .AND. iftransme) STOP ' IFTRANSME needs full spin'

  IF (nion > ng) STOP " more ions than dimensioned. enhance NG"

  IF (myn == 0) THEN
    IF (irest > itmax) STOP ' IREST > ITMAX is nonsense'
  END IF

  IF (nion2 == 0 .AND. ifsicp >= 3) &
    STOP 'Jellium not compatible with KLI or Slater-SIC'
  IF (nion2 == 2 .AND. ifsicp >= 3) &
    STOP 'EXTERNAL PsP not compatible with KLI or Slater-SIC'


  IF (ifhamdiag > 0 .AND. ifsicp == 5) STOP &
   ' H diagonalization incompatible with exact exchange (IFSICP=5)'

  directenergy = .TRUE. ! .FALSE.
  IF (ifsicp == 3 .OR. ifsicp == 4) directenergy = .TRUE.

  IF (numspin .NE. 2 .AND. ifsicp >= 3) STOP 'IFSICP>2 requires fullspin code'



  IF (itmax > 0 .AND. ifexpevol .AND. ionmdtyp /= 0) &
    STOP ' exponential evolution not yet with ionic motion'

  IF (itmax > 0 .AND. nabsorb == 0 .AND. jesc .NE. 0) &
    STOP ' JESC must be zero for NABSORB=0'

  IF (jmp > 0 .AND. nabsorb == 0) STOP 'PES requires absorbing bounds'


  IF (nion2 == 2) THEN
    IF (irotat > 0) STOP 'IROTAT>0 not possible for EXTERNAL psp (NION2=2)'
    IF (ipsptyp /= 0) STOP ' IPSPTYP=0 needed for EXTERNAL psp (NION2=2)'
    IF (ionmdtyp /= 0) STOP ' fixed ions needed for EXTERNAL psp (IONMDTYP=0)'
  END IF

  IF(abs(centfx)+abs(centfy)+abs(centfz)+abs(shiftinix)+&
     abs(shiftiniy)+abs(shiftiniz) > 0D0 .AND. irotat > 0) &
    STOP "rotational initialization incompatible with boost/shift"

! set the pointer for the energy-density functional
  IF (idenfunc == 1) THEN
    calc_lda => calc_lda_pw92
  ELSE IF (idenfunc == 2 .OR. idenfunc == 3) THEN
    calc_lda => calc_lda_gunnar
  ELSE
    STOP ' invalid value for density-functional selector IDENFUNC'
  END IF

! check number of ions if force-print is wanted
  IF(jforce > 0 .AND. nion > 9) STOP "print of ionic forces restricted to less than 10 ions"

! check initial RTA values
  IF (jrtaint > 0) THEN!MV init ONLY IF used
    IF (rtaeps .LE. 0D0) rtaeps = epswf ! IF rtaeps not specified, use static value
    IF (rtae0dmp .LE. 0D0) rtae0dmp = e0dmp ! IF rtae0dmp not specified, use static value
    IF (rtamu < 0D0) STOP 'RTAMU must be specified'
    IF (rtamuj < 0D0) STOP 'RTAMUJ must be specified'
    IF (rtasigee < 0D0) STOP "RTASIGEE must be given a value"
    IF (rtars < 0D0) STOP "RTARS must be given a value"
  END IF

  RETURN
END SUBROUTINE iparams

!-----ocoption--------------------------------------------------

SUBROUTINE ocoption(iu)

! Protocol of options on 'poptions.<NAME>'.

  USE params
  USE coulsolv, ONLY: tcoulfalr
  IMPLICIT NONE
  INTEGER :: iu
!----------------------------------------------------------------

  IF (iu == 8) OPEN (8, STATUS='unknown', FORM='FORMATTED', &
                     FILE='poptions.'//outnam)

  WRITE (iu, '(a)') 'the following options are used:'
  WRITE (iu, *)
  IF (idenfunc == 1) THEN
    WRITE (iu, '(a)') 'Perdew-Wang 92 exchange-correlation functional'
  ELSE IF (idenfunc == 2) THEN
    WRITE (iu, '(a)') 'gl 76 exchange-correlation functional'
  ELSE IF (idenfunc == 3) THEN
    WRITE (iu, '(a)') 'pure exchange functional in LDA'
  ELSE
    STOP ' invalid value for density-functional selector IDENFUNC'
  END IF

! parallel or serial:
  WRITE (iu, *) 'serial code: max. NUMBER of wf=', kstate
  IF (numspin == 2) THEN
    WRITE (iu, '(a)') ' spin code '
  ELSE IF (numspin == 1) THEN
    WRITE (iu, '(a)') ' no-spin code '
  ELSE
    STOP "invalid input parameter NUMSPIN"
  END IF
  IF (ifhamdiag > 0) THEN
    WRITE (iu, '(a)') ' static step: Hamiltonian in subspace diagonalized'
  ELSE
    WRITE (iu, '(a)') ' static step: Hamiltonian in subspace not diagonalized'
  END IF
! fft:
  WRITE (iu, '(a)') 'Fourier propagation using FFTW3 library'
  IF (tcoulfalr) THEN
    WRITE (iu, '(a)') 'falr coulomb solver'
  ELSE
    WRITE (iu, '(a)') 'exact coulomb solver'
  END IF

  IF (iffastpropag) WRITE (iu, '(a)') 'faster TV propagator'

  IF (ifexpevol) THEN
    WRITE (iu, '(a)') 'time step by exponential evolution'
  ELSE
    WRITE (iu, '(a)') 'time step by T-V splitting'
  END IF
  IF (tdipolxyz) WRITE (iu, '(a)') 'hom. electric field switched on'


! output compiled parameters

  WRITE (iu, *) '_______________________________________________'
  WRITE (iu, *)
  WRITE (iu, *) 'compiled parameters'
  WRITE (iu, *)
  WRITE (iu, *) 'maximum compiled sizes x y z', kxbox, kybox, kzbox
  WRITE (iu, *) 'maximum number of ions', ng
  WRITE (iu, *)
  WRITE (iu, *) 'units'
  WRITE (iu, *)
  WRITE (iu, *) 'h bar', hbar
  WRITE (iu, *) 'e2', e2
  WRITE (iu, *) 'm electron', ame
  WRITE (iu, *) '_______________________________________________'
  WRITE (iu, *)

! pseudos
  IF (ipsptyp == 0) THEN
    WRITE (iu, '(a)') 'soft local pseudopotentials (errf)'
  ELSE IF (ipsptyp == 1) THEN
    WRITE (iu, '(a)') 'full Goedecker pseudopotentials'
  ELSE IF (ipsptyp == 2) THEN
    WRITE (iu, '(a)') 'local Goedecker pseudopotentials'
  ELSE IF (ipsptyp == 3) THEN
    WRITE (iu, '(a)') 'full Goedecker pseudopotentials from FILE goed.asci'
  ELSE IF (ipsptyp == 4) THEN
    WRITE (iu, '(a)') 'full semicore Goedecker pseudopotentials from FILE goed.asci'
  ELSE
    STOP ' this type IPSPTYP not yet implemented'
  END IF
  IF (init_ao) THEN
    WRITE (iu, '(a)') ' initialize wavefunctions with atomic orbitals'
  ELSE
    WRITE (iu, '(a)') ' initialize wavefunctions with h.o.'
  END IF
! abs bounds
  IF (nabsorb > 0) THEN
    IF (ispherabso == 1) THEN
      WRITE (iu, '(a,i4,1pg13.5)') &
        ' spherical absorbing b.c.: nabsorb,powabso=', nabsorb, powabso
    ELSE
      WRITE (iu, '(a,i4,1pg13.5)') &
        ' cartesian absorbing b.c.: nabsorb,powabso=', nabsorb, powabso
    END IF
  ELSE
    WRITE (iu, '(a)') ' no absorbing b.c.'
  END IF
! MD TYPE
  IF (ionmdtyp == 0) THEN
    WRITE (iu, '(a)') 'no molecular dynamics '
  ELSE IF (ionmdtyp == 1) THEN
    WRITE (iu, '(a)') 'molecular dynamics: leap-frog step'
  ELSE IF (ionmdtyp == 2) THEN
    WRITE (iu, '(a)') 'molecular dynamics: velocity Verlet'
  ELSE
    STOP ' this option for MD not valid '
  END IF
  IF (tfixcmion) WRITE (iu, *) ' fix ionic c.m. during ionic motion'
  IF (ifredmas) WRITE (iu, '(a)') ' reduced ionic mass = half proton mass '
  IF (icooltyp == 1) THEN
    WRITE (iu, '(a)') 'cooling with pseudo dynamics '
  ELSE IF (icooltyp == 2) THEN
    WRITE (iu, '(a)') 'cooling with steepest descent '
  ELSE IF (icooltyp == 3) THEN
    WRITE (iu, '(a)') 'cooling with Monte Carlo '
  END IF

  IF (directenergy) THEN
    WRITE (iu, '(a)') ' energy computed directly'
  ELSE
    WRITE (iu, '(a)') ' energy computed using s.p. energies'
  END IF

  IF (ifsicp == 0) THEN
    WRITE (iu, '(a)') 'no sic'
  ELSE IF (ifsicp == 2) THEN
    WRITE (iu, '(a)') 'sic activated: ADSIC'
  ELSE IF (ifsicp == 3) THEN
    WRITE (iu, '(a)') 'sic activated: Slater'
  ELSE IF (ifsicp == 4) THEN
    WRITE (iu, '(a)') 'sic activated: KLI'
  ELSE IF (ifsicp == 5) THEN
    WRITE (iu, '(a)') 'sic activated: exact exchange'
  ELSE IF (ifsicp == 8) THEN
    WRITE (iu, '(a)') 'sic activated: DOUBLE-set SIC - REAL'
  ELSE IF (ifsicp == 9) THEN
    WRITE (iu, '(a)') 'sic activated: DOUBLE-set SIC - COMPLEX'
  ELSE
    WRITE (iu, '(a)') 'this version of SIC not available'
    STOP "invalid version of SIC"
  END IF

! dynamical options

  WRITE (iu, '(a,3i6)') ' kxbox,kybox,kzbox=', kxbox, kybox, kzbox
  WRITE (iu, '(a,3i6)') ' nelect,nion,nstate=', nelect, nion, kstate
  WRITE (iu, '(a,l2,i3,f7.2)') ' tspindip,irotat,phirot=', &
    tspindip, irotat, phirot
  WRITE (iu, '(a,3f8.2)') ' boost: centfx,centfy,centfz=', centfx, centfy, centfz
  WRITE (iu, '(a,3f8.2)') ' shift: shiftinix,y,z=', shiftinix, shiftiniy, shiftiniz
  WRITE (iu, '(a,i3,l2,f7.3)') ' irest,tstat,dt1=', irest, tstat, dt1
  WRITE (iu, '(a/a,i3,5f8.3/a,7f8.3)') ' laser:', &
    ' itft,tnode,deltat,tpeak,omega,e0=', itft, tnode, deltat, tpeak, omega, e0, &
    ' e1x,e1y,e1z,e2x,e2y,e2z,phi=', e1x, e1y, e1z, e2x, e2y, e2z, phi


  IF (iu == 8) CLOSE (8)

  RETURN
END SUBROUTINE ocoption

!-----init_output--------------------------------------------------

SUBROUTINE init_output()

! Initializes output files and writes headers

  USE params
  IMPLICIT NONE
  INTEGER :: i, j, maxnum
  REAL(DP) :: ascal
  CHARACTER(LEN=3) :: num

!------------------------------------------------------------------

  IF (myn < 10) THEN
    WRITE (num, '(i1)') myn
    maxnum = 1
  ELSE IF (myn < 100 .AND. myn > 9) THEN
    WRITE (num, '(i2)') myn
    maxnum = 2
  ELSE
    WRITE (num, '(i3)') myn
    maxnum = 3
  END IF
  outname = trim(num)//trim(outnam)

! output grid properties

  IF (myn == 0) THEN
    WRITE (6, *) '::::::::::::::::::::::::::::::::::::::::'
    WRITE (6, *) 'Dimensions of the box:'
    WRITE (6, '(a,3i5,i8)') &
      ' KXBOX,KYBOX,KZBOX,KDFULL2=', kxbox, kybox, kzbox, kdfull2
    WRITE (6, '(a,f10.2,a,f10.2)') 'x = ', (1 - nxsh)*dx, ' .. ', nxsh*dx
    WRITE (6, '(a,f10.2,a,f10.2)') 'y = ', (1 - nysh)*dy, ' .. ', nysh*dy
    WRITE (6, '(a,f10.2,a,f10.2)') 'z = ', (1 - nzsh)*dz, ' .. ', nzsh*dz
    WRITE (6, *) '::::::::::::::::::::::::::::::::::::::::'
  END IF

! deduced grid parameters

  WRITE (6, '(a,3i5,i8,g12.4)') ' INIT: nx,ny,nz,nxyz,dvol=', &
    nx, ny, nz, nxyz, dvol

! output static parameters

  WRITE (7, '(5x,a,i2,a,i2/a/)') &
    '# electr.: ', nelect, '# ions: ', nion, '=========='
  WRITE (7, '(a,i3)') 'NUMBER of wave-functions nstate = ', nstate

  WRITE (7, '(a,f4.2)') 'initialisation : osfac=', osfac
  WRITE (7, '(a,10(/t2,3(3i4,5x),3i4))') &
    'quantumnumbers : ', ((nq(i, j), i=1, 3), j=1, nstate)
  IF (temp > 0D0) THEN
    WRITE (7, '(a,f10.5,a,f4.2)') 'temperature kt= ', temp, 'Ry, occmix=', occmix
    WRITE (7, '(a,1g13.5)') 'epsoro=', epsoro
    WRITE (6, '(a,f10.5,a,f4.2)') 'temperature kt= ', temp, 'Ry, occmix=', occmix
    WRITE (6, '(a,1g13.5)') 'epsoro=', epsoro
  END IF
  WRITE (7, '(3(a,i4,3x))') 'nx=', nx, 'ny=', ny, 'nz=', nz
  WRITE (7, '(3(a,f4.2,3x))') 'dx=', dx, 'dy=', dy, 'dz=', dz
  WRITE (7, '(a)') 'damped gradient step :'
  WRITE (7, '(2(a,1pg12.4,3x))') 'epswf=', epswf, 'e0dmp=', e0dmp
  WRITE (7, '(a)') 'jellium background :'
  WRITE (7, '(2(a,1pg12.4,3x))') 'radjel=', radjel, 'surjel=', surjel

  WRITE (7, '(a,a)') 'set of fixed deformation parameters ', &
    '(IN hill-wheeler-coordinates) :'
  WRITE (7, '(a,f4.2,3x,a,f5.0)') 'bbeta=', bbeta, 'gamma=', gamma
  WRITE (7, '(//)')


  WRITE (7, *) '_______________________________________________'
  WRITE (7, *)
  WRITE (7, *) ' dynamic parameters have been read :'
  WRITE (7, *) '_______________________________________________'
  WRITE (7, *)

  IF (nion2 == 0) THEN
    WRITE (7, *)
    WRITE (7, *) 'nions2=0 : jellium code'
    ecorr = 0D0
  ELSE IF (nion2 == 0) THEN
    WRITE (7, *) 'external local PsP from FILE'
  ELSE
    WRITE (7, *) 'ionic code - number of ions', nion, nion2
  END IF
  IF (irest == 0) THEN
    WRITE (7, *) 'number of static iterations', ismax
  ELSE
    WRITE (7, *) 'restart option activated - restart at', irest
  END IF
  WRITE (7, *) 'saving wfs every ', isaves, 'static iterations'
  WRITE (7, *) 'saving wfs every ', isaved, 'dynamic steps'
  WRITE (7, *) 'number of dynamic iterations', itmax
  WRITE (7, *) 'timestep=', dt1, 'units'

  IF (irotat == 0) THEN
    WRITE (7, *) 'dipoleshift or spindipoleshift'
    IF (tspindip) THEN
      WRITE (7, *) 'shift of spinup versus spindown density by'
      WRITE (7, *) 0.5D0*centfx, 0.5D0*centfy, 0.5D0*centfz
    ELSE
      WRITE (7, *) 'shift of electronic density by', centfx, centfy, centfz
    END IF
  ELSE
    IF (irotat > 0) THEN
      WRITE (7, *) 'scissor mode excitation'
      WRITE (7, *) 'axis of rotation', irotat
      WRITE (7, *) 'angle of rotation', phirot
    END IF
  END IF
  IF (tempion /= 0D0) THEN
    WRITE (7, *) 'initial thermalization of the ions at t=', tempion
  ELSE
    WRITE (7, *) 'zero initial ionic kinetic energy'
  END IF

! laser

  IF (e0 /= 0D0) THEN
    WRITE (7, *) 'laser pulse active - parameters : '
    IF (itft == 1) THEN
      WRITE (7, *) 'ramp laser pulse, sine switching'
      WRITE (7, *) ' on/off up to/from f(t) = 1'
      IF (tpeak == 0D0) &
        STOP 'tpeak should be non-zero in the case of the ramp laser pulse'
    ELSE IF (itft == 2) THEN
      WRITE (7, *) 'Gaussian laser pulse'
      IF (deltat == 0D0) &
        STOP 'deltat should be non-zero in the case of a Gaussian laser pulse'
    ELSE IF (itft == 3) THEN
      WRITE (7, *) 'cos**2 pulse'
      IF (deltat == 0D0) &
        STOP 'deltat should be non-zero in the case of a cos**2 laser pulse'
    ELSE IF (itft == 2) THEN
      WRITE (7, *) 'cos**4 pulse'
      IF (deltat == 0D0) &
        STOP 'deltat should be non-zero in the case of a cos**4 pulse'
    END IF
    WRITE (7, *) 'parameters of first photon pulse:'
    WRITE (7, *) 'length of the pulse [fs]', deltat
    WRITE (7, *) 'peak time [fs]', tpeak
    WRITE (7, *) 'start time [fs]', tnode
    WRITE (7, *) 'frequency [Ry]', omega
    WRITE (7, *) 'field strength [Ry/a_0]', e0
    WRITE (7, *) 'orientation',e1x, e1y, e1z
    WRITE (7, *) 'phase',phi
    IF(e0_2 .NE. 0D0) THEN
      WRITE (7, *) 'parameters of second photon pulse:'
      WRITE (7, *) 'length of the pulse [fs]', deltat2
      WRITE (7, *) 'start time [fs]', tstart2
      WRITE (7, *) 'frequency [Ry]', omega2
      WRITE (7, *) 'field strength [Ry/a_0]', e0_2
      WRITE (7, *) 'orientation',e2x, e2y, e2z
      WRITE (7, *) 'phase',phase2
    END IF
  ELSE
    WRITE (7, *) 'photon pulses switched off'
  END IF
  WRITE (7, *) '_______________________________________________'
  WRITE (7, *)

  RETURN
END SUBROUTINE init_output

!-----iperio-----------------------------------------------------

SUBROUTINE iperio

! Initializes the default PsP parameters.
! Parameters are communicated via module 'params'.

  USE params
  IMPLICIT NONE

  INTEGER :: i, iel, natom, icountt, inew, iskip
  INTEGER :: l, nactual, nval
  REAL(DP) :: amfac, c1, c2, c3, c4, crr, h11, h22, h33, rloc
  REAL(DP) :: scpx, scpy, scpz
  CHARACTER(LEN=3) :: naml
  CHARACTER(LEN=2) :: namc, symb(99)
  CHARACTER(LEN=80) :: a

!----------------------------------------------------------------

  WRITE (6, *) 'Entering iperio. IPSPTYP=',ipsptyp

  DO iel = kpspm, kpsp
    amu(iel) = 0D0
  END DO

! soft PsP (error functions)

  IF (ipsptyp == 0) THEN

! hydrogen

    amu(1) = 1D0
    ch(1) = 1D0
    dr1(1) = 0.3D0
    dr2(1) = 0.45D0
    prho1(1) = 6.94527D0
    prho2(1) = -0.92056D0
    sgm1(1) = dr1(1)*0.8493218D0
    sgm2(1) = dr2(1)*0.8493218D0
    chg1(1) = ((sgm1(1)*SQRT(pi*2D0))**3)*prho1(1)*ch(1)
    chg2(1) = ((sgm2(1)*SQRT(pi*2D0))**3)*prho2(1)*ch(1)
    WRITE (6, *) 'ch(1)=', ch(1)
    WRITE (6, *) 'total charge of Psphydr=', chg1(1) + chg2(1)

! sodium

    amu(11) = 23.0D0
    ch(11) = 1D0 ! charge of pseudopotential
    dr1(11) = 0.8018D0
    dr2(11) = 1.3693D0
    prho1(11) = -0.46073D0
    prho2(11) = 0.13287D0
    sgm1(11) = dr1(11)*0.8493218D0
    sgm2(11) = dr2(11)*0.8493218D0
    chg1(11) = ((sgm1(11)*SQRT(pi*2D0))**3)*prho1(11)*ch(11)
    chg2(11) = ((sgm2(11)*SQRT(pi*2D0))**3)*prho2(11)*ch(11)
    WRITE (6, *) 'ch(11)=', ch(11)
    WRITE (6, *) 'total charge of Pspsodi=', chg1(11) + chg2(11)


! potassium

    amu(19) = 39.0D0
    ch(19) = 1D0 ! charge of pseudopotential
    dr1(19) = 0.9973D0
    dr2(19) = 1.9957D0
    prho1(19) = -0.13771D0
    prho2(19) = 0.030223D0
    sgm1(19) = dr1(19)*0.8493218D0
    sgm2(19) = dr2(19)*0.8493218D0
    chg1(19) = ((sgm1(19)*SQRT(pi*2D0))**3)*prho1(19)*ch(19)
    chg2(19) = ((sgm2(19)*SQRT(pi*2D0))**3)*prho2(19)*ch(19)
    WRITE (6, *) 'ch(19)=', ch(19)
    WRITE (6, *) 'total charge of Psppotassium=', chg1(19) + chg2(19)

! cesium

    amu(58) = 140D0
    ch(58) = 1D0 ! charge of pseudopotential
    dr1(58) = 1.12D0
    dr2(58) = 2.552D0
    prho1(58) = -0.62923D-01
    prho2(58) = 0.11554D-01
    sgm1(58) = dr1(58)*0.8493218D0
    sgm2(58) = dr2(58)*0.8493218D0
    chg1(58) = ((sgm1(58)*SQRT(pi*2D0))**3)*prho1(58)*ch(58)
    chg2(58) = ((sgm2(58)*SQRT(pi*2D0))**3)*prho2(58)*ch(58)
    WRITE (6, *) 'ch(58)=', ch(58)
    WRITE (6, *) 'total charge of Psppotassium=', chg1(58) + chg2(58)

! standard Goedecker PsP (local part initialized here)

  ELSE IF (ipsptyp == 1) THEN
! non-local PsP from S. Goedecker, M. Teter, and J. Hutter Phys. Rev. B 54, 1703 (1996)

! hydrogen

    amu(1) = 1D0
    ch(1) = 1D0
    cc1(1) = -4.180237D0
    cc2(1) = 0.725075D0
    crloc(1) = 0.2D0

! helium

    amu(2) = 4.0D0
    ch(2) = 2.0D0
    cc1(2) = -9.112023D0
    cc2(2) = 1.698368D0
    crloc(2) = 0.2D0

! bor

    amu(5) = 10.81D0
    ch(5) = 3.0D0
    cc1(5) = -5.578642D0
    cc2(5) = 0.804251D0
    crloc(5) = 0.43393D0
    r0g(5) = 0.373843D0
    r1g(5) = 0.360393D0
    radiong(5) = 2.5D0
    h0_11g(5) = 6.233928D0
    h1_11g(5) = 0D0

! carbon

    amu(6) = 12.0D0
    ch(6) = 4.0D0
! Goedecker I Goedecker II
    cc1(6) = -8.5753285 !-8.513771D0
    cc2(6) = 1.2341279 !1.228432D0
    crloc(6) = 0.3464730 !0.348830D0
    r0g(6) = 0.3045228 !0.304533D0
    r1g(6) = 0.232677D0
    radiong(6) = 2D0
    h0_11g(6) = 9.5341929 !9.522842D0
    h1_11g(6) = 0D0

! nitrogen

    amu(7) = 14.0D0
    ch(7) = 5.0D0
    cc1(7) = -12.2046419D0
    cc2(7) = 1.7558249D0
    crloc(7) = 0.2889046D0
    r0g(7) = 0.256605D0
    r1g(7) = 0.270134D0
    radiong(7) = 2.0D0
    h0_11g(7) = 13.552433D0
    h1_11g(7) = 0D0

! oxygen

    amu(8) = 16.0D0
    ch(8) = 6.0D0
    cc1(8) = -16.4822284D0 ! -16.580318
    cc2(8) = 2.3701353D0 ! 2.395701
    crloc(8) = 0.2477535D0 ! 0.247621
    r0g(8) = 0.2222028D0 ! 0.221786D0
    r1g(8) = 0.256829D0
    radiong(8) = 1.5D0
    h0_11g(8) = 18.19996387D0 ! 18.266917D0
    h1_11g(8) = 0D0

! fluor (row 2)

    amu(9) = 18D0
    ch(9) = 7.0D0
    cc1(9) = -21.307361D0
    cc2(9) = 3.072869D0
    crloc(9) = 0.218525D0
    r0g(9) = 0.195567D0
    r1g(9) = 0.2D0
    radiong(9) = 1.3D0
    h0_11g(9) = 23.58494D0
    h1_11g(9) = 0D0

! neon (row 2)

    amu(10) = 20.2D0
    ch(10) = 8.0D0
    cc1(10) = -27.692852D0
    cc2(10) = 4.005906D0
    crloc(10) = 0.19D0
    r0g(10) = 0.174268D0
    r1g(10) = 0.214913D0
    h0_11g(10) = 28.506098D0
    h0_22g(10) = -1.076245D0
    h1_11g(10) = -0.000090D0

! sodium

    amu(11) = 23.0D0
    ch(11) = 1D0
    cc1(11) = -1.238867D0
    cc2(11) = 0D0
    crloc(11) = 0.885509D0
    r0g(11) = 0.661104D0
    r1g(11) = 0.857119D0
    h0_11g(11) = 1.847271D0
    h0_22g(11) = 0.582004D0
    h1_11g(11) = 0.471133D0

! magnesium

    amu(12) = 24.312D0
    ch(12) = 2.0D0
    cc1(12) = -2.864297D0
    cc2(12) = 0D0
    crloc(12) = 0.651812D0
    r0g(12) = 0.556478D0
    r1g(12) = 0.677569D0
    h0_11g(12) = 2.970957D0
    h0_22g(12) = 1.329941D0
    h1_11g(12) = 1.049881D0

! aluminium

    amu(13) = 26.98D0
    ch(13) = 3.0D0
    cc1(13) = -8.491315D0
    cc2(13) = 0D0
    crloc(13) = 0.45D0
    r0g(13) = 0.460104D0
    r1g(13) = 0.536744D0
    h0_11g(13) = 5.088340D0
    h0_22g(13) = 2.679700D0
    h1_11g(13) = 2.193438D0

! silicon

    amu(14) = 28.09D0
    ch(14) = 4.0D0
    cc1(14) = -6.9136286D0 !-7.336103D0
    cc2(14) = 0D0
    crloc(14) = 0.44D0
    r0g(14) = 0.4243338D0
    r1g(14) = 0.4853587D0
    h0_11g(14) = 3.2081318D0
    h0_12g(14) = 0D0
    h0_22g(14) = 2.5888808D0
    h1_11g(14) = 2.6562230D0

! phosphor

    amu(15) = 30.974D0
    ch(15) = 5.0D0
    cc1(15) = -6.65422D0
    cc2(15) = 0D0
    crloc(15) = 0.43D0
    r0g(15) = 0.389803D0
    r1g(15) = 0.440796D0
    h0_11g(15) = 6.842136D0
    h0_22g(15) = 3.856693D0
    h1_11g(15) = 3.282606D0

    ! sulfur
    amu(16) = 32.06D0
    ch(16) = 6.0D0 ! Goed1 ! Goed2
    cc1(16) = -6.5960716D0 !-6.8903645 !-6.554492 !
    cc2(16) = 0D0
    crloc(16) = 0.42 !0.41D0 !0.42D0 !
    r0g(16) = 0.3626413D0 !0.3389943 !0.361757D0 !
    r1g(16) = 0.405322D0 !0.3762100 !0.405285D0 !
    h0_11g(16) = 4.2228399D0 !4.9069762 !7.905303 !
    h0_12g(16) = 0D0
    h0_22g(16) = 3.66696625D0 !4.1601818 !4.471698 !
    h1_11g(16) = 3.8853458D0 !4.4850412 !3.866579 !
    radiong(16) = 3.0

! argon (row 3)

    amu(18) = 40D0
    ch(18) = 8.0D0
    cc1(18) = -7.1D0
    cc2(18) = 0D0
    crloc(18) = 0.4D0
    r0g(18) = 0.317381D0
    r1g(18) = 0.351619D0
    h0_11g(18) = 10.249487D0
    h0_22g(18) = 5.602516D0
    h1_11g(18) = 4.978801D0

! calcium

    amu(20) = 40.08D0
    ch(20) = 2.0D0
    cc1(20) = 0D0
    cc2(20) = 0D0
    crloc(20) = 0.8D0
    r0g(20) = 0.669737D0
    r1g(20) = 0.946474D0
    r2g(20) = 0.526550D0
    h0_11g(20) = 1.645014D0
    h0_22g(20) = 1.523491D0
    h0_33g(20) = 0.295996D0
    h1_11g(20) = 0.585479D0
    h1_22g(20) = 0.126329D0
    h2_11g(20) = -3.0323321D0

! Cu

    amu(29) = 63.546D0
    ch(29) = 1D0
    cc1(29) = 0D0
    cc2(29) = 0D0
    crloc(29) = 0.58D0
    r0g(29) = 0.843283D0
    r1g(29) = 1.089543D0
    r2g(29) = 1.291602D0
    h0_11g(29) = 0.975787D0
    h0_22g(29) = -0.822070D0
    h0_33g(29) = -0.133237D0
    h1_11g(29) = 0.024580D0
    h1_22g(29) = -0.249001D0
    h2_11g(29) = -0.065292D0

! Ag

    amu(47) = 107.8682D0
    ch(47) = 1D0
    cc1(47) = -2.376061D0
    cc2(47) = 0D0
    crloc(47) = 0.65D0
    r0g(47) = 1.012705D0
    r1g(47) = 1.235842D0
    r2g(47) = 1.016159D0
    h0_11g(47) = 0.897931D0
    h0_22g(47) = -0.7483230D0
    h0_33g(47) = 0.029787D0
    h1_11g(47) = 0.130081D0
    h1_22g(47) = -0.277495D0
    h2_11g(47) = -0.038842D0

! set radius of PsP cell to embrace all contributions of an element
    DO i = 1, kpsp
      radiong(i) = 4.5D0*MAX(crloc(i), r0g(i), r1g(i), r2g(i))
    END DO

! purely local Goedecker PsP

  ELSE IF (ipsptyp == 2) THEN

! hydrogen Gianocci PsP (approx. local Goedecker)

    amu(1) = 1D0
    ch(1) = 1D0

! Na (approximate local pseudo built on Goedecker local part)

    amu(11) = 23.0D0
    ch(11) = 1D0
    cc1(11) = 1.5D0
    cc2(11) = 0D0
    crloc(11) = 0.9D0

! Ar (approximate local pseudo built on Goedecker local part)

    amu(18) = 39.95D0
    ch(18) = 8D0
    cc1(18) = 2.482D0
    cc2(18) = -1.4526D0
    crloc(18) = 0.9D0


    STOP ' IPERIO: this TYPE PsP not yet implemented '
  END IF

! reset ionic masses if asked for by IFREDMAS

  WRITE (6, *) 'resetting ionic masses...'
  IF (ifredmas) THEN
    IF (ipsptyp == 0) THEN
      DO i = kpspm, kpsp
        amu(i) = 0.5D0
      END DO
    ELSE
      amfac = 1D0/amu(6)
      DO i = kpspm, kpsp
        amu(i) = amu(i)*amfac
        amu(i) = 0.5D0
      END DO
    END IF
    WRITE (6, *) 'reduced mass is a half of hydrogen'
    scpx = 0D0
    scpy = 0D0
    scpz = 0D0
  END IF

  WRITE (6, *) 'END OF IPERIO.'

  RETURN
END SUBROUTINE iperio

!-----initions-------------------------------------------------

SUBROUTINE initions()

! read ionic positions and initialize ion-related parameters

  USE params
  USE util, ONLY: getcm, givetemperature, rotatevec3D
  IMPLICIT NONE
  CHARACTER(LEN=3) :: orderxyz
  REAL(DP) :: dd, distmax, xcm, ycm, zcm, ctmpx, ctmpy, ctmpz
  REAL(DP) :: optis, dgrid, sumion, v0, rnorm, tempv, xm
  REAL(DP) :: vxn02, vyn02, vzn02
  REAL(DP) :: vecin(3), vecout(3), vecalpha(3), totvalec = 0D0
  INTEGER :: i, ii, inx, inxg, ion, iunit, nxopti

  INTEGER :: igrid(7)
  DATA igrid/32, 48, 64, 72, 96, 128, 160/

  REAL(DP), EXTERNAL :: energ_ions

  IF (nion2 == 0) THEN
    ecorr = 0D0
    RETURN ! this is the case of jellium
  END IF

  WRITE (6, *) 'Entering initions()'


  OPEN (UNIT=9, STATUS='old', FORM='FORMATTED', FILE='for005ion.'//outnam)

  WRITE (6, *) 'Reading positions...'
  distmax = 0D0
  DO ion = 1, nion

    IF (init_ao) THEN
      READ (9, *) cx(ion), cy(ion), cz(ion), np(ion), orderxyz, radini(ion) &
        , ipol(ion)
      dd = sqrt(cx(ion)*cx(ion) + cy(ion)*cy(ion) + cz(ion)*cz(ion))
      IF (dd .gt. distmax) distmax = dd

      totvalec = totvalec + ch(np(ion))
! translate ordering of atomic states
      IF (orderxyz == 'xyz') THEN
        initord(1, ion) = 1
        initord(2, ion) = 2
        initord(3, ion) = 3
      ELSE IF (orderxyz == 'xzy') THEN
        initord(1, ion) = 1
        initord(2, ion) = 3
        initord(3, ion) = 2
      ELSE IF (orderxyz == 'yxz') THEN
        initord(1, ion) = 2
        initord(2, ion) = 1
        initord(3, ion) = 3
      ELSE IF (orderxyz == 'yzx') THEN
        initord(1, ion) = 2
        initord(2, ion) = 3
        initord(3, ion) = 1
      ELSE IF (orderxyz == 'zxy') THEN
        initord(1, ion) = 3
        initord(2, ion) = 1
        initord(3, ion) = 2
      ELSE IF (orderxyz == 'zyx') THEN
        initord(1, ion) = 3
        initord(2, ion) = 2
        initord(3, ion) = 1
      ELSE
        STOP ' this ordering of XYZ not provided '
      END IF
      WRITE (6, *) ' ion,initord=', ion, initord(:, ion)
    ELSE
      READ (9, *) cx(ion), cy(ion), cz(ion), np(ion)
      WRITE (*, *) cx(ion), cy(ion), cz(ion), np(ion), ipsptyp
      totvalec = totvalec + ch(np(ion))
    END IF

! initial kinetic momenta=0

    cpx(ion) = 0D0
    cpy(ion) = 0D0
    cpz(ion) = 0D0

  END DO
  CLOSE (UNIT=9)

  WRITE (6, *) 'total number of valence electrons', totvalec

! re-scale molecule if desired

  IF (ABS(scaleion - 1D0) > small .OR. ABS(scaleionx - 1D0) > small .OR. &
      ABS(scaleiony - 1D0) > small .OR. ABS(scaleionz - 1D0) > small) THEN

    WRITE (6, *) 'Rescaling...'

    IF (scaleion /= 1D0) THEN
      scaleionx = scaleion
      scaleiony = scaleion
      scaleionz = scaleion
    END IF

    xcm = 0D0
    ycm = 0D0
    zcm = 0D0
    DO i = 1, nion
      xcm = xcm + cx(i)
      ycm = ycm + cy(i)
      zcm = zcm + cz(i)
    END DO
    xcm = xcm/nion
    ycm = ycm/nion
    zcm = zcm/nion
    DO i = 1, nion
      cx(i) = (cx(i) - xcm)*scaleionx + xcm
      cy(i) = (cy(i) - ycm)*scaleiony + ycm
      cz(i) = (cz(i) - zcm)*scaleionz + zcm
    END DO

  END IF

! re-scaling done.

! shift cluster in space if desired

  IF (shiftionx /= 0D0 .OR. shiftiony /= 0D0 .OR. &
      shiftionz /= 0D0) THEN

    WRITE (6, *) 'Shifting cluster ions...'

    DO i = 1, nion
      cx(i) = cx(i) + shiftionx
      cy(i) = cy(i) + shiftiony
      cz(i) = cz(i) + shiftionz
    END DO

  END IF

  IF (tshiftcmtoorigin) THEN
    CALL getcm(1, 0, 0) ! c.m. now on 'rvectmp(1:3)'
    rvectmp2(1) = rvectmp(1)
    rvectmp2(2) = rvectmp(2)
    rvectmp2(3) = rvectmp(3)
    DO ii = 1, nion
      cx(ii) = cx(ii) - rvectmp(1)
      cy(ii) = cy(ii) - rvectmp(2)
      cz(ii) = cz(ii) - rvectmp(3)
    END DO
  END IF

! shift done

! Rotate cluster along rotation vector
! (/rotionx,rotiony,rotionz/).
! Rotation vector is given initially in degree and converted
! internally to radian.

  IF (abs(rotionx) + abs(rotiony) + abs(rotionz) > 0D0) THEN

    vecalpha(1) = rotionx/180D0*pi
    vecalpha(2) = rotiony/180D0*pi
    vecalpha(3) = rotionz/180D0*pi
    !WRITE(*,*) ' rotate ions by anglex,y,z=',vecalpha

! temporary shift to CM
    ctmpx = SUM(cx(1:nion))/nion
    ctmpy = SUM(cy(1:nion))/nion
    ctmpz = SUM(cz(1:nion))/nion
    cx = cx - ctmpx
    cy = cy - ctmpy
    cz = cz - ctmpz

! apply rotation
    DO ion = 1, nion

      vecin(1) = cx(ion)
      vecin(2) = cy(ion)
      vecin(3) = cz(ion)

      CALL rotatevec3D(vecin, vecout, vecalpha)

      cx(ion) = vecout(1)
      cy(ion) = vecout(2)
      cz(ion) = vecout(3)

    END DO

! shift center of mass back to original POSITION
    cx = cx + ctmpx
    cy = cy + ctmpy
    cz = cz + ctmpz

  END IF
! rotation completed

  tnonlocany = .false.
  optis = 1D10
  DO ion = 1, nion
    tblock(ion) = .FALSE.
    IF (ipsptyp == 0) THEN
      tblock(ion) = .NOT. (ABS(np(ion)) > 99)
      np(ion) = MOD(np(ion), 100)
      WRITE (6, '(a,i4,a,3g12.4,i4,1x,l1)') &
        ' ion nr.', ion, ': x,y,z,type,block=', &
        cx(ion), cy(ion), cz(ion), np(ion), tblock(ion)

      IF (np(ion) /= 11 .AND. np(ion) /= 12 .AND. np(ion) /= 1 &
          .AND. np(ion) /= 18 .AND. np(ion) /= -18 &
          .AND. np(ion) /= 19 .AND. np(ion) /= 58) THEN
        WRITE (6, '(a,1i6,a,1i6)') 'np(', ion, ')=', np(ion)
        STOP 'element not provided with erf pseudo potential'
      END IF
    END IF
    WRITE (6, '(a,i4,a,3(g12.4),i4,1x,2(1pg13.5))') &
      ' ion nr.', ion, ': x,y,z,TYPE,params=', &
      cx(ion), cy(ion), cz(ion), np(ion), amu(np(ion)), crloc(np(ion))
    IF (np(ion) > 92) STOP 'element out of range'
    IF (amu(np(ion)) == 0D0) STOP 'unknown elem. found'

    ! set flag for non-local PsP

    IF (ipsptyp == 1) THEN
      tnonloc(ion) = (ABS(h2_11g(np(ion))) + ABS(h1_22g(np(ion))) &
                      + ABS(h0_33g(np(ion))) + ABS(h1_11g(np(ion))) &
                      + ABS(h0_22g(np(ion))) + ABS(h0_11g(np(ion)))) &
                     .GT. small
      IF (tnonloc(ion)) tnonlocany = .true.
    END IF

    WRITE (7, *) 'np(', ion, ')=', np(ion)
    WRITE (6, *) 'np(', ion, ')=', np(ion)
    WRITE (7, *) 'amu(', np(ion), ')=', amu(np(ion))
    WRITE (6, *) 'amu(', np(ion), ')=', amu(np(ion))
    DO iunit = 6, 7
      WRITE (iunit, *) ' PsP parameters for ion:'
      WRITE (iunit, '(a,1pg14.6)') 'ch=', ch(np(ion)), &
        'amu=', amu(np(ion)), 'cc1=', cc1(np(ion)), &
        'cc2=', cc2(np(ion)), 'crloc=', crloc(np(ion)), &
        'crs=', crs(np(ion)), 'chs=', chs(np(ion)), &
        'h0_11g=', h0_11g(np(ion)), 'h0_12g=', h0_12g(np(ion)), &
        'h0_22g=', h0_22g(np(ion)), &
        'h0_33g=', h0_33g(np(ion)), 'h1_11g=', h1_11g(np(ion)), &
        'h1_22g=', h1_22g(np(ion)), 'h2_11g=', h2_11g(np(ion)), &
        'sgm1=', sgm1(np(ion)), 'sgm2=', sgm2(np(ion)), &
        'chg1=', chg1(np(ion)), 'chg2=', chg2(np(ion)), &
        'r0g=', r0g(np(ion)), &
        'r1g=', r1g(np(ion)), 'r2g=', r2g(np(ion)), &
        'radiong=', radiong(np(ion))
      WRITE (iunit, *) 'tblock=', tblock(ion)
      IF (ipsptyp == 1) WRITE (iunit, *) 'tnonloc=', tnonloc(ion)
    END DO

    IF (ipsptyp == 1) THEN
      dgrid = crloc(np(ion))/0.8493218D0
      IF (h0_11g(np(ion)) .NE. 0D0) dgrid = MIN(dgrid, r0g(np(ion))/0.8493218D0)
      IF (h1_11g(np(ion)) .NE. 0D0) dgrid = MIN(dgrid, r1g(np(ion))/0.8493218D0)
      IF (h2_11g(np(ion)) .NE. 0D0) dgrid = MIN(dgrid, r2g(np(ion))/0.8493218D0)
      WRITE (7, *) ' optimal grid spacing for this element=', dgrid
      WRITE (6, *) ' optimal grid spacing for this element=', dgrid
      IF (dgrid .lt. optis) optis = dgrid
    END IF
  END DO
  WRITE (7, *) 'dx=', dx, ' optimal grid spacing =', optis
  WRITE (6, *) 'dx=', dx, ' optimal grid spacing =', optis

  nxopti = (distmax + 1D0)/optis
  nxopti = (nxopti/4)*8
  inx = 0
  inxg = 0
  DO WHILE (inx <= nxopti)
    inxg = inxg + 1
    IF (inxg .lt. 7) THEN
      inx = igrid(inxg)
    ELSE
      inx = nxopti
    END IF
  END DO
  WRITE (6, *) 'distmax of ions', distmax, 'recommended nx ny nz', inx

! optionally set ionic velocities to simulate given temperature
  IF (tempion > 0D0 .AND. ionmdtyp /= 0) THEN
    CALL givetemperature(cpx, cpy, cpz, nion, tempion, amu(np(1))*1836.0D0*ame, 4)
  END IF


! Part for node "0" finished. Now distribute to all nodes.


! np(0) for use in Monte-Carlo

  IF (icooltyp == 3) CALL cenmass()
  np(0) = np(1)
  dt12 = dt1

! initialize pseudopotential background

  IF (ipsptyp == 1) THEN
  DO ion = 1, nion
    CALL calc_proj(cx(ion), cy(ion), cz(ion), cx(ion), cy(ion), cz(ion), ion)
  END DO
  END IF

! background energy

  WRITE (6, *) 'Calculating ionic energy ...'

  sumion = energ_ions()
  WRITE (7, '(a,f17.4)') 'sumion=', sumion
  WRITE (6, '(a,f17.4)') 'sumion=', sumion
  WRITE (7, *)
  WRITE (6, *)

  ecorr = sumion



! short protocol

  IF (myn == 0) THEN
    WRITE (7, '(a)') 'initial ionic positions and velocities:'
    WRITE (6, '(a)') 'initial ionic positions and velocities:'
    DO ion = 1, nion
      WRITE (7, '(a,i2,a,6f9.4)') 'ion', ion, '=', cx(ion), cy(ion), cz(ion), &
        cpx(ion), cpy(ion), cpz(ion)
      WRITE (6, '(a,i2,a,6f9.4)') 'ion', ion, '=', cx(ion), cy(ion), cz(ion), &
        cpx(ion), cpy(ion), cpz(ion)
    END DO
  END IF

  WRITE (6, *) 'Initions Done.'

  RETURN
END SUBROUTINE initions

!-----init_jellium-------------------------------------------------

SUBROUTINE init_jellium()

! Initialize jellium background

  USE params
  IMPLICIT NONE

  INTEGER :: ind
  REAL(DP):: a20fac, a22fac, gamarg, xclust
  REAL(DP):: srms, sum1
!------------------------------------------------------------------

  IF (nion2 /= 0) RETURN ! case of detailed ions

! factors for transformation from Hill-Wheeler coordinates

  gamarg = gamma*pi/180D0
  a20fac = COS(gamarg)
  a22fac = SIN(gamarg)/SQRT(2D0)

! transformation from Hill-Wheeler coordinates

  falph = bbeta*a20fac
  fbeta = bbeta*a22fac
  fhexe = beta4
  WRITE (7, '(/a,f8.3,a,f8.1,a,f8.3)') &
    ' electronic input: beta=', bbeta, ' gamma=', gamma, ' beta4=', beta4

  xclust = nion*1D0
  CALL jelbak(xclust, falph, fbeta, fhexe, srms, 0)

  sum1 = 0D0
  DO ind = 1, nxyz
    sum1 = sum1 + rhojel(ind)
  END DO
  sum1 = sum1*dvol
  WRITE (6, *) 'jellium-norm=', sum1
  WRITE (7, *) 'jellium-norm=', sum1

  RETURN
END SUBROUTINE init_jellium

!-----initwf--------------------------------------------------------

SUBROUTINE initwf(psir)

! Master routine to initialize the electronic wavefunctions

  USE params
  USE util, ONLY: shiftfield, wfovlp
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: psir(kdfull2, kstate)

  INTEGER ::i, nr, nbe, nbr
  REAL(DP)::en

!----------------------------------------------------------------------

! h2m=hbar*hbar/2.0/ame
  IF (myn == 0) THEN
    IF (nspdw > nelect/2) STOP 'NSPDW must be less or equal NELECT/2'
  END IF

! prepare book-keeping for wavefunctions

  omeg = 0.25D0*h2m

  IF (ifhamdiag > 0 .AND. nstate > nelect) THEN
    WRITE (6, '(2a,2i5)') ' IFHAMDIAG>0 ONLY allowed for NSTATE=NELECT', &
      ' Presently: nstate,nelect=', nstate, nelect
    WRITE (7, '(2a,2i5)') ' IFHAMDIAG>0 ONLY allowed for NSTATE=NELECT', &
      ' Presently: nstate,nelect=', nstate, nelect
    STOP ' IFHAMDIAG>0 ONLY allowed for NSTATE=NELECT'
  END IF

  CALL ininqb(nelect, deocc, b2occ, gamocc*pi/180D0)

! initialize H.O. wavefunctions

  IF (.NOT. init_ao) THEN
    CALL initho(psir)
    CALL initho(psir)

    WRITE (7, '(a)') 'after ordering:'
    WRITE (7, '(a,2i3)') 'nstate,nelect', nstate, nelect
    WRITE (6, '(a)') 'after ordering:'
    WRITE (6, '(a,2i3)') 'nstate,nelect', nstate, nelect
    IF (ispinsep /= 0) CALL spinsep(psir)

    WRITE (7, '(a)') 'after oscillator initialization:'
    DO i = 1, nstate
      nr = nrel2abs(i)
      WRITE (7, '(a,i3,a,f12.4,a,i2,a,3i3)') 'wf ', i, &
        ' occ', occup(i), ' sp', ispin(nr), ' knots:', nq(1, nr), nq(2, nr), nq(3, nr)
    END DO

  ELSE

! optionally atomic orbital initialization

    occup = 1D0
    nstate_all = nelect
    nstate = nelect
    CALL genermowf(psir, nmaxst)
    WRITE (6, '(a)') 'after atomic-orbital initialization:'
    WRITE (7, '(a)') 'after atomic-orbital initialization:'
    DO nbr = 1, nstate
      en = SUM(psir(:, nbr)**2)*dx*dy*dz
      WRITE (6, '(2(a,i5),3(a,1pg13.5))') 'node=', myn, ', state=', nbr, &
        ', en=', en, ', occ=', occup(nbr), ', spin=', ispin(nbr)
      WRITE (7, '(2(a,i5),3(a,1pg13.5))') 'node=', myn, ', state=', nbr, &
        ', en=', en, ', occ=', occup(nbr), ', spin=', ispin(nbr)
    END DO
  END IF

! optional shift of wavefunctions

  IF (ABS(shiftwfx) + ABS(shiftwfy) + ABS(shiftwfz) > 1D-20) THEN
    DO nbe = 1, nstate
      CALL shiftfield(psir(:, nbe), shiftwfx, shiftwfy, shiftwfz)
    END DO
  END IF

! final ortho-normalization to clean up


  CALL schmidt(psir)


  RETURN
END SUBROUTINE initwf

!-----ininqb-------------------------------------------------------

SUBROUTINE ininqb(nelectin, deoccin, betain, gamin)

! Initialization of book-keeping arrays of states nq, ispin.

! Estimates Fermi energy from particle number and considers
! all states up "Fermi-energy plus deoccin" in the ordering
! of a deformed harmonic oscillator for given deformation.
! The input parameters are:
! nelect = number of electrons
! deoccin = number of osc. shells above Fermi shell
! betain = quadrupole deformation of jellium background
! gamin = triaxiality angle of jellium background
! temp = (via 'params') temperature, temp=0 cuts to occupied only

! The output goes on the occupation fields 'nq' and 'ispin'
! via module 'params'.

  USE params
  USE util, ONLY: pair
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nelectin
  REAL(DP), INTENT(IN) :: deoccin
  REAL(DP), INTENT(IN) :: betain
  REAL(DP), INTENT(IN) :: gamin

  REAL(DP) :: esp(ksttot) ! storage for s.p. energies
  REAL(DP) :: efacto ! factor to get energies from h.o.
  REAL(DP) :: q20fac ! sqrt(5/16pi)
  REAL(DP) :: cosfac, sinfac ! weightes deduced from 'gamin'
  REAL(DP) :: speact ! actual s.p. energy in loop
  INTEGER :: n ! nr. of state
  INTEGER :: noscx, noscy, noscz ! osc. nr. in each direction
  INTEGER :: nomxup, nomxdw, nmaxdiff, mspindw
  INTEGER :: i, m, isav, j, k


  REAL(DP), ALLOCATABLE :: ph(:) ! degeneracy of wavefunction, for
  REAL(DP) :: occu(ksttot)
  REAL(DP) :: ecutdw, ecutup, efrmup, efrmdw
  INTEGER :: nelup, neldw
  REAL(DP) :: gp, partnm, epstmp = 1D-5, sav
  LOGICAL, PARAMETER :: tocc = .FALSE.   ! to initialize occupation with temperatur routine

  REAL(DP), PARAMETER :: third = 1D0/3D0

!-----------------------------------------------------------------------

  ALLOCATE (ph(2*ksttot))
  IF (numspin == 2) THEN
    ph = 1D0
  ELSE IF (numspin == 1) THEN
    ph = 0D0
    ph(1:ksttot) = 2D0
  ELSE
    STOP "invalid value for input parameter NUMSPIN"
  END IF

! prepare initial parameters
! estimate of Fermi energy relies on spherical oscillator shells.
! the shell label N is determined by

! nelect,spin = N*(N+1)*(N+2)/6

! where 'nelect,spin' is the nr. of electrons for given spin.
! this relation is resolved approximately for N.

  q20fac = SQRT(5D0/(16D0*pi))
  IF (numspin == 2) THEN
    nelup = nelectin - nspdw
    neldw = nspdw
    !test
    WRITE (6, *) 'nelup', nelup
    WRITE (6, *) 'neldw', nspdw
    !test
  ELSE
    IF (MOD(nelectin, 2) == 1) &
      STOP ' nr. of electrons must be even for spin degeneracy'
    nelup = nelectin/2
  END IF
  efacto = 0.25D0/(1D0*nelectin)**third
  efrmup = (6D0*nelup)**third
  efrmup = efrmup/(1D0 - 1D0/(efrmup*efrmup))**third - 1.5D0
  ecutup = efrmup + deoccin
  nomxup = ecutup*(1D0 + 2D0*q20fac*betain) + 0.5D0
  ecutup = efacto*ecutup
  IF (numspin == 2) THEN
    IF (neldw > 0) THEN
      efrmdw = (6D0*neldw)**third
      efrmdw = efrmdw/(1D0 - 1D0/(efrmdw*efrmdw))**third - 1.5D0
    ELSE
      efrmdw = -0.00001D0
    END IF
    ecutdw = efrmdw + deoccin
    nomxdw = ecutdw*(1D0 + 2.0D0*q20fac*betain) + 0.5D0
    ecutdw = efacto*ecutdw
  END IF
  cosfac = q20fac*COS(gamin)
  sinfac = q20fac*SQRT(2D0)*SIN(gamin)
  xfac = efacto/(1D0 - betain*(cosfac - sinfac))
  yfac = efacto/(1D0 - betain*(cosfac + sinfac))
  zfac = efacto/(1D0 + 2D0*betain*cosfac)

! loop over all possible osc. triplets, determine s.p. energy
! and keep if below ecut

! select spin up states
! distribute preliminary occupation numbers

  WRITE (6, '(a,3f8.3)') ' zfac,yfac,xfac=', zfac, yfac, xfac
  WRITE (6, '(a)') ' is n nx ny nz speact'
  nmaxdiff = (ksttot - nelect)/2
  n = 0
  lb1: DO noscz = 0, nomxup
    DO noscy = 0, nomxup
      DO noscx = 0, nomxup
        speact = noscz*zfac + noscy*yfac + noscx*xfac
        IF (speact < ecutup) THEN
          n = 1 + n
          nq(1, n) = noscx
          nq(2, n) = noscy
          nq(3, n) = noscz
          ispin(n) = 1
          esp(n) = speact
          IF (n <= nelup) THEN
            occu(n) = 1D0
          ELSE
            occu(n) = 0D0
          END IF

          IF (n > ksttot) &
            STOP ' deoccin or part.number to large for given KSTTOT'

          IF (n > (nelup + nmaxdiff)) THEN
            WRITE (6, *) 'n=', n, ' greater than nelup+dnmax=', nelup + nmaxdiff
            n = n - 1 ! new from MD
            EXIT lb1
          END IF
          WRITE (6, '(a,5i5,3f10.3)') 'init states: n,...=', &
            n, ispin(n), noscx, noscy, noscz, occu(n), speact, ecutup
        END IF
      END DO
    END DO
  END DO lb1

  WRITE (6, *) 'n(nelup)', n

  IF (n < nelup) STOP ' not enough states to reach (spin-up) particle number'
  IF (numspin == 2) THEN

! select spin down states
! distribute preliminary occupation numbers

    mspindw = 0
    lb2: DO noscz = 0, nomxdw
      DO noscy = 0, nomxdw
        DO noscx = 0, nomxdw
          speact = noscz*zfac + noscy*yfac + noscx*xfac
          IF (nelup /= neldw) speact = speact*1.01D0 ! enhance for spin asymmetry
          IF (speact <= ecutdw) THEN
            n = 1 + n
            nq(1, n) = noscx
            nq(2, n) = noscy
            nq(3, n) = noscz
            ispin(n) = 2
            mspindw = mspindw + 1
            IF (mspindw <= neldw) THEN
              occu(n) = 1D0
            ELSE
              occu(n) = 0D0
            END IF
            esp(n) = efacto*speact

            IF (n > ksttot) STOP ' deocc or part.number to large for given KSTTOT'

            IF (mspindw > (neldw + nmaxdiff)) THEN
              WRITE (6, *) 'mspindw=', mspindw, &
                ' greater than neldw+dnmax=', neldw + nmaxdiff
              mspindw = mspindw - 1
              EXIT lb2
            END IF
            WRITE (6, '(a,3i5,3f10.3)') 'check spin down:', &
              ispin(n), n, mspindw, occu(n), speact, ecutdw
            WRITE (7, *)
            WRITE (7, *) 'speact', speact, 'n', n, 'isp', ispin(n)
            WRITE (7, *) 'no', noscx, noscy, noscz, zfac, yfac, xfac
          END IF
        END DO
      END DO
    END DO lb2

    WRITE (6, *) 'mspindw', mspindw
    WRITE (6, *) 'neldw', neldw

    IF (mspindw < neldw) &
      STOP ' not enough states to reach spin-down particle number'
  END IF

  nstate = n
  WRITE (7, *) 'nstate ininqb', nstate

  IF (tocc .AND. nstate > nelect) &
    CALL pair(esp, occu, ph, nelect, nstate, gp, eferm, temp, partnm, &
              90, 4, epstmp, -1, ksttot)

! reorder states in reverse order of 'occu' (largest first), spin-wise

  DO m = 1, numspin
    DO i = 1, nstate
      IF (ispin(i) .NE. m) CYCLE
      DO j = i + 1, nstate
        IF (ispin(j) .NE. m) CYCLE
        IF (occu(j) > occu(i)) THEN
          sav = occu(i)
          occu(i) = occu(j)
          occu(j) = sav
          sav = esp(i)
          esp(i) = esp(j)
          esp(j) = sav
          isav = ispin(i)
          ispin(i) = ispin(j)
          ispin(j) = isav
          DO k = 1, 3
            isav = nq(k, i)
            nq(k, i) = nq(k, j)
            nq(k, j) = isav
          END DO
        END IF
      END DO
    END DO
  END DO
  DO i = 1, nstate
    occu(i) = ph(i)*occu(i)
    WRITE (7, '(t2, 3i2, tr5, e9.3, tr5, f9.3, tr3, i2)') &
      (nq(k, i), k=1, 3), occu(i), esp(i), ispin(i)
  END DO

  DEALLOCATE (ph)

!
! initialize book-keeping fields assoxiating states with nodes

! parallel version : knode=number of nodes
! in the rest we mean absolute = relative to all the wf
! relative = relative to the wf on our node

! nstate=number of wavefunctions on our node
! (IN the nonparallel case : all the wavefunctions !)
! nrel2abs=array to get the absolute index of the wf
! (not the one on the node !)
! nabs2rel=array to get the relative index of the wf
! ( on the node !)
! nhome=home of the wf (absolute !)
! ispin=array of spins as a function of relative index:
! up and down-wfs are alternated
! occup=array of occupation as a function of relative index
! -> always 1.0 in the present code

  DO i = 1, nstate
    nrel2abs(i) = i
    nabs2rel(i) = i
    nhome(i) = 0
    occup(i) = occu(i)
  END DO
  nstate_all = nstate

! IF uneven number : better enforce by hand a zero occupation
! number for the last wf, with an even nstate !



  RETURN
END SUBROUTINE ininqb

!-----jelbak-----------------------------------------------------jelbak

SUBROUTINE jelbak(partn, alphel, betael, hexel, sqr, iturn)

! computes the jellium background such that the radius is
! given by 'radjel'*part.numb.**1/3, the surface by
! 'surjel', and the deformation is adjusted to the
! electron deformation given on 'alphel' and 'betael'

! The y40 is expressed in terms of y20_effective to align it
! automatically with the principal axes.

! This version allows three-dimensional rotation of the
! jellium background if 'iturn' is one. This serves as an excitation
! mechanism shortly before the dynamics is started.

! Input:
! partn = charge of background
! alphel = axial deformation alpha_20
! betael = triaxial deformation alpha_22
! hexel = axial hexadecapole deformation alpha_40
! iturn = initial rotation
! Output:
! sqr = jellium radius

  USE params
  USE util, ONLY: rotatevec3D
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iturn
  REAL(DP), INTENT(IN) :: partn
  REAL(DP), INTENT(IN) :: alphel
  REAL(DP), INTENT(IN) :: betael
  REAL(DP), INTENT(IN) :: hexel
  REAL(DP), INTENT(OUT) :: sqr

  INTEGER :: i1, i2, i3, ii, iter
  REAL(DP) :: vecin(3), vecout(3), vecalpha(3), vecalp(3)
  REAL(DP) :: anglex, angley, anglez, alpha, beta, thetax, thetay, thetaz, theta0
  REAL(DP) :: alphael, beta2j, gammaj
  REAL(DP) :: onetrd, d4pi, alpfac, q20fac, q22fac, q40fac
  REAL(DP) :: alphbk, betabk, hexabk, bet2EF
  REAL(DP) :: deralp, derbet, delalp, delbet
  REAL(DP) :: y2eff, y20, y22, y40, y20fac, y20obs, y22fac, y22obs, q40obs, q20red, q40red
  REAL(DP) :: radius, argum, rho0, rhoc, precis, volel
  REAL(DP) :: dpm, qoct, sqhe, sqt, sqq, sqn
  REAL(DP) :: astep
  REAL(DP) :: r, reff, rh, rr
  REAL(DP) :: x, y, z, xac, yac, zac, xact, yact, zact, xx, yy, zz
  DATA astep/0.6D0/
!----------------------------------------------------------------------

! evaluate cases according to 'irotat'

  IF (iturn == 1) THEN
! case for scissor modes
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

    WRITE (6, '(/a,3f8.3)') &
      ' jellium background rotated by anglex,angley,anglez=', vecalpha*180D0/PI
    WRITE (7, '(/a,3f8.3)') &
      ' jellium background rotated by anglex,angley,anglez=', vecalpha*180D0/PI

  ELSE IF (iturn == 0) THEN
! case for static iteration
    anglex = 0D0
    angley = 0D0
    anglez = 0D0
  END IF

! dimensionless deformation parameters
! zero = 0.0
  onetrd = 1D0/3D0
  d4pi = 1D0/(4D0*pi)
  alpfac = 4D0*pi/5D0
  q20fac = SQRT(5D0/(16D0*pi))
  q22fac = SQRT(15D0/(8D0*pi))
  q40fac = SQRT(9D0/(4D0*pi))
!----------------------------------------------------------------------

! set initial values for iteration of deformation
  alphbk = alphel
  betabk = betael
  hexabk = hexel

! effective angle 'gamma' to determine principle axis for y40

  bet2ef = (alphbk*alphbk + 2D0*betabk*betabk)
  IF (bet2ef /= 0D0) THEN
    bet2ef = 1D0/SQRT(bet2ef)
    y20fac = alphbk*bet2ef
    y20obs = y20fac*q20fac
    y22fac = betabk*bet2ef
    y22obs = y22fac*q22fac
    q40red = 7D0*SQRT(pi/9D0) ! factor for y40 from y20^2
    q20red = SQRT(5D0/pi)/7D0 ! cofactor on y20 in y40 from
    q40obs = 8D0*q40red/q40fac
  END IF

  radius = radjel*partn**onetrd
  rhoc = 3D0*d4pi/radjel**3

  argum = radius/surjel
  IF (argum < 38D0) THEN
    rho0 = rhoc*(1D0 + EXP(-argum))
  ELSE
    rho0 = rhoc
  END IF

! iterate density to give the same deformation as the electrons
  alpha = alphel
  beta = betael

  DO iter = 1, itback

! compute the distribution and accumulate its moments
    qoct = 0D0
    sqhe = 0D0
    sqt = 0D0
    sqr = 0D0
    dpm = 0D0
    sqq = 0D0
    sqn = 0D0
    ii = 0
    thetax = 0D0
    thetay = 0D0
    thetaz = 0D0
    theta0 = 0D0

    ii = 0
    DO i3 = minz, maxz
      zact = (i3 - nzsh)*dz
      DO i2 = miny, maxy
        yact = (i2 - nysh)*dy
        DO i1 = minx, maxx
          xact = (i1 - nxsh)*dx
          ii = ii + 1
          IF (abs(rotionx) + abs(rotiony) + abs(rotionz) > 0D0) THEN
            vecalp(1) = rotionx/180D0*pi
            vecalp(2) = rotiony/180D0*pi
            vecalp(3) = rotionz/180D0*pi
            vecin(1) = xact
            vecin(2) = yact
            vecin(3) = zact
            CALL rotatevec3D(vecin, vecout, vecalp)
            xac = vecout(1)
            yac = vecout(2)
            zac = vecout(3)
          ELSE
            xac = xact
            yac = yact
            zac = zact
          END IF
          IF (iturn == 1) THEN
! CALL rotxyz(xact,yact,zact,x,y,z,anglex,angley,anglez)
            vecin(1) = xac
            vecin(2) = yac
            vecin(3) = zac
            CALL rotatevec3D(vecin, vecout, vecalpha)
            x = vecout(1)
            y = vecout(2)
            z = vecout(3)
          ELSE IF (iturn == 0) THEN
            x = xac
            y = yac
            z = zac
          END IF
          xx = x*x
          yy = y*y
          zz = z*z
          rr = xx + yy + zz
          r = SQRT(rr)
          IF (rr /= zero) THEN
            y20 = q20fac*(zz + zz - xx - yy)/rr
            y22 = q22fac*(xx - yy)/rr

            IF (bet2ef == 0D0) THEN
              y40 = q40fac*(8D0*zz*zz + 3D0*xx*xx + 3D0*yy*yy - 24D0*zz*xx - &
                            24D0*zz*yy + 6D0*xx*yy)/(8D0*rr*rr)
            ELSE
              y2eff = y20fac*y20 + y22fac*y22
              y40 = q40red*(y2eff*y2eff - q20red*y2eff - d4pi)
            END IF
          ELSE
            y20 = 0D0
            y22 = 0D0
            y40 = 0D0
          END IF
          reff = radius*(1D0 + alphbk*y20 + betabk*y22 + hexabk*y40 &
                         - (alphbk*alphbk + 2D0*betabk*betabk + hexabk*hexabk)*d4pi)
          argum = (r - reff)/surjel
          IF (argum > +38D0) THEN
            rhojel(ii) = 0D0
          ELSE IF (argum < -38D0) THEN
            rhojel(ii) = rho0
          ELSE
            rhojel(ii) = rho0/(1D0 + EXP(argum))
          END IF
          rh = rhojel(ii)
          sqn = sqn + rh
          sqr = sqr + rr*rh
          y20 = (zz + zz - xx - yy)
          sqq = sqq + y20*rh
          y22 = (xx - yy)
          sqt = sqt + y22*rh
          y2eff = y20obs*y20 + y22obs*y22
          y40 = q40obs*(y2eff*y2eff - q20red*y2eff*rr - d4pi*rr*rr)
          sqhe = sqhe + y40*rh
          thetax = thetax + (yy + zz)*rh
          thetay = thetay + (xx + zz)*rh
          thetaz = thetaz + (xx + yy)*rh
        END DO
      END DO
    END DO

    volel = dx*dy*dz
    sqn = volel*sqn
    dpm = volel*dpm
    sqr = SQRT(volel*sqr/sqn)
    sqq = volel*sqq
    qoct = volel*qoct
    sqt = volel*sqt
    sqhe = volel*sqhe
    theta0 = (thetax + thetay + thetaz)/3D0
    thetax = thetax/theta0
    thetay = thetay/theta0
    thetaz = thetaz/theta0
    alpha = alpfac*q20fac*sqq/(sqn*sqr*sqr)
    beta = alpfac*0.5D0*q22fac*sqt/(sqn*sqr*sqr)
    deralp = 1D0
    derbet = 1D0
    delalp = -astep*(alpha - alphel)/deralp
    delbet = -astep*(beta - betael)/derbet
    alphbk = alphbk + delalp
    betabk = betabk + delbet
    precis = (alpha - alphel)**2 + (beta - betael)**2 + (partn - sqn)**2

    IF ((precis < endcon) .AND. (ABS(partn - sqn) < endcon)) EXIT

    radius = radius*(partn/sqn)**onetrd
    argum = radius/surjel
    IF (argum < 38D0) THEN
      rho0 = rhoc*(1D0 + EXP(-argum))
    ELSE
      rho0 = rhoc
    END IF

  END DO

  IF (iter == itback + 1) THEN
    WRITE (7, '(a/a,g11.3)') &
      ' ---> background deformation did not converge!', &
      ' residual error IN"alpha"=', precis
    STOP ' no convergence in "jelbak"'
  END IF

  WRITE (7, '(a,i4,a,f12.7,a,2(/3(a,f12.7)),5(/2(a,f12.7)))') &
    ' iteration nr.', iter, ' PRECISION=', precis, ':', &
    ' sqn =', sqn, ' dpm =', dpm, ' sqr =', sqr, &
    ' sqq =', sqq, ' sqt =', sqt, ' sqhe =', sqhe, &
    ' alpha =', alpha, ' beta =', beta, &
    ' alphbk=', alphbk - delalp, ' betabk=', betabk - delbet, &
    ' alphel=', alphael, ' betael=', betael, ' deralp=', deralp,&
    ' derbet=', derbet,  ' delalp=', delalp, ' delbet=', delbet
  WRITE (7, *) 'thetax=', thetax
  WRITE (7, *) 'thetay=', thetay
  WRITE (7, *) 'thetaz=', thetaz
  beta2j = SQRT(alpha*alpha + 2D0*beta*beta)
  gammaj = ATAN(1.4142136D0*beta/alpha)*180D0/pi
  WRITE (7, *) 'effective jellium deformations: beta2j=', beta2j, &
    ' gammaj=', gammaj

  RETURN
END SUBROUTINE jelbak

!-----initho------------------------------------------------------initho

SUBROUTINE initho(psir)

! Initializes harmonic oscillator wavefunctions

  USE params
  USE util, ONLY: wfovlp
  IMPLICIT NONE


  REAL(DP), INTENT(OUT) :: psir(kdfull2, kstate)

  REAL(DP) :: valx(nx2), valy(ny2), valz(nz2)
!~ REAL(DP), ALLOCATABLE :: phix(:)
  INTEGER :: i, ii, inx, iny, inz, ix, iy, iz, nb
  REAL(DP) :: an, bk1, bxx, bxy, bxz, hom, homx, homy, homz
  REAL(DP) :: vx, vy, vz

  REAL(DP), PARAMETER :: third = 1D0/3D0
  REAL(DP), PARAMETER :: sixth = 1D0/6D0

!----------------------------------------------------------------------

! check workspace

  an = REAL(2*nelect, DP)
  homx = omeg*an**(-third)*xfac
  homy = omeg*an**(-third)*yfac
  homz = omeg*an**(-third)*zfac
  bxx = (2D0*h2m/homx)**3
  bxy = (2D0*h2m/homy)**3
  bxz = (2D0*h2m/homz)**3
  bxx = osfac*(bxx**sixth)
  bxy = osfac*(bxy**sixth)
  bxz = osfac*(bxz**sixth)

  DO nb = 1, nstate

! number of knots

! nq is relative to the proc

! this way we implicitly parallelize the computation of phix

    inx = nq(1, nrel2abs(nb))
    iny = nq(2, nrel2abs(nb))
    inz = nq(3, nrel2abs(nb))

    DO i = 1, 3
      ipar(i, nb) = 1 - 2*MOD(nq(i, nrel2abs(nb)), 2)
    END DO

    CALL clust(inx, bxx, 0D0, xval, valx, nx2)
    CALL clust(iny, bxy, 0D0, yval, valy, ny2)
    CALL clust(inz, bxz, 0D0, zval, valz, nz2)

! composition of the factorised wave-function
! occupies only upper box (1/8 part), but is returned on 'psir'

    ii = 0
    WRITE (6, *) nx2, ny2, nz2, kstate
    WRITE (6, *) knode, nstate
    DO iz = 1, nz2
      vz = valz(iz)
      DO iy = 1, ny2
        vy = valy(iy)
        DO ix = 1, nx2
          vx = valx(ix)
          ii = ii + 1
          psir(ii, nb) = vx*vy*vz
        END DO
      END DO
    END DO

  END DO


  CALL schmidt(psir)

  RETURN
END SUBROUTINE initho

!-----clust-------------------------------------------------------

SUBROUTINE clust(IN, b, z, x, val, n1)

! Initializes 1D oscillator wavefunction

! Input:
! in number of knots
! b inverse width of harmonic oscillator solution
! z location of the nucleus resp. cluster (see initho)
! x value of grid-point
! n1 number of grid-points
! Output:
! val 1D wave-function

! The function is
! sqrt(2)^n*h(sqrt(2)*x)
! where h(x) is the Hermite polynomial.

  USE params
  IMPLICIT NONE

  INTEGER, INTENT(IN OUT) :: IN
  REAL(DP), INTENT(IN) :: b
  REAL(DP), INTENT(IN) :: z
  REAL(DP), INTENT(IN) :: x(*)
  REAL(DP), INTENT(OUT) :: val(*)
  INTEGER, INTENT(IN) :: n1

  INTEGER :: i, j
  REAL(DP) :: argum, coef, fak, tau, v
!-----------------------------------------------------------------------

  fak = 1
  DO i = 2, IN
    fak = fak*i
  END DO
  coef = 1D0/SQRT(b*1.772454D0*(2D0**IN)*fak)

  DO j = 1, n1

    tau = (x(j) - z)/b
    argum = -0.5D0*tau*tau
    IF (argum < -38D0) THEN
      val(j) = 0D0
    ELSE
      v = coef*EXP(argum)
      SELECT CASE (IN)
      CASE (0)
        val(j) = v
      CASE (1)
        val(j) = v*2*tau
      CASE (2)
        val(j) = v*(4D0*tau*tau - 2D0)
      CASE (3)
        val(j) = v*(8D0*tau**3 - 12D0*tau)
      CASE (4)
        val(j) = v*(16D0*tau**4 - 48D0*tau**2 + 12D0)
      CASE (5)
        val(j) = v*(32D0*tau**5 - 160D0*tau**3 + 120D0*tau)
      CASE (6)
        val(j) = v*8D0*(8D0*tau**6 - 60D0*tau**4 + 90D0*tau**2 - 15D0)
      CASE (7)
        val(j) = v*16D0*(8D0*tau**7 - 84D0*tau**5 + 210D0*tau**3 - 105D0*tau)
      CASE (8)
        val(j) = v*(1680D0 - 13440D0*tau**2 + 13440D0*tau**4 - &
                    3584D0*tau**6 + 256D0*tau**8)
      CASE (9)
        val(j) = v*(30240D0*tau - 80640D0*tau**3 + 483840D0*tau**5 - &
                    9216D0*tau**7 + 512D0*tau**9)
      CASE (10)
        val(j) = v*(-30240D0 + 302400D0*tau**2 - 403200D0*tau**4 + &
                    161280D0*tau**6 - 23040D0*tau**8 + 1024D0*tau**10)
      CASE (11)
        val(j) = v*(-665280D0*tau + 2217600D0*tau**3 - 1774080D0*tau**5 + &
                    506880D0*tau**7 - 56320D0*tau**9 + 2048D0*tau**11)
      CASE (12)
        val(j) = v*(665280D0 - 7983360D0*tau**2 + 13305600D0*tau**4 - &
                    7096320D0*tau**6 + 1520640D0*tau**8 - 135168D0*tau**10 + &
                    4096D0*tau**12)
      CASE (13)
        val(j) = v*(17297280D0*tau - 69189120D0*tau**3 + 69189120D0*tau**5 - &
                    26357760D0*tau**7 + 4392960D0*tau**9 - 319488D0*tau**11 + &
                    8192D0*tau**13)
      CASE (14)
        val(j) = v*(-17297280D0 + 242161920D0*tau**2 - 484323840D0*tau**4 + &
                    322882560D0*tau**6 - 92252160D0*tau**8 + 12300288D0*tau**10 - &
                    745472D0*tau**12 + 16384D0*tau**14)
      CASE (15)
        val(j) = v*(-518918400D0*tau + 2421619200D0*tau**3 - 2905943040D0*tau**5 + &
                    1383782400D0*tau**7 - 307507200D0*tau**9 + 33546240D0*tau**11 - &
                    1720320D0*tau**13 + 32768D0*tau**15)
      CASE (16)
        val(j) = v*(518918400D0 - 8302694400D0*tau**2 + 19372953600D0*tau**4 - &
                    15498362880D0*tau**6 + 5535129600D0*tau**8 - 984023040D0*tau**10 + &
                    89456640D0*tau**12 - 3932160D0*tau**14 + 65536D0*tau**16)
      CASE (17)
        val(j) = v*(17643225600D0*tau - 94097203200D0*tau**3 + 131736084480D0*tau**5 - &
                    75277762560D0*tau**7 + 20910489600D0*tau**9 - 3041525760D0*tau**11 + &
                    233963520D0*tau**13 - 8912896D0*tau**15 + 131072D0*tau**17)
      CASE (18)
        val(j) = v*(-17643225600D0 + 317578060800D0*tau**2 - 846874828800D0*tau**4 + &
                    790416506880D0*tau**6 - 338749931520D0*tau**8 + 75277762560D0*tau**10 &
                    - 9124577280D0*tau**12 + 601620480D0*tau**14 - 20054016D0*tau**16 + &
                    262144D0*tau**18)
      CASE (19)
        val(j) = v*(-670442572800D0*tau + 4022655436800D0*tau**3 - 6436248698880D0*tau**5 + &
                    4290832465920D0*tau**7 - 1430277488640D0*tau**9 + 260050452480D0*tau**11 - &
                    26671841280D0*tau**13 + 1524105216D0*tau**15 - 44826624D0*tau**17 + &
                    524288D0*tau**19)
      CASE (20)
        val(j) = v*(670442572800D0 - 13408851456000D0*tau**2 + 40226554368000D0*tau**4 - &
                    42908324659200D0*tau**6 + 21454162329600D0*tau**8 - 5721109954560D0*tau**10 + &
                    866834841600D0*tau**12 - 76205260800D0*tau**14 + 3810263040D0*tau**16 - &
                    99614720D0*tau**18 + 1048576D0*tau**20)
      CASE DEFAULT
        WRITE (6, '(a,1pg12.4)') ' wrong radial quantum number in clust: IN=', IN
        STOP 'wrong radial quantum number in clust'
      END SELECT
    END IF

  END DO

  RETURN
END SUBROUTINE clust

!-----rotxyz-----------------------------------------------------------

SUBROUTINE rotxyz(xin, yin, zin, xout, yout, zout, anglex, angley, anglez)

! Rotation in three dimensions in three separate steps,
! first, a rotation by an angle 'anglex' about the x-axis,
! second, a rotation by an angle 'angley' about the y-axis, and
! third, a rotation by an angle 'anglez' about the z-axis.
! (xin,yin,zin) is the input vector and
! (xout,yout,zout) the RESULT.

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: xin
  REAL(DP), INTENT(IN) :: yin
  REAL(DP), INTENT(IN) :: zin
  REAL(DP), INTENT(OUT) :: xout
  REAL(DP), INTENT(OUT) :: yout
  REAL(DP), INTENT(OUT) :: zout
  REAL(DP), INTENT(IN) ::anglex
  REAL(DP), INTENT(IN) ::angley
  REAL(DP), INTENT(IN) ::anglez

  REAL(DP) :: x1, y1, z1, x2, y2, z2
  REAL(DP) :: cphi, sphi
!----------------------------------------------------------------------

! rotation about x-axis, intermediate result on x1,y1,z1

  x1 = xin
  cphi = COS(anglex)
  sphi = SIN(anglex)
  y1 = yin*cphi + zin*sphi
  z1 = -yin*sphi + zin*cphi

! rotation about y-axis, intermediate result on x2,y2,z2

  y2 = y1
  cphi = COS(angley)
  sphi = SIN(angley)
  z2 = z1*cphi + x1*sphi
  x2 = -z1*sphi + x1*cphi

! rotation about z-axis, final result on xout,yout,zout

  zout = z2
  cphi = COS(anglez)
  sphi = SIN(anglez)
  xout = x2*cphi + y2*sphi
  yout = -x2*sphi + y2*cphi

  RETURN
END SUBROUTINE rotxyz

!-----spinsep-----------------------------------------------------------

SUBROUTINE spinsep(psir)

! to induce local spin current by shifting spin-up and down in
! different directions

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(IN OUT) :: psir(kdfull2, kstate)

  INTEGER :: ii, ix, iy, iz, nb
  REAL(DP) :: sgeps, x1, y1, z1

!----------------------------------------------------------------------

  IF (numspin == 2) THEN
    WRITE (*, *) ' SPINSEP invoked'
    DO nb = 1, nstate
      sgeps = (3 - 2*ispin(nb))*0.01D0
      ii = 0
      DO iz = minz, maxz
        z1 = (iz - nzsh)*dz
        DO iy = miny, maxy
          y1 = (iy - nysh)*dy
          DO ix = minx, maxx
            x1 = (ix - nxsh)*dx
            ii = ii + 1
            psir(ii, nb) = psir(ii, nb)*(1D0 + sgeps*(x1 + y1 + z1))
          END DO
        END DO
      END DO
    END DO
  END IF
  RETURN
END SUBROUTINE spinsep

!-----checkoptions-------------------------------------------------

SUBROUTINE checkoptions()

! Check consistency of compile-time options

  USE params
  IMPLICIT NONE

!------------------------------------------------------------------


  IF (.NOT. directenergy .AND. ifsicp == 4) STOP " KLI requires directenergy=.true."
  IF (directenergy .AND. idenfunc .NE. 1) STOP ' directenergy=.true. requires Perdew&Wang functional '
  IF (directenergy .AND. ifsicp == 5) &
    STOP ' directenergy=.true. not yet prepared for exact exchange '


  RETURN
END SUBROUTINE checkoptions

!-----init_grid-----------------------------------------------------

SUBROUTINE init_grid()

! Unitialize Coulomb solver, kinetic energy and other grid properties

  USE params
  USE kinetic
  USE coulsolv
  IMPLICIT NONE


  INTEGER :: nrank = 0

  INTEGER::coeff, level

  INTEGER :: ix, iy, iz
  REAL(DP) :: x1, y1, z1
!------------------------------------------------------------------

  level = nrank

  IF (level >= 1) THEN
    coeff = 2**level
    dx = dx*coeff
    dy = dy*coeff
    dz = dz*coeff
    kxbox = kxbox/coeff
    kybox = kybox/coeff
    kzbox = kzbox/coeff
    WRITE (6, *) 'level', level, dx, kxbox

    nx2 = kxbox
    ny2 = kybox
    nz2 = kzbox
    kdfull2 = nx2*ny2*nz2

    kxmax = kxbox/2 + 1; kymax = kybox/2 + 1; kzmax = kzbox/2 + 1
    nx2 = kxbox; ny2 = kybox; nz2 = kzbox; nx2fine = 2*nx2 - 1; ny2fine = 2*ny2 - 1
    nxyz = nx2*ny2*nz2; nyf = nx2; nxyf = nx2*ny2; nxyfine = nx2fine*ny2fine
    kdfull2 = kxbox*kybox*kzbox
    kdfull2fine = (2*kxbox - 1)*(2*kybox - 1)*(2*kzbox - 1)
    kdfull2fine = (2*kxbox)*(2*kybox)*(2*kzbox)
    nx = nx2/2; ny = ny2/2; nz = nz2/2
    nxfine = kxbox; nyfine = kybox; nzfine = kzbox
  END IF

! bounds of loops
  minx = 1; maxx = nx2
  miny = 1; maxy = ny2
  minz = 1; maxz = nz2
! bounds of esc. el.
  nbnx = 2; nbxx = nx2 - 1
  nbny = 2; nbxy = ny2 - 1
  nbnz = 2; nbxz = nz2 - 1
! mid-point for x,y,z-values
  nxsh = nx2/2; nysh = ny2/2; nzsh = nz2/2

  ALLOCATE (xval(nx2), yval(ny2), zval(nz2)) ! grid coordinates
  ALLOCATE (xt2(nx2), yt2(ny2), zt2(nz2)) ! coordinates**2

! check and initialize

  IF (myn == 0) THEN

    WRITE (6, *) 'kxbox,kybox,kzbox:', kxbox, kybox, kzbox
    WRITE (7, *) 'kxbox,kybox,kzbox:', kxbox, kybox, kzbox
    IF (kxbox*kybox*kzbox == 0) &
      STOP ' you must specify the box sizes KXBOX,KYBOX,KZBOX'
    IF (kstate == 0) STOP ' you must specify the maximum nr. of states KSTATE'
    IF (nabsorb > MIN(kxbox, kybox, kzbox)/4) STOP " NABSO too large"
  END IF

  dvol = dx*dy*dz

  DO ix = 1, nx2
    x1 = (ix - nxsh)*dx
    xval(ix) = x1
    xt2(ix) = x1*x1
  END DO

  DO iy = 1, ny2
    y1 = (iy - nysh)*dy
    yval(iy) = y1
    yt2(iy) = y1*y1
  END DO

  DO iz = 1, nz2
    z1 = (iz - nzsh)*dz
    zval(iz) = z1
    zt2(iz) = z1*z1
  END DO

! init kinetic energy array

  CALL init_grid_fft(dx, dy, dz, nx2, ny2, nz2, dt1, h2m)
  CALL init_coul(dx, dy, dz, nx2, ny2, nz2)


  RETURN
END SUBROUTINE init_grid

!-----ithion--------------------------------------------------------

SUBROUTINE ithion

! Initial thermalization of the ions (ONLY executed IF needed)
! we ASSIGN 0.5*kb*t per degree of freedom

  USE params
  IMPLICIT NONE

  INTEGER :: ion
  REAL(DP) :: bk, ek, ekin, pmoy, pn, vmoy, xm
!---------------------------------------------------------------------

  bk = 1D0 ! USE temperature in units of Ry
  ekin = 0.5D0*bk*tempion
  WRITE (7, '(a,f9.4)') 'wanted temp', tempion
  WRITE (7, '(a,f9.4)') 'wanted energ per degree of freedom', ekin
  xm = 0.5D0*1836.0D0*amu(np(1))*ame
  WRITE (7, *) 'warning : ithion not yet able to treat unhomogeneous systems'

! attention to masses !

  ekin = ekin/3/nion*(3*nion - 6D0)
  vmoy = SQRT(3D0*2D0*ekin/xm)
  WRITE (7, '(a,f9.4)') 'corresponding speed', vmoy
  pmoy = xm*vmoy

! we choose random directions of the speeds

  CALL spheric(pmoy, cpx, cpy, cpz, nion)
  CALL conslw(cx, cy, cz, cpx, cpy, cpz, nion)

! we rescale the speeds

  DO ion = 1, nion
    ek = cpx(ion)*cpx(ion)
    ek = ek + cpy(ion)*cpy(ion)
    ek = ek + cpz(ion)*cpz(ion)
    pn = SQRT(ek)
    cpx(ion) = cpx(ion)/pn*pmoy
    cpy(ion) = cpy(ion)/pn*pmoy
    cpz(ion) = cpz(ion)/pn*pmoy
  END DO
  CALL conslw(cx, cy, cz, cpx, cpy, cpz, nion)
  ekion = 0D0
  DO ion = 1, nion
    ek = 0D0
    ek = ek + cpx(ion)*cpx(ion)
    ek = ek + cpy(ion)*cpy(ion)
    ek = ek + cpz(ion)*cpz(ion)
    xm = 0.5D0*1836D0*amu(np(1))*ame
    ek = ek/2D0/xm
    WRITE (7, '(a,i2,a,f9.4)') 'ion', ion, 'ek=', ek
    ekion = ekion + ek
    WRITE (7, '(a,i2,3f9.4)') 'ion', ion, cpx(ion), cpy(ion), cpz(ion)
  END DO
  WRITE (7, '(a,f9.4)') 'kin.energ after renormalization', ekion
  ekion = ekion/(3*nion - 6D0)
  WRITE (7, '(a,f9.4)') 'av.kin.energ per net degree of freedom', ekion
  ekion = ekion*2D0/bk
  WRITE (7, '(a,f9.4)') 'corresponding temperature', ekion

  RETURN
END SUBROUTINE ithion

!-----transf ----------------------------------------------------------

SUBROUTINE transf

! Transform ions to the main axis of inertia

  USE params
  IMPLICIT NONE

  INTEGER :: i, i1, i2, iax, idummy, ion, ion1, ion2, icoo, iswap
  INTEGER :: j, k
  REAL(DP) :: aviner, delmx, delx, dely, dist3, gamu, zpos
  REAL(DP) :: trafo(3, 3), cenmas(3), tiner(3, 3), dminer(3)
  REAL(DP) :: pos(ng, 3), help(3)
  CHARACTER(LEN=3) :: ext

!-------------------------------------------------------------------------

  ext = outnam

  OPEN (9, FILE='center.dat', STATUS='UNKNOWN')
  OPEN (10, FILE=ext//'.pdb', STATUS='UNKNOWN')

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

! transform to center of mass coordinates

  DO ion = 1, nion
    DO icoo = 1, 3
      pos(ion, icoo) = pos(ion, icoo) - cenmas(icoo)
    END DO
  END DO

  WRITE (9, *) 'center of mass coordinates:'
  WRITE (9, *) ' '

  DO ion = 1, nion
    cx(ion) = pos(ion, 1)
    cy(ion) = pos(ion, 2)
    cz(ion) = pos(ion, 3)
    WRITE (9, '(3f10.4)') cx(ion), cy(ion), cz(ion)
  END DO
  WRITE (9, *) ' '

! compute tensor of inertia (tiner)
! cenmas = origin = fix point

  DO i1 = 1, 3
    DO i2 = 1, 3
      tiner(i1, i2) = zero
    END DO
  END DO

  DO ion = 1, nion
    tiner(1, 1) = tiner(1, 1) + pos(ion, 2)*pos(ion, 2) + pos(ion, 3)*pos(ion, 3)
    tiner(2, 2) = tiner(2, 2) + pos(ion, 3)*pos(ion, 3) + pos(ion, 1)*pos(ion, 1)
    tiner(3, 3) = tiner(3, 3) + pos(ion, 1)*pos(ion, 1) + pos(ion, 2)*pos(ion, 2)
    tiner(1, 2) = tiner(1, 2) - pos(ion, 1)*pos(ion, 2)
    tiner(1, 3) = tiner(1, 3) - pos(ion, 1)*pos(ion, 3)
    tiner(2, 3) = tiner(2, 3) - pos(ion, 3)*pos(ion, 2)
  END DO

  tiner(2, 1) = tiner(1, 2)
  tiner(3, 1) = tiner(1, 3)
  tiner(3, 2) = tiner(2, 3)

! diagonalize tensor of inertia

  CALL jacobi(tiner, 3, 3, dminer, trafo, idummy)

! transform to system of inertia
! trafo = matrix with normalized eigenvectors in columns,
! 1/trafo = trafo (t) transforms to new coordinates

  WRITE (9, '(5x,a,/,3g15.6,/)') 'Moments of inertia (x,y,z) :', &
    dminer(1), dminer(2), dminer(3)

  DO ion = 1, nion

    DO i = 1, 3
      help(i) = 0D0
      DO j = 1, 3
        help(i) = help(i) + trafo(j, i)*pos(ion, j)
      END DO
    END DO

    DO i = 1, 3
      pos(ion, i) = help(i)
    END DO

  END DO

  iswap = 1

  IF (iswap == 1) THEN

! find and fix z-axis = axis with outstanding moment of inertia

    aviner = 1D0/3D0*(dminer(1) + dminer(2) + dminer(3))
    delmx = dminer(3) - aviner
    iax = 3
! check whether x or y are outstanding
    delx = dminer(1) - aviner
    IF (ABS(delx) > ABS(delmx)) THEN
      iax = 1
      delmx = delx
    END IF
    dely = dminer(2) - aviner
    IF (ABS(dely) > ABS(delmx)) THEN
      iax = 2
      delmx = dely
    END IF
! change coordinates if necessary
    IF (iax /= 3) THEN
      DO ion = 1, nion
        zpos = pos(ion, 3)
        pos(ion, 3) = pos(ion, iax)
        pos(ion, iax) = zpos
      END DO
      WRITE (9, '(5x,a,i5,/)') 'new axis : ', iax
    END IF
    IF (delmx > zero) THEN
      WRITE (9, '(5x,a,/)') ' shape: oblate'
    ELSE
      WRITE (9, '(5x,a,/)') ' shape: prolate'
    END IF

  END IF !iswap

  WRITE (9, *) 'after transformation on the main axis of inertia:'
  WRITE (9, *) ' '
  WRITE (10, '(a10,a3)') 'COMPND ', ext

  DO ion = 1, nion
    cx(ion) = pos(ion, 1)
    cy(ion) = pos(ion, 2)
    cz(ion) = pos(ion, 3)
    WRITE (9, '(3f10.4)') cx(ion), cy(ion), cz(ion)
    WRITE (10, '(a6,a4,i1,a1,a2,a11,a1,a6,f6.3,a2, f6.3,a2,f6.3,a2,a4,a2,a4)') &
      'HETATM', ' ', ion, ' ', 'Na', ' ', '1', ' ', cx(ion), ' ', cy(ion), ' ', cz(ion), &
      ' ', '1.00', ' ', '0.00'
  END DO

  WRITE (10, '(a3)') 'END'
  WRITE (9, *) ' '

  DO ion1 = 1, nion
    k = 0
    DO ion2 = 1, nion
      IF (ion1 /= ion2 .AND. ion1 < ion2) THEN
        dist3 = SQRT((cx(ion1) - cx(ion2))**2 + (cy(ion1) - cy(ion2))**2 &
                     + (cz(ion1) - cz(ion2))**2)
        IF (k == 0) THEN
          WRITE (9, '(a16,i2,a4,i2,a2,f8.3)') &
            'distance between', ion1, ' and', ion2, ' =', dist3
        ELSE
          WRITE (9, '(a22,i2,a2,f8.3)') ' and', ion2, ' =', dist3
        END IF
        k = k + 1
      END IF
    END DO
  END DO

  CLOSE (9)
  CLOSE (10)
  RETURN
END SUBROUTINE transf

!---------------------------------------------------------------------

SUBROUTINE jacobi(a, n, np, d, v, nrot)
  USE params, ONLY: DP
  IMPLICIT NONE

! Compute eigenvalues 'd' and eigenvectors of a symmetric matrix
! 'a' and returns the matrix 'v' whose columns contain the
! normalized eigenvectors of 'a'
! Adapted from numerical recipes (p. 456):

  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: np
  REAL(DP), INTENT(IN OUT) :: a(np, np)
  REAL(DP), INTENT(OUT) :: d(np)
  REAL(DP), INTENT(OUT) :: v(np, np)
  INTEGER, INTENT(OUT) :: nrot

  INTEGER, PARAMETER :: nmax = 500
  INTEGER :: i, ip, iq, j
  REAL(DP) :: c, g, h, s, sm, t, tau, theta, tresh, b(nmax), z(nmax)

  DO ip = 1, n
    DO iq = 1, n
      v(ip, iq) = 0D0
    END DO
    v(ip, ip) = 1D0
  END DO

  DO ip = 1, n
    b(ip) = a(ip, ip)
    d(ip) = b(ip)
    z(ip) = 0D0
  END DO

  nrot = 0
  DO i = 1, 50
    sm = 0D0
    DO ip = 1, n - 1
      DO iq = ip + 1, n
        sm = sm + ABS(a(ip, iq))
      END DO
    END DO
    IF (sm == 0D0) RETURN
    IF (i < 4) THEN
      tresh = 0.2D0*sm/n**2
    ELSE
      tresh = 0D0
    END IF
    DO ip = 1, n - 1
      DO iq = ip + 1, n
        g = 1D2*ABS(a(ip, iq))
        IF ((i > 4) .AND. (ABS(d(ip)) + g == ABS(d(ip))) &
            .AND. (ABS(d(iq)) + g == ABS(d(iq)))) THEN
          a(ip, iq) = 0D0
        ELSE IF (ABS(a(ip, iq)) > tresh) THEN
          h = d(iq) - d(ip)
          IF (ABS(h) + g == ABS(h)) THEN
            t = a(ip, iq)/h
          ELSE
            theta = 0.5D0*h/a(ip, iq)
            t = 1D0/(ABS(theta) + SQRT(1D0 + theta**2))
            IF (theta < 0D0) t = -t
          END IF
          c = 1./SQRT(1 + t**2)
          s = t*c
          tau = s/(1.+c)
          h = t*a(ip, iq)
          z(ip) = z(ip) - h
          z(iq) = z(iq) + h
          d(ip) = d(ip) - h
          d(iq) = d(iq) + h
          a(ip, iq) = 0D0
          DO j = 1, ip - 1
            g = a(j, ip)
            h = a(j, iq)
            a(j, ip) = g - s*(h + g*tau)
            a(j, iq) = h + s*(g - h*tau)
          END DO
          DO j = ip + 1, iq - 1
            g = a(ip, j)
            h = a(j, iq)
            a(ip, j) = g - s*(h + g*tau)
            a(j, iq) = h + s*(g - h*tau)
          END DO
          DO j = iq + 1, n
            g = a(ip, j)
            h = a(iq, j)
            a(ip, j) = g - s*(h + g*tau)
            a(iq, j) = h + s*(g - h*tau)
          END DO
          DO j = 1, n
            g = v(j, ip)
            h = v(j, iq)
            v(j, ip) = g - s*(h + g*tau)
            v(j, iq) = h + s*(g - h*tau)
          END DO
          nrot = nrot + 1
        END IF
      END DO
    END DO
    DO ip = 1, n
      b(ip) = b(ip) + z(ip)
      d(ip) = b(ip)
      z(ip) = 0D0
    END DO
  END DO

  STOP 'too many iterations in jacobi'

  RETURN
END SUBROUTINE jacobi


SUBROUTINE genermowf(psiom, nmxst)

! Generates orbital molecular wave function.
!
! Output:
! psiom = wavefunctions
! nmxst = maximum quantum state for an atom
! other parameters via module 'params'

  USE params
  IMPLICIT NONE

  REAL(DP), INTENT(OUT) :: psiom(kdfull2, kstate)
  INTEGER, INTENT(OUT) :: nmxst(1:ng) ! maximum nr. for each atom

  INTEGER :: i, ind, ion, ix, iy, iz, is
  INTEGER :: nadd, nadded, natlevel, nbr, ncy, ncycle, ndiff, nmaxval, nstx, nsty, nstz, numstate, numlocstate
  REAL(DP) :: x1, y1, z1
  INTEGER, ALLOCATABLE :: nactst(:) ! keep track of actual atomic state
  INTEGER, ALLOCATABLE :: ncount(:) ! keep track of actual atomic state
  INTEGER, ALLOCATABLE :: ispinsav(:) ! keep track of spin
  INTEGER :: nnodes(1:3, 1:10)
  DATA nnodes/0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, &
    1, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 2/

!---------------------------------------------------------------------

! reset wavefunctions 'psiom'

  ALLOCATE (nactst(1:ng))
  ALLOCATE (ncount(0:knode))
  nactst = 0
  ncount = 0

  WRITE (6, *) ' GENERLCGO entered. NION,NSTATE=', nion, nstate
  WRITE (6, *) ' CH(ion):', (ch(np(ion)), ion=1, nion)
  FLUSH (6)

  DO nbr = 1, nstate
    DO i = 1, kdfull2
      psiom(i, nbr) = 0D0
    END DO
  END DO

! book-keeping of number of states

  nmaxval = 0
  DO ion = 1, nion
    nmxst(ion) = ch(np(ion))
    nmaxval = nmaxval + ch(np(ion))
  END DO
  ndiff = nmaxval - nelect
  ncycle = ndiff/nion + 1
  IF (ndiff > 0) THEN
    nadd = +1
  ELSE
    nadd = -1
  END IF

  IF (ndiff /= 0) THEN
    nadded = 0
    cycleloop: DO ncy = 1, ncycle
      DO ion = nion, 1, -1
        nmxst(ion) = nmxst(ion) - nadd
        nadded = nadded + nadd
        IF (nadded == ndiff) EXIT cycleloop
      END DO
    END DO cycleloop
  END IF
  WRITE (6, '(a,5i5)') &
    ' ndiff,nmaxval,nstate,ncycle,nadd=', ndiff, nmaxval, nstate, ncycle, nadd
  WRITE (6, '(a,100i3)') ' nmxst:', (nmxst(ion), ion=1, nion)
  WRITE (6, '(a,100i3)') ' ipol:', (ipol(ion), ion=1, nion)
  WRITE (6, *) ' ipolcheck:', (ipol(ion), ion=1, nion)

! loop through ions and fill electron states successively

  numstate = 0
  numlocstate = 0 ! in case of parallelism
  ionloop: DO ion = 1, nion
    DO natlevel = 1, nmxst(ion)
      IF (numstate == ksttot) THEN
        WRITE (6, *) 'numstate = ksttot, increase kstate or nproc'
        WRITE (7, *) 'numstate = ksttot, increase kstate or nproc'
        EXIT ionloop
      END IF
      numstate = 1 + numstate
      ncount(nhome(numstate)) = ncount(nhome(numstate)) + 1
      IF (numspin == 2) THEN
        ispin(numstate) = 2 - MOD(numstate, 2)
        IF (ipol(ion) .eq. -1) ispin(numstate) = 2 - mod(numstate + 1, 2)
        nactst(ion) = nactst(ion) + MOD(natlevel, 2)
      ELSE
        nactst(ion) = nactst(ion) + 1
      END IF
      IF (nhome(numstate) == myn) THEN
        numlocstate = numlocstate + 1

        WRITE (6, *) 'ion,natlev,numst,ispin', &
          ion, natlevel, numstate, ispin(numstate)
        WRITE (7, *) 'ksttot,nstate,ion,natlev,numst,ispin', &
          ksttot, nstate, ion, natlevel, numstate, ispin(numstate)
        IF (nactst(ion) > 10) STOP 'GENERMOWF: too high atomic level'
        IF (numlocstate > nstate) THEN
          WRITE (6, *) 'numlocstate > nstate, increase kstate or nproc'
          WRITE (7, *) 'numlocstate > nstate, increase kstate or nproc'
          EXIT ionloop
        END IF
! select nodes
        nstx = nnodes(initord(1, ion), nactst(ion))
        nsty = nnodes(initord(2, ion), nactst(ion))
        nstz = nnodes(initord(3, ion), nactst(ion))
! compute raw wf
        ind = 0
        DO iz = 1, nz2
          z1 = ((iz - nzsh)*dz - cz(ion))/radini(ion)
          DO iy = 1, ny2
            y1 = ((iy - nysh)*dy - cy(ion))/radini(ion)
            DO ix = 1, nx2
              x1 = ((ix - nxsh)*dx - cx(ion))/radini(ion)
              ind = ind + 1
              psiom(ind, numlocstate) = (x1**nstx - 0.5D0*MAX(0, nstx - 1)) &
                                        *(y1**nsty - 0.5D0*MAX(0, nsty - 1))&
                                        *(z1**nstz - 0.5D0*MAX(0, nstz - 1)) &
                                        *EXP(-(x1*x1 + y1*y1 + z1*z1)*0.5D0)
            END DO
          END DO
        END DO

        WRITE (6, '(a,i3,a,a,5i3,1pg13.5)') &
          ' electron state nr.', numstate, ': ', &
          ' ion,nstx,nsty,nstz=', ion, nstx, nsty, nstz, ispin(numstate), &
          SUM(psiom(:, numlocstate)**2)*dvol
      END IF
    END DO
  END DO ionloop

  DEALLOCATE (nactst)

! sort wavefunctions in blocks of spins
  OPEN (UNIT=19,STATUS='SCRATCH',FORM='unformatted')
  ALLOCATE(ispinsav(1:numstate))
  ispinsav = ispin(1:numstate)

  DO i=1,numstate
    WRITE (19) ispin(i),occup(i),psiom(:,i)
  END DO
  REWIND(19)

  i = 0
  DO is=1,2
    DO nbr = 1,numstate
      IF(ispinsav(nbr) == is) THEN
        i = i + 1
        READ (19) ispin(i),occup(i),psiom(:,i)
        IF(ispin(i) .NE. is) STOP 'GENERMOWF: error in spin sorting'
      ELSE
        READ (19) ind                        ! read a dummy
      END IF
    END DO
    REWIND(19)
  END DO

  IF(i .NE. numstate) &
    STOP 'GENERMOWF: wrong nr. of wavefunctions after spin sorting'

  CLOSE (19)

  FLUSH (6)

  RETURN
END SUBROUTINE genermowf

