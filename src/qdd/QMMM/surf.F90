
INTEGER, PARAMETER :: klmax = 2 ! ????
REAL(DP), PARAMETER :: varelcorelimit = 1D-8 ! limiting value for V_elArcore [Ry]

!INTEGER,PARAMETER :: kfermi=kxbox*17330 ! length of interpolation table
INTEGER, PARAMETER :: kfermi = 64*17330 ! length of interpolation table
! uses here fixed length -- check!

! surface part:
INTEGER, PARAMETER :: ngpar = 3000 ! max NUMBER of surface particules
INTEGER :: maxpar
INTEGER, PARAMETER :: maxnlayers = 3
INTEGER, PARAMETER :: kdsub = (2*nxsg + 1)*(2*nysg + 1)*(2*nzsg + 1)

! switches
INTEGER :: ivdw = 0, isurf = 0

INTEGER :: idielec = 0 ! switch TO dielectricum
REAL(DP) :: xdielec = 0D0 ! x below which dielectric zone is activated
REAL(DP) :: epsdi = 1D0 ! dielectric constant IN the dielectric zone

! sub-grid for V_elArcore

INTEGER :: nsg_arelcore, imobtmp
INTEGER :: isubgcenter(3*ngpar), ioutofbox(3*ngpar)
REAL(DP) :: c_dipmod = 3.387D0
REAL(DP), ALLOCATABLE :: rfieldtmp(:)
REAL(DP) :: rscaltmp, chgtmp, sigtmp

REAL(DP) :: dte, told, testeps, temperature
INTEGER :: iregc(ngpar), iregk(ngpar), irege(ngpar)
INTEGER :: imobc(ngpar), imobe(ngpar), imobk(ngpar)
INTEGER :: itypc(ngpar), itype(ngpar), itypk(ngpar)
INTEGER :: ippointer1(3*ngpar), ippointer2(3*ngpar)
INTEGER :: ipointer(3*ngpar), ipointerfix(3*ngpar), ipointermob(3*ngpar)
INTEGER :: isrtyp(5, 5) = 0, isrtypall = 0
INTEGER :: ndfull2, nmob, nmobc, nmobe, nmobk, ionthefly
INTEGER :: iposevry, ienevry, iobsevry, ithermevry !,ifile
INTEGER :: icool, icoolevry, icoolsteps, itindex
!INTEGER :: jenergion = 10

REAL(DP) :: xc(ngpar), yc(ngpar), zc(ngpar) ! POSITION of fixed cores / gsm particles
REAL(DP) :: xe(ngpar), ye(ngpar), ze(ngpar) ! POSITION of fixed shells
REAL(DP) :: xk(ngpar), yk(ngpar), zk(ngpar) ! POSITION of fixed cations
REAL(DP) :: xcinit(ngpar), ycinit(ngpar), zcinit(ngpar)
REAL(DP) :: xeinit(ngpar), yeinit(ngpar), zeinit(ngpar)
REAL(DP) :: xkinit(ngpar), ykinit(ngpar), zkinit(ngpar)
REAL(DP) :: xcold(ngpar), ycold(ngpar), zcold(ngpar)
REAL(DP) :: xeold(ngpar), yeold(ngpar), zeold(ngpar)
REAL(DP) :: xkold(ngpar), ykold(ngpar), zkold(ngpar)
REAL(DP) :: pxc(ngpar), pyc(ngpar), pzc(ngpar)
REAL(DP) :: pxe(ngpar), pye(ngpar), pze(ngpar), pxk(ngpar), pyk(ngpar), pzk(ngpar)
REAL(DP) :: pxcold(ngpar), pycold(ngpar), pzcold(ngpar)
REAL(DP) :: pxeold(ngpar), pyeold(ngpar), pzeold(ngpar)
REAL(DP) :: pxkold(ngpar), pykold(ngpar), pzkold(ngpar)
REAL(DP) :: fxc(ngpar), fyc(ngpar), fzc(ngpar)
REAL(DP) :: fxcold(ngpar), fycold(ngpar), fzcold(ngpar)
REAL(DP) :: fxkold(ngpar), fykold(ngpar), fzkold(ngpar)
REAL(DP) :: fxk(ngpar), fyk(ngpar), fzk(ngpar), fxe(ngpar), fye(ngpar), fze(ngpar)
REAL(DP) :: fxeold(ngpar), fyeold(ngpar), fzeold(ngpar)
REAL(DP) :: chgc(ngpar), chge(ngpar), chgk(ngpar)
REAL(DP) :: forcesx(3*ngpar), forcesy(3*ngpar), forcesz(3*ngpar)
REAL(DP), ALLOCATABLE :: fxt(:), fyt(:), fzt(:)
REAL(DP) :: rtransvec(3, 3)
REAL(DP) :: celldipole(3), cellquad(3, 3), cellmult(klmax, 2*klmax + 1)
REAL(DP) :: shiftx = 0D0, shifty = 0D0, shiftz = 0D0, fermiac, fermibc, fermicc
REAL(DP) :: fermiae, fermibe, fermice, fermiak, fermibk, fermick
REAL(DP) :: fermia2c, fermib2c, fermic2c
REAL(DP) :: fermia2e, fermib2e, fermic2e, fermia2k, fermib2k, fermic2k
REAL(DP) :: rlattvec(3), xatom(ngpar, 3)
INTEGER :: ipsort(ngpar)
INTEGER :: natoms, ncells1x, ncells1y, ncells1z, ncells2x, ncells2y, ncells2z
INTEGER :: ibh = 0, cbh = 0

! probably obsolete ??
REAL(DP) :: fermia, fermib, fermic

! for vdW
REAL(DP) :: evdw, espvdw, esub
REAL(DP), ALLOCATABLE :: potvdw(:), frho(:, :)

REAL(DP), ALLOCATABLE :: potstat(:), potesfixed(:)
REAL(DP), ALLOCATABLE :: potesmob(:), phim(:), phimv(:), phimd(:)

REAL(DP) :: forx, fory, forz
REAL(DP) :: sigkv, sigvv, sigkk, sigcc, sigck, sigcv, bkv, bvv, bkk, bcc, bck, bcv
REAL(DP) :: srcut, crcut, rmob, ckv, cvv, ckk, ccc, cck, ccv, ccn, sigcn, bcn, enerinfty
REAL(DP) :: sigkn, bkn, ckn, ecn, fcn, ekn, fkn, sigkn2, bkn2, ekn2, fkn2
REAL(DP) :: ccn2, ccn3, ccn4, ccn5, ccn6, ccn7, ccn8, ccnd, dcn, ccn9, ccn10
REAL(DP) :: ckn2, ckn3, ckn4, ckn5, ckn6, ckn7, ckn8, cknd, dkn
REAL(DP) :: ckn9, ckn10, ccncu2, ccncu3, ccncu4, ccncu5, ccncu6
REAL(DP) :: ccncu7, ccncu8, ccncu9, ccncu10, ccncud, ckncu2, ckncu3, ckncu4, ckncu5, ckncu6
REAL(DP) :: ckncu7, ckncu8, ckncu9, ckncu10, ckncud
REAL(DP) :: ccncul2, ccncul3, ccncul4, ccncul5, ccncul6
REAL(DP) :: ccncul7, ccncul8, ccncul9, ccncul10, ccnculd
REAL(DP) :: ckncul2, ckncul3, ckncul4, ckncul5, ckncul6
REAL(DP) :: ckncul7, ckncul8, ckncul9, ckncul10, cknculd
REAL(DP) :: ccel6, ccel8, ccel10, ckel6, ckel8, ckel10
REAL(DP) :: sigmac = 1.01116D0, sigmav = 1.01116D0, sigmak = 1D0, sigmacc, sigmacv
REAL(DP) :: sigmack, sigmavv, sigmakk, sigmakv, zsurf
REAL(DP) :: chgcore, chgval, chgkat
REAL(DP) :: cspr = 6.7585128D0 ! =6.119**2*e2/11.08
REAL(DP) :: me = 0.00238562D0 ! =4.38/1836.0
REAL(DP) :: mkat = 1D1, mion = 39.95D0
INTEGER :: nrare = 0 ! NUMBER of raregas atoms
INTEGER :: nc = 0 ! NUMBER of fixed cores IN substrate (O cores IN MgO(001) ?)
INTEGER :: nk = 0 ! NUMBER of fixed cations IN substrate (Mg cations IN MGO(001) ?)
INTEGER :: ne = 0 ! NUMBER of fixed shells IN substrate
INTEGER :: nlayers, nside, iararlj, ipog
INTEGER :: iforcecl2co = 0

REAL(DP) :: surftemp = 0D0, scaledist = 1D0
REAL(DP) :: chgc0, chge0, chgk0
REAL(DP) :: unfixcrad = -1D0, unfixerad = -1D0, unfixkrad = -1D0
REAL(DP) :: unfixclateralrad = -1D0, unfixelateralrad = -1D0, unfixklateralrad = -1D0
REAL(DP) :: unfixclateralrady = -1D0, unfixelateralrady = -1D0, unfixklateralrady = -1D0
REAL(DP) :: unfixclateralradx = -1D0, unfixelateralradx = -1D0, unfixklateralradx = -1D0
REAL(DP) :: fixcbelow = -1D6, fixebelow = -1D6, fixkbelow = -1D6
REAL(DP) :: fixcbelowy = -1D6, fixebelowy = -1D6, fixkbelowy = -1D6
REAL(DP) :: fixcbelowx = -1D6, fixebelowx = -1D6, fixkbelowx = -1D6
REAL(DP) :: r1x, r1y, r1z, r2ax, r2ay, r2az, r2bx, r2by, r2bz
INTEGER :: ipotfixed, isystem = 2, jsavesurf = 10000000
INTEGER :: iaxis = 0, ilayerc(ngpar), ilayere(ngpar), ilayerk(ngpar)
INTEGER :: nlay, nmovc, nmove, nmovk, nimovc, nimove, nimovk
INTEGER :: ncr1, ncr2a, ncr2b, ner1, ner2a, ner2b, nkr1, nkr2a, nkr2b
INTEGER :: ixyzcm(ngpar), ixyzem(ngpar), ixyzkm(ngpar), iunfixall = 0
INTEGER :: ifixall, iusecell = 0
INTEGER :: iuselast = 0, iprintonlyifmob = 0
INTEGER :: ixyzcim(ngpar), ixyzeim(ngpar), ixyzkim(ngpar)
INTEGER :: ifmdshort, ifadiadip

REAL(DP) :: epots, ecoul, ekinion, ekinkat, ekinel, eshort, energyg
REAL(DP) :: ekincsurf, ekinesurf, ekinksurf, ekinsurf, ttemp, ttempc, ttempe, ttempk

!
! interpolation tables for radial potentials
! ---> eats up too much space !
! ==> check need and find more efficient ways
REAL(DP) :: varelcore0(kfermi), varelcore1(kfermi)
REAL(DP) :: vfermi0(kfermi), vfermi1(kfermi)

