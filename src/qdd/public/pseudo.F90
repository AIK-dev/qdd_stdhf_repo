! This file is part of 'params.F90'. It contains arrays and settings
! for the PsP of a selection of elements.

INTEGER, PARAMETER :: kpsp = 99, kpspm = 0

REAL(DP) :: ch(kpspm:kpsp) = 0D0, amu(kpspm:kpsp) = 0D0, &
            cc1(kpspm:kpsp) = 0D0, cc2(kpspm:kpsp) = 0D0
REAL(DP) :: crloc(kpspm:kpsp) = 0D0, crs(kpspm:kpsp) = 0D0, chs(kpspm:kpsp) = 0D0, &
            chg1(kpspm:kpsp) = 0D0, chg2(kpspm:kpsp) = 0D0
REAL(DP) :: sgm1(kpspm:kpsp) = 0D0, sgm2(kpspm:kpsp) = 0D0
REAL(DP) :: dr1(kpspm:kpsp), dr2(kpspm:kpsp)
REAL(DP) :: prho1(kpspm:kpsp), prho2(kpspm:kpsp)
REAL(DP) :: r0g(kpspm:kpsp) = 0D0, r1g(kpspm:kpsp) = 0D0, r2g(kpspm:kpsp) = 0D0, &
            h0_11g(kpspm:kpsp) = 0D0, h0_22g(kpspm:kpsp) = 0D0
REAL(DP) :: h0_33g(kpspm:kpsp) = 0D0, h1_11g(kpspm:kpsp) = 0D0, &
            h1_22g(kpspm:kpsp) = 0D0
REAL(DP) :: h2_11g(kpspm:kpsp) = 0D0, radiong(kpspm:kpsp) = 0D0
REAL(DP) :: h0_12g(kpspm:kpsp) = -1D20 ! default signal "automatic"

INTEGER, PARAMETER :: knl = 64000 ! storage for PsP projectors

LOGICAL :: tnonlocany ! flag to invoke non-local part
LOGICAL, ALLOCATABLE :: tnonloc(:) ! flag to invoke non-local part
INTEGER, ALLOCATABLE :: ifin(:), icount(:, :), np(:)
INTEGER, ALLOCATABLE :: ifinfine(:), icountfine(:, :)
INTEGER, ALLOCATABLE :: icountfinesp(:), icountsp(:)

REAL(DP), ALLOCATABLE :: p0_1(:, :), p0_2(:, :), p1_1(:, :), p1_1x(:, :)
REAL(DP), ALLOCATABLE :: p1_1y(:, :), p1_1z(:, :)
REAL(DP), ALLOCATABLE :: p0_3(:, :), p1_2(:, :), p1_2x(:, :)
REAL(DP), ALLOCATABLE :: p1_2y(:, :), p1_2z(:, :)
REAL(DP), ALLOCATABLE :: p2_1(:, :), p2_xy(:, :), p2_xz(:, :)
REAL(DP), ALLOCATABLE :: p2_yz(:, :), p2_xy2(:, :), p2_z2(:, :)
REAL(DP), ALLOCATABLE :: p0_1fine(:, :), p0_2fine(:, :), p1_1fine(:, :), p1_1xfine(:, :)
REAL(DP), ALLOCATABLE :: p1_1yfine(:, :), p1_1zfine(:, :)
REAL(DP), ALLOCATABLE :: p0_3fine(:, :), p1_2fine(:, :), p1_2xfine(:, :)
REAL(DP), ALLOCATABLE :: p1_2yfine(:, :), p1_2zfine(:, :)
REAL(DP), ALLOCATABLE :: p2_1fine(:, :), p2_xyfine(:, :), p2_xzfine(:, :)
REAL(DP), ALLOCATABLE :: p2_yzfine(:, :), p2_xy2fine(:, :), p2_z2fine(:, :)
