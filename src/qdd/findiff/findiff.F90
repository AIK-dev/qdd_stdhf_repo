
MODULE kinetic
  USE params
  IMPLICIT NONE

CONTAINS
  SUBROUTINE kinprop(q0, q1)
!
! INTERFACE TO finite-difference propagator
!
    USE params
    COMPLEX(DP), INTENT(IN OUT) :: q1(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: q2(kdfull2) ! workspace, obsolete here

    CALL d3mixpropag(q1, dt1)

  END SUBROUTINE kinprop

  SUBROUTINE calc_ekin(psin, ekinout)

! calculates kinetic energy for single particle state with
! COMPLEX wavefunction 'psin'.

    USE params
    USE util, ONLY: wfovlp
    IMPLICIT NONE

    COMPLEX(DP), INTENT(IN) :: psin(kdfull2)
    REAL(DP), INTENT(OUT) :: ekinout
    INTEGER :: ii
    REAL(DP) :: sum0, sumk

    COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: psi2

!------------------------------------------------------------------

    ALLOCATE (psi2(kdfull2))

! exp.value of kinetic energy

    CALL ckin3d(psin, psi2)
!sum0 = 0D0
!sumk = 0D0
!DO ii=1,nxyz
! sumk = REAL(psin(ii),DP)*REAL(psi2(ii),DP) + AIMAG(psin(ii))*AIMAG(psi2(ii)) &
! + sumk
! sum0 = REAL(psin(ii),DP)*REAL(psin(ii),DP) + AIMAG(psin(ii))*AIMAG(psin(ii)) &
! + sum0
!END DO
    ekinout = REAL(wfovlp(psin, psi2), DP)

  END SUBROUTINE calc_ekin

#if(findiff)

!-----RKIN3D_3R------------------------------------------------------------

  SUBROUTINE rkin3d(psi, dxpsi)

!RKIN3D COMPUTES THE LAPLACIAN OF 'PSI' AND PUT THE RESULT IN 'DXPSI'
!'PSI' AND 'DXPSI' ARE REAL NUMBERS
!THIS COMPUTATION IS BASED ON THE FINITES DIFFERENCES METHOD IN THREE POINTS

!-----------------------------------------------------------------------

! D2PSI(X,Y,Z) D2PSI(X,Y,Z) D2PSI(X,Y,Z)
! DXPSI(X,Y,Z) = ------------ + ------------ + ------------
! DX2 DY2 DZ2

!-----------------------------------------------------------------------

    REAL(DP), INTENT(IN) :: psi(kdfull2)
    REAL(DP), INTENT(OUT) :: dxpsi(kdfull2)

!-----------------------------------------------------------------------

    INTEGER :: ix, iy, iz ! LOOP INDICES
    INTEGER :: ista, iend, ind ! INDICE OF SELECTION OF THE DIRECTION

!-----------------------------------------------------------------------

! RESET ACCUMULATOR
    dxpsi = 0

! KINETIC ENERGY IN X DIRECTION
    ind = 1
    DO iz = 1, kzbox
      DO iy = 1, kybox
        ista = (iz - 1)*kxbox*kybox + (iy - 1)*kxbox + 1
        iend = (iz - 1)*kxbox*kybox + (iy - 1)*kxbox + kxbox
        CALL rkin1d_3r(psi(ista:iend:ind), dx, kxbox, dxpsi(ista:iend:ind))
      END DO
    END DO

! KINETIC ENERGY IN Y DIRECTION
    ind = kxbox
    DO iz = 1, kzbox
      DO ix = 1, kxbox
        ista = (iz - 1)*kxbox*kybox + (ix - 1) + 1
        iend = (iz - 1)*kxbox*kybox + (ix - 1) + kxbox*(kybox - 1) + 1
        CALL rkin1d_3r(psi(ista:iend:ind), dy, kybox, dxpsi(ista:iend:ind))
      END DO
    END DO

! KINETIC ENERGY IN Z DIRECTION
    ind = kxbox*kybox
    DO iy = 1, kybox
      DO ix = 1, kxbox
        ista = (iy - 1)*kxbox + (ix - 1) + 1
        iend = (iy - 1)*kxbox + (ix - 1) + kxbox*kybox*(kzbox - 1) + 1
        CALL rkin1d_3r(psi(ista:iend:ind), dz, kzbox, dxpsi(ista:iend:ind))
      END DO
    END DO

    RETURN
  END SUBROUTINE rkin3d

!-----D2_3R1D ----------------------------------------------------------

  SUBROUTINE rkin1d_3r(psi, deltax, nmax, dxpsi)

! PSI = INPUT ARRAY IN ONE DIRECTION (A 3D ARRAY IN CALLING ROUTINE)
! DELTAX = MESH SPACING
! NMAX = NUMBER OF MESH POINT IN GIVEN DIRECTION
! DXPSI = OUTPUT CONTAINING 2. DERIVATIVE + INPUT VALUE

!-----------------------------------------------------------------------

!RKIN3D_3R COMPUTES THE SECOND DERIVATIVE OF 'PSI' IN ONE DIRECTION AND PUT THE RESULT IN 'DXPSI'
!'PSI' AND 'DXPSI' ARE REAL NUMBERS
!THIS COMPUTATION IS BASED ON THE FINITES DIFFERENCES METHOD IN THREE POINTS

!-----------------------------------------------------------------------
!
! DXPSI(X) = PSI(X+H) - 2*PSI(X) + PSI(X-H)
! ------------------------------
! (H^2)
!
!-----------------------------------------------------------------------

    REAL(DP), INTENT(IN) :: psi(*)
    REAL(DP), INTENT(IN) :: deltax
    INTEGER, INTENT(IN) :: nmax
    REAL(DP), INTENT(INOUT) :: dxpsi(*)

!-----------------------------------------------------------------------

    INTEGER :: i ! LOOP INDICE
    REAL(DP):: d2i ! FINITE DIFFERENCE COEFFICIENT

!-----------------------------------------------------------------------
    d2i = -one/(deltax*deltax)
!-----------------------------------------------------------------------

    dxpsi(1) = d2i*(psi(2) - 2*psi(1)) + dxpsi(1)

    DO i = 2, nmax - 1
      dxpsi(i) = dxpsi(i) + d2i*(psi(i + 1) + psi(i - 1) - 2.0D0*psi(i))
    END DO

    dxpsi(nmax) = d2i*(psi(nmax - 1) - 2*psi(nmax)) + dxpsi(nmax)

    RETURN
  END SUBROUTINE rkin1d_3r

!-----CKIN3D_3R------------------------------------------------------------

  SUBROUTINE ckin3d(psi, dxpsi)

!CKIN3D COMPUTES THE LAPLACIAN OF 'PSI' AND PUT THE RESULT IN 'DXPSI'
!'PSI' AND 'DXPSI' ARE COMPLEX NUMBERS
!THIS COMPUTATION IS BASED ON THE FINITES DIFFERENCES METHOD IN THREE POINTS

!-----------------------------------------------------------------------

! D2PSI(X,Y,Z) D2PSI(X,Y,Z) D2PSI(X,Y,Z)
! DXPSI(X,Y,Z) = ------------ + ------------ + ------------
! DX2 DY2 DZ2

!-----------------------------------------------------------------------

    COMPLEX(DP), INTENT(IN) :: psi(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: dxpsi(kdfull2)

!-----------------------------------------------------------------------

    INTEGER :: ix, iy, iz ! LOOP INDICES
    INTEGER :: ista, iend, ind ! INDICE OF SELECTION OF THE DIRECTION

!-----------------------------------------------------------------------

! RESET ACCUMULATOR
    dxpsi = 0

! KINETIC ENERGY IN X DIRECTION
    ind = 1
    DO iz = 1, kzbox
      DO iy = 1, kybox
        ista = (iz - 1)*kxbox*kybox + (iy - 1)*kxbox + 1
        iend = (iz - 1)*kxbox*kybox + (iy - 1)*kxbox + kxbox
        CALL ckin1d_3r(psi(ista:iend:ind), dx, kxbox, dxpsi(ista:iend:ind))
      END DO
    END DO

! KINETIC ENERGY IN Y DIRECTION
    ind = kxbox
    DO iz = 1, kzbox
      DO ix = 1, kxbox
        ista = (iz - 1)*kxbox*kybox + (ix - 1) + 1
        iend = (iz - 1)*kxbox*kybox + (ix - 1) + kxbox*(kybox - 1) + 1
        CALL ckin1d_3r(psi(ista:iend:ind), dy, kybox, dxpsi(ista:iend:ind))
      END DO
    END DO

! KINETIC ENERGY IN Z DIRECTION
    ind = kxbox*kybox
    DO iy = 1, kybox
      DO ix = 1, kxbox
        ista = (iy - 1)*kxbox + (ix - 1) + 1
        iend = (iy - 1)*kxbox + (ix - 1) + kxbox*kybox*(kzbox - 1) + 1
        CALL ckin1d_3r(psi(ista:iend:ind), dz, kzbox, dxpsi(ista:iend:ind))
      END DO
    END DO

    RETURN
  END SUBROUTINE ckin3d
!-----d2_3r1D ----------------------------------------------------------

  SUBROUTINE ckin1d_3r(psi, deltax, nmax, dxpsi)

! PSI = INPUT ARRAY IN ONE DIRECTION (A 3D ARRAY IN CALLING ROUTINE)
! DELTAX = MESH SPACING
! NMAX = NUMBER OF MESH POINT IN GIVEN DIRECTION
! DXPSI = OUTPUT CONTAINING 2. DERIVATIVE + INPUT VALUE

!-----------------------------------------------------------------------

!CKIN3D_3R COMPUTES THE SECOND DERIVATIVE OF 'PSI' IN ONE DIRECTION AND PUT THE RESULT IN 'DXPSI'
!'PSI' AND 'DXPSI' ARE COMPLEX NUMBERS
!THIS COMPUTATION IS BASED ON THE FINITES DIFFERENCES METHOD IN THREE POINTS

!-----------------------------------------------------------------------
!
! DXPSI(X) = PSI(X+H) - 2*PSI(X) + PSI(X-H)
! ------------------------------
! (H^2)
!
!-----------------------------------------------------------------------

    COMPLEX(DP), INTENT(IN) :: psi(*)
    REAL(DP), INTENT(IN) :: deltax
    INTEGER, INTENT(IN) :: nmax
    COMPLEX(DP), INTENT(INOUT) :: dxpsi(*)

    INTEGER :: i ! LOOP INDICE
    REAL(DP):: d2i ! FINITE DIFFERENCE COEFFICIENT

!-----------------------------------------------------------------------
    d2i = -one/(deltax*deltax)
!-----------------------------------------------------------------------

    dxpsi(1) = d2i*(psi(2) - 2*psi(1)) + dxpsi(1)

    DO i = 2, nmax - 1
      dxpsi(i) = dxpsi(i) + d2i*(psi(i + 1) + psi(i - 1) - 2.0D0*psi(i))
    END DO

    dxpsi(nmax) = d2i*(psi(nmax - 1) - 2*psi(nmax)) + dxpsi(nmax)

    RETURN
  END SUBROUTINE ckin1d_3r

!-----d3mixpropag_3----------------------------------------------------

  SUBROUTINE d3mixpropag(psi, deltim)

!USE params
! include 'pot3D.inc'

! Propagates T OPERATOR for mixed step:

! 1 - 0.5 * i * deltim * T^
! neupsi = ------------------------ * psi
! 1 + 0.5 * i * deltim * T^

! IN 3 Dimensions with 3- point PRECISION
! with boundary conditions zero at doth ends (reflecting)

! . The "Laplacian"
! OPERATOR is here dimensionless with diagonal element -2
! and off-diagonals 1.

! list parameters:

    COMPLEX(DP), INTENT(IN OUT) :: psi(kdfull2)
    REAL(DP), INTENT(IN OUT) :: deltim

! ! old FUNCTION
! COMPLEX neupsi (MINX:MAXX, MINY:MAXY, MINZ:MAXZ)
! ! new FUNCTION

! Timestep
! variables:

!-----------------------------------------------------------------------

    INTEGER :: ix, iy, iz ! LOOP INDICES
    INTEGER :: ista, iend, ind ! INDICE OF SELECTION OF THE DIRECTION

!-----------------------------------------------------------------------

! WRITE(6,*) ' enter D3MIXSTEP: norm=',wfnorm(psi)
    DO iz = minz, maxz
      DO iy = miny, maxy
        ! CALL kinprop_1d3(psi(minx,iy,iz),2*maxx+1,1,dx,deltim )
        ista = (iz - 1)*kxbox*kybox + (iy - 1)*kxbox + 1
        iend = (iz - 1)*kxbox*kybox + (iy - 1)*kxbox + kxbox
        CALL kinprop_1d3(psi(ista:iend:ind), dx, kxbox, deltim)
      END DO
    END DO ! x- direction
! WRITE(6,*) ' after x D3MIXSTEP: norm=',wfnorm(psi)

    DO iz = minz, maxz
      DO ix = minx, maxx
        !CALL kinprop_1d3(psi(ix,miny,iz),2*maxy+1,2*maxx+1, dy,deltim)
        ista = (iz - 1)*kxbox*kybox + (ix - 1) + 1
        iend = (iz - 1)*kxbox*kybox + (ix - 1) + kxbox*(kybox - 1) + 1
        CALL kinprop_1d3(psi(ista:iend:ind), dy, kybox, deltim)
      END DO
    END DO ! y- direction
! WRITE(6,*) ' after y D3MIXSTEP: norm=',wfnorm(psi)

    DO iy = miny, maxy
      DO ix = minx, maxx
        !CALL kinprop_1d3(psi(ix,iy,minz),2*maxz+1,(2*maxx+1)*(2*maxy+1),dz,deltim)
        ista = (iy - 1)*kxbox + (ix - 1) + 1
        iend = (iy - 1)*kxbox + (ix - 1) + kxbox*kybox*(kzbox - 1) + 1
        CALL kinprop_1d3(psi(ista:iend:ind), dz, kzbox, deltim)
      END DO
    END DO ! z- direction
! WRITE(6,*) ' after z D3MIXSTEP: norm=',wfnorm(psi)

    RETURN
  END SUBROUTINE d3mixpropag

!-----inv3p_ini-------------------------------------------------------------

  SUBROUTINE inv3p_ini(deltim)

! initialize Array of fixed values for invers_3r.
! ONLY for same delta IN all directions, yet

!USE params
! include 'pot3D.inc'

    REAL(DP), INTENT(IN OUT) :: deltim
    COMPLEX(DP), ALLOCATABLE :: invnum(:)
!inverse diagonal elements
    COMPLEX(DP) :: diag
!diagonalelement of forward-matrix
!COMMON /invnum3c/ invnum,diag

    COMPLEX(DP) :: fac !diagonal elements
    INTEGER :: n !loop

    ALLOCATE (invnum(2*(maxx + maxy + maxz) + 3))

!---------------------------------------------------------------------------

    WRITE (*, *) 'Inititalisierung der Koeffizienten fuer fin. Differenzen'

    IF (dx /= dy .OR. dy /= dz .OR. dz /= dy) &
      STOP 'ONLY same gridstep IN every direction'

    fac = CMPLX(-2.0D0, 2.0D0*dx*dx/deltim, DP)
    invnum(1) = one/fac

    DO n = 2, 2*(maxx + maxy + maxz) + 3
      invnum(n) = one/(fac - invnum(n - 1))
    END DO

    diag = CMPLX(2.0D0, 2.0D0*dx*dx/deltim, DP)

! WRITE(6,*) ' init: diag,invnum=',diag,invnum
    DEALLOCATE (invnum)
    RETURN
  END SUBROUTINE inv3p_ini

!-----kinprop_1d3--------------------------------------------------------------

  SUBROUTINE kinprop_1d3(psi, deltax, ndiml, deltim)

! 1
! Solve = ------------------------- * Rhs
! 1 + 0.5 * i * deltim * T^

! T^ is kinetic energy OPERATOR gives IN 3 point PRECISION
! (finite differenc)
! with boundary conditions ZERO at both ends (reflecting)

!USE params

    COMPLEX(DP), INTENT(IN OUT) :: psi(*)
    INTEGER, INTENT(IN) :: ndiml
    REAL(DP), INTENT(IN) :: deltax
    REAL(DP), INTENT(IN) :: deltim
    COMPLEX(DP), ALLOCATABLE :: invnum(:)
    COMPLEX(DP) :: diag
!COMMON /invnum3c/ invnum,diag

!INTEGER :: ! array SIZE
!COMPLEX(DP) :: ! wave FUNCTION TO be propagated
    COMPLEX(DP) :: solve
    COMPLEX(DP) :: psip, psim, psi0
!REAL(DP):: ! mesh SIZE

!INTEGER :: ! Work array SIZE
    INTEGER:: ndimx

    COMPLEX(DP), ALLOCATABLE:: reff(:) ! effective r.h.s.
!c $ offdiag, ! constant off-diagonal el.
!c $ offdiag2 ! squared off-diagonal el.

    INTEGER :: i, n, ninc, inc ! loop index

!----------------------------------------------------------------------------

    ALLOCATE (invnum(2*(maxx + maxy + maxz) + 3), reff(ndimx))

! DIRECT step (1-i dt/2 H) on psi

    ndimx = 2*(maxx + maxy + maxz) + 3
    ninc = 1
    psi0 = psi(ninc)
    ninc = ninc + inc
    psip = psi(ninc)
    reff(1) = diag*psi0 - psip

    DO i = 2, ndiml - 1
      psim = psi0
      psi0 = psip
      ninc = ninc + inc
      psip = psi(ninc)
      reff(i) = diag*psi0 - (psip + psim)
    END DO

    psim = psi0
    psi0 = psip
    reff(ndiml) = diag*psi0 - psim

! forward loop (eliminate lower off- diagonal elements)

    DO n = 2, ndiml
      reff(n) = reff(n) - reff(n - 1)*invnum(n - 1)
    END DO

! backward loop: eliminate upper off-diagonal and solve

    ninc = (ndiml - 1)*inc + 1
    psi(ninc) = reff(ndiml)*invnum(ndiml)
    solve = psi(ninc)
    DO n = ndiml - 1, 1, -1
      ninc = ninc - inc
      solve = (reff(n) - solve)*invnum(n)
      psi(ninc) = solve
    END DO

    DEALLOCATE (invnum, reff)

    RETURN
  END SUBROUTINE kinprop_1d3
#endif

#if(numerov)
!-----ckin3D_5r------------------------------------------------------------

  SUBROUTINE rkin3d(psi, dxpsi)

!RKIN3D COMPUTES THE LAPLACIAN OF 'PSI' AND PUT THE RESULT IN 'DXPSI'
!'PSI' AND 'DXPSI' ARE REAL NUMBERS
!THIS COMPUTATION IS BASED ON THE FINITES DIFFERENCES METHOD IN FIVE POINTS

!-----------------------------------------------------------------------

! D2PSI(X,Y,Z) D2PSI(X,Y,Z) D2PSI(X,Y,Z)
! DXPSI(X,Y,Z) = ------------ + ------------ + ------------
! DX2 DY2 DZ2

!-----------------------------------------------------------------------

    REAL(DP), INTENT(IN) :: psi(kdfull2)
    REAL(DP), INTENT(OUT) :: dxpsi(kdfull2)

!-----------------------------------------------------------------------

    INTEGER :: ix, iy, iz ! LOOP INDICES
    INTEGER :: ista, iend, ind ! INDICE OF SELECTION OF THE DIRECTION

!-----------------------------------------------------------------------

! RESET ACCUMULATOR
    dxpsi = 0

! KINETIC ENERGY IN X DIRECTION
    ind = 1
    DO iz = 1, kzbox
      DO iy = 1, kybox
        ista = (iz - 1)*kxbox*kybox + (iy - 1)*kxbox + 1
        iend = (iz - 1)*kxbox*kybox + (iy - 1)*kxbox + kxbox
        CALL rkin1d_5r(psi(ista:iend:ind), dx, kxbox, dxpsi(ista:iend:ind))
      END DO
    END DO

! KINETIC ENERGY IN Y DIRECTION
    ind = kxbox
    DO iz = 1, kzbox
      DO ix = 1, kxbox
        ista = (iz - 1)*kxbox*kybox + (ix - 1) + 1
        iend = (iz - 1)*kxbox*kybox + (ix - 1) + kxbox*(kybox - 1) + 1
        CALL rkin1d_5r(psi(ista:iend:ind), dy, kybox, dxpsi(ista:iend:ind))
      END DO
    END DO

! KINETIC ENERGY IN Z DIRECTION
    ind = kxbox*kybox
    DO iy = 1, kybox
      DO ix = 1, kxbox
        ista = (iy - 1)*kxbox + (ix - 1) + 1
        iend = (iy - 1)*kxbox + (ix - 1) + kxbox*kybox*(kzbox - 1) + 1
        CALL rkin1d_5r(psi(ista:iend:ind), dz, kzbox, dxpsi(ista:iend:ind))
      END DO
    END DO

    RETURN
  END SUBROUTINE rkin3d
!-----CKIN1D_5R ----------------------------------------------------------

  SUBROUTINE rkin1d_5r(psi, deltax, nmax, dxpsi)

! PSI = INPUT ARRAY IN ONE DIRECTION (A 3D ARRAY IN CALLING ROUTINE)
! DELTAX = MESH SPACING
! NMAX = NUMBER OF MESH POINT IN GIVEN DIRECTION
! DXPSI = OUTPUT CONTAINING 2. DERIVATIVE + INPUT VALUE

!-----------------------------------------------------------------------

!RKIN3D_5R COMPUTES THE SECOND DERIVATIVE OF 'PSI' IN ONE DIRECTION AND PUT THE RESULT IN 'DXPSI'
!'PSI' AND 'DXPSI' ARE REAL NUMBERS
!THIS COMPUTATION IS BASED ON THE FINITES DIFFERENCES METHOD IN FIVE POINTS

!-----------------------------------------------------------------------
!
! DXPSI(X) = - PSI(X-2H) + 16*PSI(X-H) - 30*PSI(X) + 16*PSI(X+H) - PSI(X+2H)
! ---------------------------------------------------------------
! 12*(H^2)
!
!-----------------------------------------------------------------------

    REAL(DP), INTENT(IN) :: psi(*)
    REAL(DP), INTENT(IN) :: deltax
    INTEGER, INTENT(IN) :: nmax
    REAL(DP), INTENT(INOUT) :: dxpsi(*)

!-----------------------------------------------------------------------
    INTEGER :: i ! LOOP INDICE
    REAL(DP):: d2i ! FINITE DIFFERENCE COEFFICIENT
!-----------------------------------------------------------------------
    d2i = -one/(deltax*deltax)/12.0D0
!-----------------------------------------------------------------------
    i = 1
    dxpsi(i) = dxpsi(i) + d2i*(-psi(i + 2) + 16.0D0*psi(i + 1) - 30.0D0*psi(i))

    i = 2
    dxpsi(i) = dxpsi(i) + d2i*(-psi(i + 2) + 16.0D0*psi(i + 1) - 30.0D0*psi(i) + 16.0D0*psi(i - 1))

    DO i = 3, nmax - 2
      dxpsi(i) = dxpsi(i) + d2i*(-psi(i + 2) + 16.0D0*psi(i + 1) - 30.0D0*psi(i) + 16.0D0*psi(i - 1) - psi(i - 2))
    END DO

    i = nmax - 1
    dxpsi(i) = dxpsi(i) + d2i*(16.0D0*psi(i + 1) - 30.0D0*psi(i) + 16.0D0*psi(i - 1) - psi(i - 2))

    i = nmax
    dxpsi(i) = dxpsi(i) + d2i*(-30.0D0*psi(i) + 16.0D0*psi(i - 1) - psi(i - 2))

    RETURN
  END SUBROUTINE rkin1d_5r
!-----rkin3D_5r------------------------------------------------------------

  SUBROUTINE ckin3d(psi, dxpsi)

!CKIN3D COMPUTES THE LAPLACIAN OF 'PSI' AND PUT THE RESULT IN 'DXPSI'
!'PSI' AND 'DXPSI' ARE COMPLEX NUMBERS
!THIS COMPUTATION IS BASED ON THE FINITES DIFFERENCES METHOD IN FIVE POINTS

!-----------------------------------------------------------------------

! D2PSI(X,Y,Z) D2PSI(X,Y,Z) D2PSI(X,Y,Z)
! DXPSI(X,Y,Z) = ------------ + ------------ + ------------
! DX2 DY2 DZ2

!-----------------------------------------------------------------------

    COMPLEX(DP), INTENT(IN) :: psi(kdfull2)
    COMPLEX(DP), INTENT(OUT) :: dxpsi(kdfull2)

!-----------------------------------------------------------------------

    INTEGER :: ix, iy, iz ! LOOP INDICES
    INTEGER :: ista, iend, ind ! INDICE OF SELECTION OF THE DIRECTION

!-----------------------------------------------------------------------

! RESET ACCUMULATOR
    dxpsi = 0

! KINETIC ENERGY IN X DIRECTION
    ind = 1
    DO iz = 1, kzbox
      DO iy = 1, kybox
        ista = (iz - 1)*kxbox*kybox + (iy - 1)*kxbox + 1
        iend = (iz - 1)*kxbox*kybox + (iy - 1)*kxbox + kxbox
        CALL ckin1d_5r(psi(ista:iend:ind), dx, kxbox, dxpsi(ista:iend:ind))
      END DO
    END DO

! KINETIC ENERGY IN Y DIRECTION
    ind = kxbox
    DO iz = 1, kzbox
      DO ix = 1, kxbox
        ista = (iz - 1)*kxbox*kybox + (ix - 1) + 1
        iend = (iz - 1)*kxbox*kybox + (ix - 1) + kxbox*(kybox - 1) + 1
        CALL ckin1d_5r(psi(ista:iend:ind), dy, kybox, dxpsi(ista:iend:ind))
      END DO
    END DO

! KINETIC ENERGY IN Z DIRECTION
    ind = kxbox*kybox
    DO iy = 1, kybox
      DO ix = 1, kxbox
        ista = (iy - 1)*kxbox + (ix - 1) + 1
        iend = (iy - 1)*kxbox + (ix - 1) + kxbox*kybox*(kzbox - 1) + 1
        CALL ckin1d_5r(psi(ista:iend:ind), dz, kzbox, dxpsi(ista:iend:ind))
      END DO
    END DO

    RETURN
  END SUBROUTINE ckin3d
!-----rkin1D_5r ----------------------------------------------------------

  SUBROUTINE ckin1d_5r(psi, deltax, nmax, dxpsi)

! PSI = INPUT ARRAY IN ONE DIRECTION (A 3D ARRAY IN CALLING ROUTINE)
! DELTAX = MESH SPACING
! NMAX = NUMBER OF MESH POINT IN GIVEN DIRECTION
! DXPSI = OUTPUT CONTAINING 2. DERIVATIVE + INPUT VALUE

!-----------------------------------------------------------------------

!CKIN3D_5R COMPUTES THE SECOND DERIVATIVE OF 'PSI' IN ONE DIRECTION AND PUT THE RESULT IN 'DXPSI'
!'PSI' AND 'DXPSI' ARE COMPLEX NUMBERS
!THIS COMPUTATION IS BASED ON THE FINITES DIFFERENCES METHOD IN FIVE POINTS

!-----------------------------------------------------------------------
!
! DXPSI(X) = - PSI(X-2H) + 16*PSI(X-H) - 30*PSI(X) + 16*PSI(X+H) - PSI(X+2H)
! ---------------------------------------------------------------
! 12*(H^2)
!
!-----------------------------------------------------------------------

    COMPLEX(DP), INTENT(IN) :: psi(*)
    REAL(DP), INTENT(IN) :: deltax
    INTEGER, INTENT(IN) :: nmax
    COMPLEX(DP), INTENT(OUT) :: dxpsi(*)
!-----------------------------------------------------------------------
    INTEGER :: i ! LOOP INDICE
    REAL(DP):: d2i ! FINITE DIFFERENCE COEFFICIENT
!-----------------------------------------------------------------------
    d2i = -one/(deltax*deltax)/12.0D0
!-----------------------------------------------------------------------
    i = 1
    dxpsi(i) = dxpsi(i) + d2i*(-psi(i + 2) + 16.0D0*psi(i + 1) - 30.0D0*psi(i))

    i = 2
    dxpsi(i) = dxpsi(i) + d2i*(-psi(i + 2) + 16.0D0*psi(i + 1) - 30.0D0*psi(i) + 16.0D0*psi(i - 1))

    DO i = 3, nmax - 2
      dxpsi(i) = dxpsi(i) + d2i*(-psi(i + 2) + 16.0D0*psi(i + 1) - 30.0D0*psi(i) + 16.0D0*psi(i - 1) - psi(i - 2))
    END DO

    i = nmax - 1
    dxpsi(i) = dxpsi(i) + d2i*(16.0D0*psi(i + 1) - 30.0D0*psi(i) + 16.0D0*psi(i - 1) - psi(i - 2))

    i = nmax
    dxpsi(i) = dxpsi(i) + d2i*(-30.0D0*psi(i) + 16.0D0*psi(i - 1) - psi(i - 2))

    RETURN
  END SUBROUTINE ckin1d_5r

!-----d3mixpropag_5----------------------------------------------------

  SUBROUTINE d3mixpropag(psi, deltim)

!USE params
! include 'pot3D.inc'

! Propagates T OPERATOR for mixed step:

! 1 - 0.5 * i * deltim * T^
! neupsi = ------------------------ * psi
! 1 + 0.5 * i * deltim * T^

! IN 3 Dimensions with 5- point PRECISION (Numerov)
! with boundary conditions zero at doth ends (reflecting)

! . The "Laplacian"
! OPERATOR is here dimensionless with diagonal element -2
! and off-diagonals 1.

! list parameters:

    COMPLEX(DP), INTENT(IN OUT) :: psi(minx:maxx, miny:maxy, minz:maxz)
!COMPLEX(DP), INTENT(IN OUT) :: neupsi(minx:maxx, miny:maxy,minz:maxz)
    REAL(DP), INTENT(IN OUT) :: deltim

! ! old FUNCTION

! ! new FUNCTION

! Timestep

! variables:

    INTEGER :: ix, iy, iz
! loop indices

!-----------------------------------------------------------------------

    DO iz = minz, maxz
      DO iy = miny, maxy
        CALL kinprop_1d5(psi(minx, iy, iz), 2*maxx + 1, 1, dx, deltim)
      END DO
    END DO ! x- direction

    DO iz = minz, maxz
      DO ix = minx, maxx
        CALL kinprop_1d5(psi(ix, miny, iz), 2*maxy + 1, 2*maxx + 1, dy, deltim)
      END DO
    END DO ! y- direction

    DO iy = miny, maxy
      DO ix = minx, maxx
        CALL kinprop_1d5(psi(ix, iy, minz), 2*maxz + 1, &
                         (2*maxx + 1)*(2*maxy + 1), dz, deltim)
      END DO
    END DO ! z- direction

    RETURN
  END SUBROUTINE d3mixpropag

!-----inv5p_ini-------------------------------------------------------------

  SUBROUTINE inv5p_ini(deltim)

! initialize Array of fixed values for kinprop__1d5.
! ONLY for same delta IN all directions, yet

!USE params
! include 'pot3D.inc'

    REAL(DP), INTENT(IN) :: deltim
!COMPLEX(DP) :: invnum(2*(maxx+maxy+maxz)+3)
    COMPLEX(DP), ALLOCATABLE :: invnum(:)

!inverse diagonal Elements
    COMPLEX(DP) :: diag, offd
!matrix-elements for forward-matrix
!COMMON /invnum5c/ invnum,diag,offd

    REAL(DP):: fac1, fac2 !diagonal elements
    COMPLEX(DP) :: fac
    INTEGER :: n, dim !loop

!---------------------------------------------------------------------------
    ALLOCATE (invnum(2*(maxx + maxy + maxz) + 3))

    WRITE (*, *) 'Inititalisierung der Koeffizienten fuer Numerov'

    IF (dx /= dy .OR. dy /= dz .OR. dz /= dy) &
      STOP 'ONLY same gridstep IN every direction'

    fac1 = (-72.0D0*deltim*deltim + 10.0D0*dx**4)/(36.0D0*deltim*deltim + dx**4)
    fac2 = 72.0D0*deltim*dx*dx/(36.0D0*deltim*deltim + dx**4)
    fac = CMPLX(fac1, fac2, DP)

    invnum(1) = one/fac

    DO n = 2, 2*(maxx + maxy + maxz) + 3
      invnum(n) = one/(fac - invnum(n - 1))
    END DO

    fac1 = (72.0D0*deltim*deltim + 10.0D0*dx**4)/(36.0D0*deltim*deltim + dx**4)
    fac2 = 48.0D0*deltim*dx*dx/(36.0D0*deltim*deltim + dx**4)
    diag = CMPLX(fac1, fac2, DP)

    fac1 = (-36.0D0*deltim*deltim + dx**4)/(36.0D0*deltim*deltim + dx**4)
    fac2 = 12.0D0*deltim*dx*dx/(dx**4 + 36.0D0*deltim*deltim)
    offd = CMPLX(fac1, fac2, DP)

    DEALLOCATE (invnum)

    RETURN
  END SUBROUTINE inv5p_ini

!-----kinprop_1d5--------------------------------------------------------------

  SUBROUTINE kinprop_1d5(psi, ndiml, inc, deltax, deltim)

! 1 - 0.5 * i * deltim * T
! Solve = ------------------------- * Rhs
! 1 + 0.5 * i * deltim * T^

! T^ is kinetic energy OPERATOR gives IN 3 point PRECISION
! (finite differenc)
! with boundary conditions ZERO at both ends (reflecting)

!USE params
! include 'work3D.inc'

    COMPLEX(DP), INTENT(IN OUT) :: psi(ndiml)
    INTEGER, INTENT(IN) :: ndiml
    INTEGER, INTENT(IN) :: inc
    REAL(DP), INTENT(IN OUT) :: deltax
    REAL(DP), INTENT(IN OUT) :: deltim
!COMPLEX(DP) :: invnum(2*(maxx+maxy+maxz)+3)
    COMPLEX(DP), ALLOCATABLE :: invnum(:), reff(:)
    COMPLEX(DP) :: diag, offd
!COMMON /invnum5c/ invnum,diag,offd

!INTEGER :: ! array SIZE
!COMPLEX(DP) :: ! wave FUNCTION TO be propagated
    COMPLEX(DP) :: solve
    COMPLEX(DP) :: psip, psim, psi0
!REAL(DP):: ! mesh SIZE

!INTEGER :: ! Work array SIZE
    INTEGER :: ndimx

!COMPLEX(DP) :: reff (ndimx) ! effective r.h.s.
!COMPLEX,ALLOCATABLE :: invnum

!c $ offdiag, ! constant off-diagonal el.
!c $ offdiag2 ! squared off-diagonal el.

    INTEGER :: i, n, ninc ! loop index
    REAL(DP):: fac1, fac2, fac3, fac4

!----------------------------------------------------------------------------

! DIRECT step (1-i dt/2 H) on psi
    ndimx = 2*(maxx + maxy + maxz) + 3
    ALLOCATE (invnum(ndimx), reff(ndimx))

    ninc = 1
    psi0 = psi(ninc)
    ninc = ninc + inc
    psip = psi(ninc)
    reff(1) = diag*psi0 + offd*psip

    DO i = 2, ndiml - 1
      psim = psi0
      psi0 = psip
      ninc = ninc + inc
      psip = psi(ninc)
      reff(i) = diag*psi0 + offd*(psip + psim)
    END DO

    psim = psi0
    psi0 = psip
    reff(ndiml) = diag*psi0 + offd*psim

! forward loop (eliminate lower off- diagonal elements)

    DO n = 2, ndiml
      reff(n) = reff(n) - reff(n - 1)*invnum(n - 1)
    END DO

! backward loop: eliminate upper off-diagonal and solve

    ninc = (ndiml - 1)*inc + 1
    psi(ninc) = reff(ndiml)*invnum(ndiml)
    solve = psi(ninc)
    DO n = ndiml - 1, 1, -1
      ninc = ninc - inc
      solve = (reff(n) - solve)*invnum(n)
      psi(ninc) = solve
    END DO

    DEALLOCATE (invnum)

    RETURN
  END SUBROUTINE kinprop_1d5
#endif
END MODULE kinetic

