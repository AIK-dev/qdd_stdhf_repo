MODULE fishpack
! FISHPK12 FROM PORTLIB 12/30/83
CONTAINS

  SUBROUTINE hw3crt(XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M, MBDCND, BDYS,&
  & BDYF, ZS, ZF, N, NBDCND, BDZS, BDZF, lambda, LDIMF,&
  & MDIMF, F, PERTRB, IERROR, W)
!
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * *
! * F I S H P A K *
! * *
! * *
! * A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF *
! * *
! * SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS *
! * *
! * (VERSION 3.1 , OCTOBER 1980) *
! * *
! * BY *
! * *
! * JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET *
! * *
! * OF *
! * *
! * THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH *
! * *
! * BOULDER, COLORADO (80307) U.S.A. *
! * *
! * WHICH IS SPONSORED BY *
! * *
! * THE NATIONAL SCIENCE FOUNDATION *
! * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! * * * * * * * * * PURPOSE * * * * * * * * * * * * * * * * * *
!
! SUBROUTINE HW3CRT SOLVES THE STANDARD SEVEN-POINT FINITE
! DIFFERENCE APPROXIMATION TO THE HELMHOLTZ EQUATION IN CARTESIAN
! COORDINATES:
!
! (D/DX)(DU/DX) + (D/DY)(DU/DY) + (D/DZ)(DU/DZ)
!
! + LAMBDA*U = F(X,Y,Z) .
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! * * * * * * * * PARAMETER DESCRIPTION * * * * * * * * * *
!
!
! * * * * * * ON INPUT * * * * * *
!
! XS,XF
! THE RANGE OF X, I.E. XS .LE. X .LE. XF .
! XS MUST BE LESS THAN XF.
!
! L
! THE NUMBER OF PANELS INTO WHICH THE INTERVAL (XS,XF) IS
! SUBDIVIDED. HENCE, THERE WILL BE L+1 GRID POINTS IN THE
! X-DIRECTION GIVEN BY X(I) = XS+(I-1)DX FOR I=1,2,...,L+1,
! WHERE DX = (XF-XS)/L IS THE PANEL WIDTH. L MUST BE AT
! LEAST 5 .
!
! LBDCND
! INDICATES THE TYPE OF BOUNDARY CONDITIONS AT X = XS AND X = XF.
!
! = 0 IF THE SOLUTION IS PERIODIC IN X, I.E.
! U(L+I,J,K) = U(I,J,K).
! = 1 IF THE SOLUTION IS SPECIFIED AT X = XS AND X = XF.
! = 2 IF THE SOLUTION IS SPECIFIED AT X = XS AND THE DERIVATIVE
! OF THE SOLUTION WITH RESPECT TO X IS SPECIFIED AT X = XF.
! = 3 IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO X IS
! SPECIFIED AT X = XS AND X = XF.
! = 4 IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO X IS
! SPECIFIED AT X = XS AND THE SOLUTION IS SPECIFIED AT X=XF.
!
! BDXS
! A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE
! DERIVATIVE OF THE SOLUTION WITH RESPECT TO X AT X = XS.
! WHEN LBDCND = 3 OR 4,
!
! BDXS(J,K) = (D/DX)U(XS,Y(J),Z(K)), J=1,2,...,M+1,
! K=1,2,...,N+1.
!
! WHEN LBDCND HAS ANY OTHER VALUE, BDXS IS A DUMMY VARIABLE.
! BDXS MUST BE DIMENSIONED AT LEAST (M+1)*(N+1).
!
! BDXF
! A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE
! DERIVATIVE OF THE SOLUTION WITH RESPECT TO X AT X = XF.
! WHEN LBDCND = 2 OR 3,
!
! BDXF(J,K) = (D/DX)U(XF,Y(J),Z(K)), J=1,2,...,M+1,
! K=1,2,...,N+1.
!
! WHEN LBDCND HAS ANY OTHER VALUE, BDXF IS A DUMMY VARIABLE.
! BDXF MUST BE DIMENSIONED AT LEAST (M+1)*(N+1).
!
! YS,YF
! THE RANGE OF Y, I.E. YS .LE. Y .LE. YF.
! YS MUST BE LESS THAN YF.
!
! M
! THE NUMBER OF PANELS INTO WHICH THE INTERVAL (YS,YF) IS
! SUBDIVIDED. HENCE, THERE WILL BE M+1 GRID POINTS IN THE
! Y-DIRECTION GIVEN BY Y(J) = YS+(J-1)DY FOR J=1,2,...,M+1,
! WHERE DY = (YF-YS)/M IS THE PANEL WIDTH. M MUST BE AT
! LEAST 5 .
!
! MBDCND
! INDICATES THE TYPE OF BOUNDARY CONDITIONS AT Y = YS AND Y = YF.
!
! = 0 IF THE SOLUTION IS PERIODIC IN Y, I.E.
! U(I,M+J,K) = U(I,J,K).
! = 1 IF THE SOLUTION IS SPECIFIED AT Y = YS AND Y = YF.
! = 2 IF THE SOLUTION IS SPECIFIED AT Y = YS AND THE DERIVATIVE
! OF THE SOLUTION WITH RESPECT TO Y IS SPECIFIED AT Y = YF.
! = 3 IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y IS
! SPECIFIED AT Y = YS AND Y = YF.
! = 4 IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y IS
! SPECIFIED AT Y = YS AND THE SOLUTION IS SPECIFIED AT Y=YF.
!
! BDYS
! A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE
! DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y AT Y = YS.
! WHEN MBDCND = 3 OR 4,
!
! BDYS(I,K) = (D/DY)U(X(I),YS,Z(K)), I=1,2,...,L+1,
! K=1,2,...,N+1.
!
! WHEN MBDCND HAS ANY OTHER VALUE, BDYS IS A DUMMY VARIABLE.
! BDYS MUST BE DIMENSIONED AT LEAST (L+1)*(N+1).
!
! BDYF
! A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE
! DERIVATIVE OF THE SOLUTION WITH RESPECT TO Y AT Y = YF.
! WHEN MBDCND = 2 OR 3,
!
! BDYF(I,K) = (D/DY)U(X(I),YF,Z(K)), I=1,2,...,L+1,
! K=1,2,...,N+1.
!
! WHEN MBDCND HAS ANY OTHER VALUE, BDYF IS A DUMMY VARIABLE.
! BDYF MUST BE DIMENSIONED AT LEAST (L+1)*(N+1).
!
! ZS,ZF
! THE RANGE OF Z, I.E. ZS .LE. Z .LE. ZF.
! ZS MUST BE LESS THAN ZF.
!
! N
! THE NUMBER OF PANELS INTO WHICH THE INTERVAL (ZS,ZF) IS
! SUBDIVIDED. HENCE, THERE WILL BE N+1 GRID POINTS IN THE
! Z-DIRECTION GIVEN BY Z(K) = ZS+(K-1)DZ FOR K=1,2,...,N+1,
! WHERE DZ = (ZF-ZS)/N IS THE PANEL WIDTH. N MUST BE AT LEAST 5.
!
! NBDCND
! INDICATES THE TYPE OF BOUNDARY CONDITIONS AT Z = ZS AND Z = ZF.
!
! = 0 IF THE SOLUTION IS PERIODIC IN Z, I.E.
! U(I,J,N+K) = U(I,J,K).
! = 1 IF THE SOLUTION IS SPECIFIED AT Z = ZS AND Z = ZF.
! = 2 IF THE SOLUTION IS SPECIFIED AT Z = ZS AND THE DERIVATIVE
! OF THE SOLUTION WITH RESPECT TO Z IS SPECIFIED AT Z = ZF.
! = 3 IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Z IS
! SPECIFIED AT Z = ZS AND Z = ZF.
! = 4 IF THE DERIVATIVE OF THE SOLUTION WITH RESPECT TO Z IS
! SPECIFIED AT Z = ZS AND THE SOLUTION IS SPECIFIED AT Z=ZF.
!
! BDZS
! A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE
! DERIVATIVE OF THE SOLUTION WITH RESPECT TO Z AT Z = ZS.
! WHEN NBDCND = 3 OR 4,
!
! BDZS(I,J) = (D/DZ)U(X(I),Y(J),ZS), I=1,2,...,L+1,
! J=1,2,...,M+1.
!
! WHEN NBDCND HAS ANY OTHER VALUE, BDZS IS A DUMMY VARIABLE.
! BDZS MUST BE DIMENSIONED AT LEAST (L+1)*(M+1).
!
! BDZF
! A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE
! DERIVATIVE OF THE SOLUTION WITH RESPECT TO Z AT Z = ZF.
! WHEN NBDCND = 2 OR 3,
!
! BDZF(I,J) = (D/DZ)U(X(I),Y(J),ZF), I=1,2,...,L+1,
! J=1,2,...,M+1.
!
! WHEN NBDCND HAS ANY OTHER VALUE, BDZF IS A DUMMY VARIABLE.
! BDZF MUST BE DIMENSIONED AT LEAST (L+1)*(M+1).
!
! lambda
! THE CONSTANT LAMBDA IN THE HELMHOLTZ EQUATION. IF
! LAMBDA .GT. 0, A SOLUTION MAY NOT EXIST. HOWEVER, HW3CRT WILL
! ATTEMPT TO FIND A SOLUTION.
!
! F
! A THREE-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE
! RIGHT SIDE OF THE HELMHOLTZ EQUATION AND BOUNDARY VALUES (IF
! ANY). FOR I=2,3,...,L, J=2,3,...,M, AND K=2,3,...,N
!
! F(I,J,K) = F(X(I),Y(J),Z(K)).
!
! ON THE BOUNDARIES F IS DEFINED BY
!
! LBDCND F(1,J,K) F(L+1,J,K)
! ------ --------------- ---------------
!
! 0 F(XS,Y(J),Z(K)) F(XS,Y(J),Z(K))
! 1 U(XS,Y(J),Z(K)) U(XF,Y(J),Z(K))
! 2 U(XS,Y(J),Z(K)) F(XF,Y(J),Z(K)) J=1,2,...,M+1
! 3 F(XS,Y(J),Z(K)) F(XF,Y(J),Z(K)) K=1,2,...,N+1
! 4 F(XS,Y(J),Z(K)) U(XF,Y(J),Z(K))
!
! MBDCND F(I,1,K) F(I,M+1,K)
! ------ --------------- ---------------
!
! 0 F(X(I),YS,Z(K)) F(X(I),YS,Z(K))
! 1 U(X(I),YS,Z(K)) U(X(I),YF,Z(K))
! 2 U(X(I),YS,Z(K)) F(X(I),YF,Z(K)) I=1,2,...,L+1
! 3 F(X(I),YS,Z(K)) F(X(I),YF,Z(K)) K=1,2,...,N+1
! 4 F(X(I),YS,Z(K)) U(X(I),YF,Z(K))
!
! NBDCND F(I,J,1) F(I,J,N+1)
! ------ --------------- ---------------
!
! 0 F(X(I),Y(J),ZS) F(X(I),Y(J),ZS)
! 1 U(X(I),Y(J),ZS) U(X(I),Y(J),ZF)
! 2 U(X(I),Y(J),ZS) F(X(I),Y(J),ZF) I=1,2,...,L+1
! 3 F(X(I),Y(J),ZS) F(X(I),Y(J),ZF) J=1,2,...,M+1
! 4 F(X(I),Y(J),ZS) U(X(I),Y(J),ZF)
!
! F MUST BE DIMENSIONED AT LEAST (L+1)*(M+1)*(N+1).
!
! NOTE:
!
! IF THE TABLE CALLS FOR BOTH THE SOLUTION U AND THE RIGHT SIDE F
! ON A BOUNDARY, THEN THE SOLUTION MUST BE SPECIFIED.
!
! LDIMF
! THE ROW (OR FIRST) DIMENSION OF THE ARRAYS F,BDYS,BDYF,BDZS,
! AND BDZF AS IT APPEARS IN THE PROGRAM CALLING HW3CRT. THIS
! PARAMETER IS USED TO SPECIFY THE VARIABLE DIMENSION OF THESE
! ARRAYS. LDIMF MUST BE AT LEAST L+1.
!
! MDIMF
! THE COLUMN (OR SECOND) DIMENSION OF THE ARRAY F AND THE ROW (OR
! FIRST) DIMENSION OF THE ARRAYS BDXS AND BDXF AS IT APPEARS IN
! THE PROGRAM CALLING HW3CRT. THIS PARAMETER IS USED TO SPECIFY
! THE VARIABLE DIMENSION OF THESE ARRAYS.
! MDIMF MUST BE AT LEAST M+1.
!
! W
! A ONE-DIMENSIONAL ARRAY THAT MUST BE PROVIDED BY THE USER FOR
! WORK SPACE. THE LENGTH OF W MUST BE AT LEAST 30 + L + M + 5*N
! + MAX(L,M,N) + 7*(INT((L+1)/2) + INT((M+1)/2))
!
!
! * * * * * * ON OUTPUT * * * * * *
!
! F
! CONTAINS THE SOLUTION U(I,J,K) OF THE FINITE DIFFERENCE
! APPROXIMATION FOR THE GRID POINT (X(I),Y(J),Z(K)) FOR
! I=1,2,...,L+1, J=1,2,...,M+1, AND K=1,2,...,N+1.
!
! PERTRB
! IF A COMBINATION OF PERIODIC OR DERIVATIVE BOUNDARY CONDITIONS
! IS SPECIFIED FOR A POISSON EQUATION (LAMBDA = 0), A SOLUTION
! MAY NOT EXIST. PERTRB IS A CONSTANT, CALCULATED AND SUBTRACTED
! FROM F, WHICH ENSURES THAT A SOLUTION EXISTS. PWSCRT THEN
! COMPUTES THIS SOLUTION, WHICH IS A LEAST SQUARES SOLUTION TO
! THE ORIGINAL APPROXIMATION. THIS SOLUTION IS NOT UNIQUE AND IS
! UNNORMALIZED. THE VALUE OF PERTRB SHOULD BE SMALL COMPARED TO
! THE RIGHT SIDE F. OTHERWISE, A SOLUTION IS OBTAINED TO AN
! ESSENTIALLY DIFFERENT PROBLEM. THIS COMPARISON SHOULD ALWAYS
! BE MADE TO INSURE THAT A MEANINGFUL SOLUTION HAS BEEN OBTAINED.
!
! IERROR
! AN ERROR FLAG THAT INDICATES INVALID INPUT PARAMETERS. EXCEPT
! FOR NUMBERS 0 AND 12, A SOLUTION IS NOT ATTEMPTED.
!
! = 0 NO ERROR
! = 1 XS .GE. XF
! = 2 L .LT. 5
! = 3 LBDCND .LT. 0 .OR. LBDCND .GT. 4
! = 4 YS .GE. YF
! = 5 M .LT. 5
! = 6 MBDCND .LT. 0 .OR. MBDCND .GT. 4
! = 7 ZS .GE. ZF
! = 8 N .LT. 5
! = 9 NBDCND .LT. 0 .OR. NBDCND .GT. 4
! = 10 LDIMF .LT. L+1
! = 11 MDIMF .LT. M+1
! = 12 LAMBDA .GT. 0
!
! SINCE THIS IS THE ONLY MEANS OF INDICATING A POSSIBLY INCORRECT
! CALL TO HW3CRT, THE USER SHOULD TEST IERROR AFTER THE CALL.
!
!
! * * * * * * * PROGRAM SPECIFICATIONS * * * * * * * * * * * *
!
! DIMENSION OF BDXS(MDIMF,N+1),BDXF(MDIMF,N+1),BDYS(LDIMF,N+1),
! ARGUMENTS BDYF(LDIMF,N+1),BDZS(LDIMF,M+1),BDZF(LDIMF,M+1),
! F(LDIMF,MDIMF,N+1),W(SEE ARGUMENT LIST)
!
! LATEST DECEMBER 1, 1978
! REVISION
!
! SUBPROGRAMS HW3CRT,POIS3D,POS3D1,TRID,RFFTI,RFFTF,RFFTF1,
! REQUIRED RFFTB,RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF,
! COSQF1,COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI,
! CFFTI1,CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB,
! CFFTF,CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF,
!
! SPECIAL NONE
! CONDITIONS
!
! COMMON VALUE
! BLOCKS
!
! I/O NONE
!
! PRECISION SINGLE
!
! SPECIALIST ROLAND SWEET
!
! LANGUAGE FORTRAN
!
! HISTORY WRITTEN BY ROLAND SWEET AT NCAR IN JULY,1977
!
! ALGORITHM THIS SUBROUTINE DEFINES THE FINITE DIFFERENCE
! EQUATIONS, INCORPORATES BOUNDARY DATA, AND
! ADJUSTS THE RIGHT SIDE OF SINGULAR SYSTEMS AND
! THEN CALLS POIS3D TO SOLVE THE SYSTEM.
!
! SPACE 7862(DECIMAL) = 17300(OCTAL) LOCATIONS ON THE
! REQUIRED NCAR CONTROL DATA 7600
!
! TIMING AND THE EXECUTION TIME T ON THE NCAR CONTROL DATA
! ACCURACY 7600 FOR SUBROUTINE HW3CRT IS ROUGHLY PROPORTIONAL
! TO L*M*N*(LOG2(L)+LOG2(M)+5), BUT ALSO DEPENDS ON
! INPUT PARAMETERS LBDCND AND MBDCND. SOME TYPICAL
! VALUES ARE LISTED IN THE TABLE BELOW.
! THE SOLUTION PROCESS EMPLOYED RESULTS IN A LOSS
! OF NO MORE THAN THREE SIGNIFICANT DIGITS FOR L,M AN
! N AS LARGE AS 32. MORE DETAILED INFORMATION ABOUT
! ACCURACY CAN BE FOUND IN THE DOCUMENTATION FOR
! SUBROUTINE POIS3D WHICH IS THE ROUTINE THAT ACTUALL
! SOLVES THE FINITE DIFFERENCE EQUATIONS.
!
!
! L(=M=N) LBDCND(=MBDCND=NBDCND) T(MSECS)
! ------- ---------------------- --------
!
! 16 0 300
! 16 1 302
! 16 3 348
! 32 0 1925
! 32 1 1929
! 32 3 2109
!
! PORTABILITY AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN.
! THE MACHINE DEPENDENT CONSTANT PI IS DEFINED IN
! FUNCTION PIMACH.
!
! REQUIRED COS,SIN,ATAN
! RESIDENT
! ROUTINES
!
! REFERENCE NONE
!
! REQUIRED COS,SIN,ATAN
! RESIDENT
! ROUTINES
!
! REFERENCE NONE
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
    USE params, ONLY: DP, PI
    IMPLICIT REAL(DP) (A - H, O - Z)

    REAL(DP), INTENT(IN) :: XS ! min value of x
    REAL(DP), INTENT(IN) :: XF ! max value of x
    INTEGER, INTENT(IN) :: L ! There are L+1 grid points IN x direction
    INTEGER, INTENT(IN) :: LBDCND ! Boundary condition on xs, xf
    REAL(DP), INTENT(IN) :: BDXS(MDIMF, *) ! values of dF/dx at x=xs boundary
    REAL(DP), INTENT(IN) :: BDXF(MDIMF, *) ! values of dF/dx x=xf boundary

    REAL(DP), INTENT(IN) :: YS ! min value of y
    REAL(DP), INTENT(IN) :: YF ! max value of y
    INTEGER, INTENT(IN) :: M ! there are M+1 grid points IN y direction
    INTEGER, INTENT(IN) :: MBDCND ! boundary condition on ys, yf
    REAL(DP), INTENT(IN):: BDYS(LDIMF, *) ! values of dF/dy at y=ys boundary
    REAL(DP), INTENT(IN):: BDYF(LDIMF, *) ! values of dF/dY at y=yf boundary

    REAL(DP), INTENT(IN) :: ZS ! min value of z
    REAL(DP), INTENT(IN) :: ZF ! max value of z
    INTEGER, INTENT(IN) :: N ! there are N+1 grid points IN z direction
    INTEGER, INTENT(IN) :: NBDCND ! boundary condition on zs, zf
    REAL(DP), INTENT(IN):: BDZS(LDIMF, *) ! values of dF/dy at z=zs boundary
    REAL(DP), INTENT(IN):: BDZF(LDIMF, *) ! values of dF/dY at z=zf boundary

    REAL(DP), INTENT(IN) :: lambda ! lambda IN the Helmoltz equation
    INTEGER, INTENT(IN) :: LDIMF ! row (or first) DIMENSION of the arrays F, BDYS, BDYF, BDZS, and BDZF
    INTEGER, INTENT(IN) :: MDIMF ! column (or second) DIMENSION of the array ARRAY F and
    ! the row (or first) DIMENSION of the arrays BDXS and BDXF
    REAL(DP), INTENT(INOUT):: F(LDIMF, MDIMF, *) ! On input : right side of the equation and boundary values, IF any.
    ! On output : computed solution of the equation
    REAL(DP), INTENT(OUT) :: PERTRB ! indicates IF the solution matches the problem.
    INTEGER, INTENT(OUT) :: IERROR ! error flag
    REAL(DP), INTENT(INOUT):: W(:)
    ! 30 + L + M + 5*N + MAX(L,M,N) + 7*(INT((L+1)/2) + INT((M+1)/2))
! CHECK FOR INVALID INPUT.
!
    IERROR = 0
    IF (XF .LE. XS) IERROR = 1
    IF (L .LT. 5) IERROR = 2
    IF (LBDCND .LT. 0 .OR. LBDCND .GT. 4) IERROR = 3
    IF (YF .LE. YS) IERROR = 4
    IF (M .LT. 5) IERROR = 5
    IF (MBDCND .LT. 0 .OR. MBDCND .GT. 4) IERROR = 6
    IF (ZF .LE. ZS) IERROR = 7
    IF (N .LT. 5) IERROR = 8
    IF (NBDCND .LT. 0 .OR. NBDCND .GT. 4) IERROR = 9
    IF (LDIMF .LT. L + 1) IERROR = 10
    IF (MDIMF .LT. M + 1) IERROR = 11

    IF (IERROR .NE. 0) RETURN
!
! Prepare y DIMENSION bounds
!
    DY = (YF - YS)/REAL(M, DP)
    TWBYDY = 2D0/DY
    C2 = 1D0/(DY**2)
    MP1 = M + 1
    MP = MBDCND + 1

    SELECT CASE (MBDCND)
    CASE (0) ! periodic BC
      MSTART = 1
      MSTOP = M
    CASE (1) ! F fixed on ys, yf
      MSTART = 2
      MSTOP = M
    CASE (2) ! F fixed on ys, dF/dy fixed on yf
      MSTART = 2
      MSTOP = MP1
    CASE (3) ! dF/dy fixed on ys, yf
      MSTART = 1
      MSTOP = MP1
    CASE (4) ! dF/dy fixed on ys, F fixed on yf
      MSTART = 1
      MSTOP = M
    END SELECT

    MUNK = MSTOP - MSTART + 1

! Prepare z DIMENSION bounds

    DZ = (ZF - ZS)/FLOAT(N)
    TWBYDZ = 2D0/DZ
    C3 = 1D0/(DZ**2)
    NP1 = N + 1
    NP = NBDCND + 1

    SELECT CASE (NBDCND)
    CASE (0) ! periodic BC
      NSTART = 1
      NSTOP = N
    CASE (1) ! F fixed on zs, zf
      NSTART = 2
      NSTOP = N
    CASE (2) ! F fixed on zs, dF/dz fixed on zf
      NSTART = 2
      NSTOP = NP1
    CASE (3) ! dF/dz fixed on zs, zf
      NSTART = 1
      NSTOP = NP1
    CASE (4) ! dF/dz fixed on zs, F fixed on zf
      NSTART = 1
      NSTOP = N
    END SELECT

    NUNK = NSTOP - NSTART + 1

!
! Prepare x DIMENSION
!
    LP1 = L + 1
    DX = (XF - XS)/FLOAT(L)
    C1 = 1D0/(DX**2)
    TWBYDX = 2D0/DX
    LP = LBDCND + 1
    LSTART = 1
    LSTOP = L
!
! Enter boundary DATA for x-boundaries.
!
    SELECT CASE (LBDCND)
    CASE (0) ! periodic BC
      LSTART = 1
      LSTOP = L
    CASE (1) ! F fixed on xs, xf
      LSTART = 2
      LSTOP = L
      DO J = MSTART, MSTOP
        DO K = NSTART, NSTOP
          F(2, J, K) = F(2, J, K) - C1*F(1, J, K)
          F(L, J, K) = F(L, J, K) - C1*F(LP1, J, K)
        END DO
      END DO
    CASE (2) ! F fixed on xs, dF/dx fixed on xf
      LSTART = 2
      LSTOP = LP1
      DO J = MSTART, MSTOP
        DO K = NSTART, NSTOP
          F(2, J, K) = F(2, J, K) - C1*F(1, J, K)
          F(LP1, J, K) = F(LP1, J, K) - TWBYDX*BDXF(J, K)
        END DO
      END DO
    CASE (3)! dF/dx fixed on xs, xf
      LSTART = 1
      LSTOP = LP1
      DO J = MSTART, MSTOP
        DO K = NSTART, NSTOP
          F(1, J, K) = F(1, J, K) + TWBYDX*BDXS(J, K)
          F(LP1, J, K) = F(LP1, J, K) - TWBYDX*BDXF(J, K)
        END DO
      END DO
    CASE (4)! dF/dx fixed on xs, F fixed on xf
      LSTART = 1
      LSTOP = L
      DO J = MSTART, MSTOP
        DO K = NSTART, NSTOP
          F(1, J, K) = F(1, J, K) + TWBYDX*BDXS(J, K)
          F(L, J, K) = F(L, J, K) - C1*F(LP1, J, K)
        END DO
      END DO
    END SELECT

    LUNK = LSTOP - LSTART + 1
!
! Enter boundary DATA for y-boundaries.
!
    SELECT CASE (MBDCND)
      !CASE(0) : nothing TO be done.
    CASE (1) ! F fixed on ys, yf
      DO I = LSTART, LSTOP
        DO K = NSTART, NSTOP
          F(I, 2, K) = F(I, 2, K) - C2*F(I, 1, K)
          F(I, M, K) = F(I, M, K) - C2*F(I, MP1, K)
        END DO
      END DO
    CASE (2) ! F fixed on ys, dF/dy fixed on yf
      DO I = LSTART, LSTOP
        DO K = NSTART, NSTOP
          F(I, 2, K) = F(I, 2, K) - C2*F(I, 1, K)
          F(I, MP1, K) = F(I, MP1, K) - TWBYDY*BDYF(I, K)
        END DO
      END DO
    CASE (3) ! dF/dy fixed on ys, yf
      DO I = LSTART, LSTOP
        DO K = NSTART, NSTOP
          F(I, 1, K) = F(I, 1, K) + TWBYDY*BDYS(I, K)
          F(I, MP1, K) = F(I, MP1, K) - TWBYDY*BDYF(I, K)
        END DO
      END DO
    CASE (4) ! dF/dy fixed on ys, F fixed on yf
      DO I = LSTART, LSTOP
        DO K = NSTART, NSTOP
          F(I, 1, K) = F(I, 1, K) + TWBYDY*BDYS(I, K)
          F(I, M, K) = F(I, M, K) - C2*F(I, MP1, K)
        END DO
      END DO
    END SELECT
!
! Enter boundary DATA for z-boundaries.
!
    SELECT CASE (NBDCND)
      !CASE(0) nothing TO be done
    CASE (1) ! F fixed on zs, zf
      DO I = LSTART, LSTOP
        DO J = MSTART, MSTOP
          F(I, J, 2) = F(I, J, 2) - C3*F(I, J, 1)
          F(I, J, N) = F(I, J, N) - C3*F(I, J, NP1)
        END DO
      END DO
    CASE (2) ! F fixed on zs, dF/dz fixed on zf
      DO I = LSTART, LSTOP
        DO J = MSTART, MSTOP
          F(I, J, 2) = F(I, J, 2) - C3*F(I, J, 1)
          F(I, J, NP1) = F(I, J, NP1) - TWBYDZ*BDZF(I, J)
        END DO
      END DO
    CASE (3) ! dF/dz fixed on zs, zf
      DO I = LSTART, LSTOP
        DO J = MSTART, MSTOP
          F(I, J, 1) = F(I, J, 1) + TWBYDZ*BDZS(I, J)
          F(I, J, NP1) = F(I, J, NP1) - TWBYDZ*BDZF(I, J)
        END DO
      END DO
    CASE (4) ! dF/dz fixed on zs, F fixed on zf
      DO I = LSTART, LSTOP
        DO J = MSTART, MSTOP
          F(I, J, 1) = F(I, J, 1) + TWBYDZ*BDZS(I, J)
          F(I, J, N) = F(I, J, N) - C3*F(I, J, NP1)
        END DO
      END DO
    END SELECT
!
! Define A,B,C coefficients IN w-array.
!
    IWB = NUNK + 1
    IWC = IWB + NUNK
    IWW = IWC + NUNK

    DO K = 1, NUNK
      I = IWC + K - 1
      W(K) = C3
      W(I) = C3
      I = IWB + K - 1
      W(I) = -2.0D0*C3 + lambda
    END DO

    SELECT CASE (NBDCND)
      !CASE(0,1) : periodic, or F fixed at both ends : nothing TO be done
    CASE (2) ! F fixed on zs, dF/dz fixed on zf
      W(IWB - 1) = 2.0D0*C3
    CASE (3) ! dF/dz fixed on zs, zf
      W(IWC) = 2.0D0*C3
      W(IWB - 1) = 2.0D0*C3
    CASE (4) ! dF/dz fixed on zs, F fixed on zf
      W(IWC) = 2.0D0*C3
    END SELECT

    PERTRB = 0.0D0

!
! For singular problems adjust DATA TO insure a solution will EXIST.
! singular problems means F is fixed on no boundary (neither LBDCND, MBDCND nor NBDCND are equal TO 1,2 or 4)
!
!~ IF( (LBDCND \= 1) .AND. (LBDCND \= 2) .AND. (LBDCND \= 4) .AND. &
!~ &(MBDCND \= 1) .AND. (MBDCND \= 2) .AND. (MBDCND \= 4) .AND. &
!~ &(NBDCND \= 1) .AND. (NBDCND \= 2) .AND. (NBDCND \= 4) ) THEN

    GO TO(156, 172, 172, 156, 172), LP
156 GO TO(157, 172, 172, 157, 172), MP
157 GO TO(158, 172, 172, 158, 172), NP
158 IF (lambda) 172, 160, 159
159 IERROR = 12
    GO TO 172
160 CONTINUE
    MSTPM1 = MSTOP - 1
    LSTPM1 = LSTOP - 1
    NSTPM1 = NSTOP - 1
    XLP = (2 + LP)/3
    YLP = (2 + MP)/3
    ZLP = (2 + NP)/3
    S1 = 0.0D0
    DO 164 K = 2, NSTPM1
      DO 162 J = 2, MSTPM1
        DO 161 I = 2, LSTPM1
          S1 = S1 + F(I, J, K)
161       CONTINUE
          S1 = S1 + (F(1, J, K) + F(LSTOP, J, K))/XLP
162       CONTINUE
          S2 = 0.0D0
          DO 163 I = 2, LSTPM1
            S2 = S2 + F(I, 1, K) + F(I, MSTOP, K)
163         CONTINUE
            S2 = (S2 + (F(1, 1, K) + F(1, MSTOP, K) + F(LSTOP, 1, K) + F(LSTOP, MSTOP, K))/&
            & XLP)/YLP
            S1 = S1 + S2
164         CONTINUE
            S = (F(1, 1, 1) + F(LSTOP, 1, 1) + F(1, 1, NSTOP) + F(LSTOP, 1, NSTOP) +&
            & F(1, MSTOP, 1) + F(LSTOP, MSTOP, 1) + F(1, MSTOP, NSTOP) +&
            & F(LSTOP, MSTOP, NSTOP))/(XLP*YLP)
            DO 166 J = 2, MSTPM1
              DO 165 I = 2, LSTPM1
                S = S + F(I, J, 1) + F(I, J, NSTOP)
165             CONTINUE
166             CONTINUE
                S2 = 0.0D0
                DO 167 I = 2, LSTPM1
                  S2 = S2 + F(I, 1, 1) + F(I, 1, NSTOP) + F(I, MSTOP, 1) + F(I, MSTOP, NSTOP)
167               CONTINUE
                  S = S2/YLP + S
                  S2 = 0.0D0
                  DO 168 J = 2, MSTPM1
                    S2 = S2 + F(1, J, 1) + F(1, J, NSTOP) + F(LSTOP, J, 1) + F(LSTOP, J, NSTOP)
168                 CONTINUE
                    S = S2/XLP + S
                    PERTRB = (S/ZLP + S1)/((FLOAT(LUNK + 1) - XLP)*(FLOAT(MUNK + 1) - YLP)*&
                    & (FLOAT(NUNK + 1) - ZLP))
                    DO 171 I = 1, LUNK
                      DO 170 J = 1, MUNK
                        DO 169 K = 1, NUNK
                          F(I, J, K) = F(I, J, K) - PERTRB
169                       CONTINUE
170                       CONTINUE
171                       CONTINUE
172                       CONTINUE
                          !END IF
                          NPEROD = 0
                          IF (NBDCND .EQ. 0) GO TO 173
                          NPEROD = 1
                          W(1) = 0.0D0
                          W(IWW - 1) = 0.0D0
173                       CONTINUE
                          CALL POIS3D(LBDCND, LUNK, C1, MBDCND, MUNK, C2, NPEROD, NUNK, W, W(IWB:),&
                          & W(IWC:), LDIMF, MDIMF, F(LSTART, MSTART, NSTART), IR, W(IWW:))
!
! FILL IN SIDES FOR PERIODIC BOUNDARY CONDITIONS.
!
                          IF (LP .NE. 1) GO TO 180
                          IF (MP .NE. 1) GO TO 175
                          DO 174 K = NSTART, NSTOP
                            F(1, MP1, K) = F(1, 1, K)
174                         CONTINUE
                            MSTOP = MP1
175                         IF (NP .NE. 1) GO TO 177
                            DO 176 J = MSTART, MSTOP
                              F(1, J, NP1) = F(1, J, 1)
176                           CONTINUE
                              NSTOP = NP1
177                           DO 179 J = MSTART, MSTOP
                                DO 178 K = NSTART, NSTOP
                                  F(LP1, J, K) = F(1, J, K)
178                               CONTINUE
179                               CONTINUE
180                               CONTINUE
                                  IF (MP .NE. 1) GO TO 185
                                  IF (NP .NE. 1) GO TO 182
                                  DO 181 I = LSTART, LSTOP
                                    F(I, 1, NP1) = F(I, 1, 1)
181                                 CONTINUE
                                    NSTOP = NP1
182                                 DO 184 I = LSTART, LSTOP
                                      DO 183 K = NSTART, NSTOP
                                        F(I, MP1, K) = F(I, 1, K)
183                                     CONTINUE
184                                     CONTINUE
185                                     CONTINUE
                                        IF (NP .NE. 1) GO TO 188
                                        DO 187 I = LSTART, LSTOP
                                          DO 186 J = MSTART, MSTOP
                                            F(I, J, NP1) = F(I, J, 1)
186                                         CONTINUE
187                                         CONTINUE
188                                         CONTINUE
                                            RETURN
                                            END SUBROUTINE
!~ FUNCTION PIMACH (DUM)
!~ !
!~ ! THIS SUBPROGRAM SUPPLIES THE VALUE OF THE CONSTANT PI CORRECT TO
!~ ! MACHINE PRECISION WHERE
!~ !
!~ ! PI=3.1415926535897932384626433832795028841971693993751058209749446
!~ !
!~ PIMACH = 3.14159265358979D0
!~ RETURN
!~ END
                                            SUBROUTINE POIS3D(LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, C, LDIMF,&
                                            & MDIMF, F, IERROR, W)
!
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * *
! * F I S H P A K *
! * *
! * *
! * A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF *
! * *
! * SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS *
! * *
! * (VERSION 3.1 , OCTOBER 1980) *
! * *
! * BY *
! * *
! * JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET *
! * *
! * OF *
! * *
! * THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH *
! * *
! * BOULDER, COLORADO (80307) U.S.A. *
! * *
! * WHICH IS SPONSORED BY *
! * *
! * THE NATIONAL SCIENCE FOUNDATION *
! * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! * * * * * * * * * PURPOSE * * * * * * * * * * * * * * * * * *
!
! SUBROUTINE POIS3D SOLVES THE LINEAR SYSTEM OF EQUATIONS
!
! C1*(X(I-1,J,K)-2.*X(I,J,K)+X(I+1,J,K))
! + C2*(X(I,J-1,K)-2.*X(I,J,K)+X(I,J+1,K))
! + A(K)*X(I,J,K-1)+B(K)*X(I,J,K)+C(K)*X(I,J,K+1) = F(I,J,K)
!
! FOR I=1,2,...,L , J=1,2,...,M , AND K=1,2,...,N .
!
! THE INDICES K-1 AND K+1 ARE EVALUATED MODULO N, I.E.
! X(I,J,0) = X(I,J,N) AND X(I,J,N+1) = X(I,J,1). THE UNKNOWNS
! X(0,J,K), X(L+1,J,K), X(I,0,K), AND X(I,M+1,K) ARE ASSUMED TO TAKE
! ON CERTAIN PRESCRIBED VALUES DESCRIBED BELOW.
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! * * * * * * * * PARAMETER DESCRIPTION * * * * * * * * * *
!
!
! * * * * * * ON INPUT * * * * * *
!
! LPEROD INDICATES THE VALUES THAT X(0,J,K) AND X(L+1,J,K) ARE
! ASSUMED TO HAVE.
!
! = 0 IF X(0,J,K) = X(L,J,K) AND X(L+1,J,K) = X(1,J,K).
! = 1 IF X(0,J,K) = X(L+1,J,K) = 0.
! = 2 IF X(0,J,K) = 0 AND X(L+1,J,K) = X(L-1,J,K).
! = 3 IF X(0,J,K) = X(2,J,K) AND X(L+1,J,K) = X(L-1,J,K).
! = 4 IF X(0,J,K) = X(2,J,K) AND X(L+1,J,K) = 0.
!
! L THE NUMBER OF UNKNOWNS IN THE I-DIRECTION. L MUST BE AT
! LEAST 3.
!
! C1 THE REAL CONSTANT THAT APPEARS IN THE ABOVE EQUATION.
!
! MPEROD INDICATES THE VALUES THAT X(I,0,K) AND X(I,M+1,K) ARE
! ASSUMED TO HAVE.
!
! = 0 IF X(I,0,K) = X(I,M,K) AND X(I,M+1,K) = X(I,1,K).
! = 1 IF X(I,0,K) = X(I,M+1,K) = 0.
! = 2 IF X(I,0,K) = 0 AND X(I,M+1,K) = X(I,M-1,K).
! = 3 IF X(I,0,K) = X(I,2,K) AND X(I,M+1,K) = X(I,M-1,K).
! = 4 IF X(I,0,K) = X(I,2,K) AND X(I,M+1,K) = 0.
!
! M THE NUMBER OF UNKNOWNS IN THE J-DIRECTION. M MUST BE AT
! LEAST 3.
!
! C2 THE REAL CONSTANT WHICH APPEARS IN THE ABOVE EQUATION.
!
! NPEROD = 0 IF A(1) AND C(N) ARE NOT ZERO.
! = 1 IF A(1) = C(N) = 0.
!
! N THE NUMBER OF UNKNOWNS IN THE K-DIRECTION. N MUST BE AT
! LEAST 3.
!
!
! A,B,C ONE-DIMENSIONAL ARRAYS OF LENGTH N THAT SPECIFY THE
! COEFFICIENTS IN THE LINEAR EQUATIONS GIVEN ABOVE.
!
! IF NPEROD = 0 THE ARRAY ELEMENTS MUST NOT DEPEND UPON THE
! INDEX K, BUT MUST BE CONSTANT. SPECIFICALLY,THE
! SUBROUTINE CHECKS THE FOLLOWING CONDITION
!
! A(K) = C(1)
! C(K) = C(1)
! B(K) = B(1)
!
! FOR K=1,2,...,N.
!
! LDIMF THE ROW (OR FIRST) DIMENSION OF THE THREE-DIMENSIONAL
! ARRAY F AS IT APPEARS IN THE PROGRAM CALLING POIS3D.
! THIS PARAMETER IS USED TO SPECIFY THE VARIABLE DIMENSION
! OF F. LDIMF MUST BE AT LEAST L.
!
! MDIMF THE COLUMN (OR SECOND) DIMENSION OF THE THREE-DIMENSIONAL
! ARRAY F AS IT APPEARS IN THE PROGRAM CALLING POIS3D.
! THIS PARAMETER IS USED TO SPECIFY THE VARIABLE DIMENSION
! OF F. MDIMF MUST BE AT LEAST M.
!
! F A THREE-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF
! THE RIGHT SIDE OF THE LINEAR SYSTEM OF EQUATIONS GIVEN
! ABOVE. F MUST BE DIMENSIONED AT LEAST L X M X N.
!
! W A ONE-DIMENSIONAL ARRAY THAT MUST BE PROVIDED BY THE
! USER FOR WORK SPACE. THE LENGTH OF W MUST BE AT LEAST
! 30 + L + M + 2*N + MAX(L,M,N) +
! 7*(INT((L+1)/2) + INT((M+1)/2)).
!
!
! * * * * * * ON OUTPUT * * * * * *
!
! F CONTAINS THE SOLUTION X.
!
! IERROR AN ERROR FLAG THAT INDICATES INVALID INPUT PARAMETERS.
! EXCEPT FOR NUMBER ZERO, A SOLUTION IS NOT ATTEMPTED.
! = 0 NO ERROR
! = 1 IF LPEROD .LT. 0 OR .GT. 4
! = 2 IF L .LT. 3
! = 3 IF MPEROD .LT. 0 OR .GT. 4
! = 4 IF M .LT. 3
! = 5 IF NPEROD .LT. 0 OR .GT. 1
! = 6 IF N .LT. 3
! = 7 IF LDIMF .LT. L
! = 8 IF MDIMF .LT. M
! = 9 IF A(K) .NE. C(1) OR C(K) .NE. C(1) OR B(I) .NE.B(1)
! FOR SOME K=1,2,...,N.
! = 10 IF NPEROD = 1 AND A(1) .NE. 0 OR C(N) .NE. 0
!
! SINCE THIS IS THE ONLY MEANS OF INDICATING A POSSIBLY
! INCORRECT CALL TO POIS3D, THE USER SHOULD TEST IERROR
! AFTER THE CALL.
!
!
! * * * * * * * PROGRAM SPECIFICATIONS * * * * * * * * * * * *
!
! DIMENSION OF A(N),B(N),C(N),F(LDIMF,MDIMF,N),
! ARGUMENTS W(SEE ARGUMENT LIST)
!
! LATEST DECEMBER 1, 1978
! REVISION
!
! SUBPROGRAMS POIS3D,POS3D1,TRID,RFFTI,RFFTF,RFFTF1,RFFTB,
! REQUIRED RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF,COSQF1
! COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI,CFFTI1,
! CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB,CFFTF,
! CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF,PIMACH,
!
! SPECIAL NONE
! CONDITIONS
!
! COMMON VALUE
! BLOCKS
!
! I/O NONE
!
! PRECISION SINGLE
!
! SPECIALIST ROLAND SWEET
!
! LANGUAGE FORTRAN
!
! HISTORY WRITTEN BY ROLAND SWEET AT NCAR IN JULY,1977
!
! ALGORITHM THIS SUBROUTINE SOLVES THREE-DIMENSIONAL BLOCK
! TRIDIAGONAL LINEAR SYSTEMS ARISING FROM FINITE
! DIFFERENCE APPROXIMATIONS TO THREE-DIMENSIONAL
! POISSON EQUATIONS USING THE FOURIER TRANSFORM
! PACKAGE SCLRFFTPAK WRITTEN BY PAUL SWARZTRAUBER.
!
! SPACE 6561(DECIMAL) = 14641(OCTAL) LOCATIONS ON THE
! REQUIRED NCAR CONTROL DATA 7600
!
! TIMING AND THE EXECUTION TIME T ON THE NCAR CONTROL DATA
! ACCURACY 7600 FOR SUBROUTINE POIS3D IS ROUGHLY PROPORTIONAL
! TO L*M*N*(LOG2(L)+LOG2(M)+5), BUT ALSO DEPENDS ON
! INPUT PARAMETERS LPEROD AND MPEROD. SOME TYPICAL
! VALUES ARE LISTED IN THE TABLE BELOW WHEN NPEROD=0.
! TO MEASURE THE ACCURACY OF THE ALGORITHM A
! UNIFORM RANDOM NUMBER GENERATOR WAS USED TO CREATE
! A SOLUTION ARRAY X FOR THE SYSTEM GIVEN IN THE
! 'PURPOSE' WITH
!
! A(K) = C(K) = -0.5*B(K) = 1, K=1,2,...,N
!
! AND, WHEN NPEROD = 1
!
! A(1) = C(N) = 0
! A(N) = C(1) = 2.
!
! THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN SYS-
! TEM AND, USING DOUBLE PRECISION, A RIGHT SIDE Y WAS
! COMPUTED. USING THIS ARRAY Y SUBROUTINE POIS WAS
! CALLED TO PRODUCE AN APPROXIMATE SOLUTION Z. THEN
! THE RELATIVE ERROR, DEFINED AS
!
! E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K)))
!
! WHERE THE TWO MAXIMA ARE TAKEN OVER I=1,2,...,L,
! J=1,2,...,M AND K=1,2,...,N, WAS COMPUTED. THE
! VALUE OF E IS GIVEN IN THE TABLE BELOW FOR SOME
! TYPICAL VALUES OF L,M AND N.
!
!
! L(=M=N) LPEROD MPEROD T(MSECS) E
! ------ ------ ------ -------- ------
!
! 16 0 0 272 1.E-13
! 15 1 1 287 4.E-13
! 17 3 3 338 2.E-13
! 32 0 0 1755 2.E-13
! 31 1 1 1894 2.E-12
! 33 3 3 2042 7.E-13
!
!
! PORTABILITY AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN.
! THE MACHINE DEPENDENT CONSTANT PI IS DEFINED IN
! FUNCTION PIMACH.
!
! REQUIRED COS,SIN,ATAN
! RESIDENT
! ROUTINES
!
! REFERENCE NONE
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
                                              USE params, ONLY: DP, PI
                                              IMPLICIT REAL(DP) (A - H, O - Z)

                                              INTEGER :: LPEROD
                                              INTEGER :: L
                                              REAL(DP):: C1
                                              INTEGER :: MPEROD
                                              INTEGER :: M
                                              REAL(DP):: C2
                                              INTEGER :: NPEROD
                                              INTEGER :: N
                                              REAL(DP):: A(:)
                                              REAL(DP):: B(:)
                                              REAL(DP):: C(:)
                                              INTEGER :: LDIMF
                                              INTEGER :: MDIMF
                                              REAL(DP):: F(LDIMF, MDIMF, *)
                                              INTEGER IERROR
                                              REAL(DP) :: W(:)

                                              REAL(DP) :: SAVE(6)
!~ DIMENSION A(1), B(1), C(1), F(LDIMF,MDIMF,1), W(1), SAVE(6)

                                              LP = LPEROD + 1
                                              MP = MPEROD + 1
                                              NP = NPEROD + 1
!
! CHECK FOR INVALID INPUT.
!
                                              IERROR = 0
                                              IF (LP .LT. 1 .OR. LP .GT. 5) IERROR = 1
                                              IF (L .LT. 3) IERROR = 2
                                              IF (MP .LT. 1 .OR. MP .GT. 5) IERROR = 3
                                              IF (M .LT. 3) IERROR = 4
                                              IF (NP .LT. 1 .OR. NP .GT. 2) IERROR = 5
                                              IF (N .LT. 3) IERROR = 6
                                              IF (LDIMF .LT. L) IERROR = 7
                                              IF (MDIMF .LT. M) IERROR = 8
                                              IF (NP .NE. 1) GO TO 103
                                              DO 101 K = 1, N
                                                IF (A(K) .NE. C(1)) GO TO 102
                                                IF (C(K) .NE. C(1)) GO TO 102
                                                IF (B(K) .NE. B(1)) GO TO 102
101                                             CONTINUE
                                                GO TO 104
102                                             IERROR = 9
103                                             IF (NPEROD .EQ. 1 .AND. (A(1) .NE. 0. .OR. C(N) .NE. 0.)) IERROR = 10
104                                             IF (IERROR .NE. 0) GO TO 122
                                                IWYRT = L + 1
                                                IWT = IWYRT + M
                                                IWD = IWT + MAX0(L, M, N) + 1
                                                IWBB = IWD + N
                                                IWX = IWBB + N
                                                IWY = IWX + 7*((L + 1)/2) + 15
                                                GO TO(105, 114), NP
!
! REORDER UNKNOWNS WHEN NPEROD = 0.
!
105                                             NH = (N + 1)/2
                                                NHM1 = NH - 1
                                                NODD = 1
                                                IF (2*NH .EQ. N) NODD = 2
                                                DO 111 I = 1, L
                                                  DO 110 J = 1, M
                                                    DO 106 K = 1, NHM1
                                                      NHPK = NH + K
                                                      NHMK = NH - K
                                                      W(K) = F(I, J, NHMK) - F(I, J, NHPK)
                                                      W(NHPK) = F(I, J, NHMK) + F(I, J, NHPK)
106                                                   CONTINUE
                                                      W(NH) = 2.0D0*F(I, J, NH)
                                                      GO TO(108, 107), NODD
107                                                   W(N) = 2.0D0*F(I, J, N)
108                                                   DO 109 K = 1, N
                                                        F(I, J, K) = W(K)
109                                                     CONTINUE
110                                                     CONTINUE
111                                                     CONTINUE
                                                        SAVE(1) = C(NHM1)
                                                        SAVE(2) = A(NH)
                                                        SAVE(3) = C(NH)
                                                        SAVE(4) = B(NHM1)
                                                        SAVE(5) = B(N)
                                                        SAVE(6) = A(N)
                                                        C(NHM1) = 0.0D0
                                                        A(NH) = 0.0D0
                                                        C(NH) = 2.0D0*C(NH)
                                                        GO TO(112, 113), NODD
112                                                     B(NHM1) = B(NHM1) - A(NH - 1)
                                                        B(N) = B(N) + A(N)
                                                        GO TO 114
113                                                     A(N) = C(NH)
114                                                     CONTINUE
                                                     CALL POS3D1(LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, W, W(IWYRT:), W(IWT:), &
                                                                    W(IWD:), W(IWX:), W(IWY:), C1, C2, W(IWBB:))
                                                        GO TO(115, 122), NP
115                                                     DO 121 I = 1, L
                                                          DO 120 J = 1, M
                                                            DO 116 K = 1, NHM1
                                                              NHMK = NH - K
                                                              NHPK = NH + K
                                                              W(NHMK) = 0.5D0*(F(I, J, NHPK) + F(I, J, K))
                                                              W(NHPK) = 0.5D0*(F(I, J, NHPK) - F(I, J, K))
116                                                           CONTINUE
                                                              W(NH) = 0.5D0*F(I, J, NH)
                                                              GO TO(118, 117), NODD
117                                                           W(N) = 0.5D0*F(I, J, N)
118                                                           DO 119 K = 1, N
                                                                F(I, J, K) = W(K)
119                                                             CONTINUE
120                                                             CONTINUE
121                                                             CONTINUE
                                                                C(NHM1) = SAVE(1)
                                                                A(NH) = SAVE(2)
                                                                C(NH) = SAVE(3)
                                                                B(NHM1) = SAVE(4)
                                                                B(N) = SAVE(5)
                                                                A(N) = SAVE(6)
122                                                             CONTINUE
                                                                RETURN
                                                                END SUBROUTINE

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                                                      SUBROUTINE POS3D1(LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, XRT, YRT, T, D, &
                                                                                  WX, WY, C1, C2, BB)
                                                                  USE params, ONLY: DP, PI
                                                                  IMPLICIT REAL(DP) (A - H, O - Z)

                                                                  INTEGER, INTENT(IN)::LP
                                                                  INTEGER, INTENT(IN)::L
                                                                  INTEGER, INTENT(IN)::MP
                                                                  INTEGER, INTENT(IN)::M
                                                                  INTEGER, INTENT(IN)::N
                                                                  REAL(DP), INTENT(IN)::A(:)
                                                                  REAL(DP), INTENT(IN)::B(:)
                                                                  REAL(DP), INTENT(IN)::C(:)
                                                                  REAL(DP), INTENT(INOUT)::F(LDIMF, MDIMF, *)
                                                                  REAL(DP), INTENT(OUT)::XRT(:)
                                                                  REAL(DP), INTENT(OUT)::YRT(:)
                                                                  REAL(DP), INTENT(INOUT)::T(:) ! SIZE of T is MAX(L,M,N)
                                                                  REAL(DP), INTENT(OUT)::D(:)
                                                                  REAL(DP), INTENT(OUT)::WX(:)
                                                                  REAL(DP), INTENT(OUT)::WY(:)
                                                                  REAL(DP), INTENT(OUT)::BB(:)

!~ REAL(DP),INTENT(OUT),XRT(:)
!~ DIMENSION A(1) ,B(1) ,C(1) ,&
!~ & F(LDIMF,MDIMF,1) ,XRT(1) ,YRT(1) ,&
!~ & T(1) ,D(1) ,WX(1) ,WY(1) ,&
!~ & BB(1)

                                                                  INTEGER :: ifac(10) ! work array for fftpack

!
! GENERATE TRANSFORM ROOTS
!
                                                                  LRDEL = ((LP - 1)*(LP - 3)*(LP - 5))/3
                                                                  SCALX = L + LRDEL
                                                                  DX = PI/(2.0D0*SCALX)
                                                                  GO TO(108, 103, 101, 102, 101), LP
101                                                               DI = 0.5D0
                                                                  SCALX = 2.0D0*SCALX
                                                                  GO TO 104
102                                                               DI = 1.0D0
                                                                  GO TO 104
103                                                               DI = 0.0D0
104                                                               DO 105 I = 1, L
                                                                    XRT(I) = -4.0D0*C1*(SIN((FLOAT(I) - DI)*DX))**2
105                                                                 CONTINUE
                                                                    SCALX = 2.0D0*SCALX
                                                                    GO TO(112, 106, 110, 107, 111), LP
106                                                                 CALL SINTI(L, WX, ifac)
                                                                    GO TO 112
107                                                                 CALL COSTI(L, WX, ifac)
                                                                    GO TO 112
108                                                                 XRT(1) = 0.0D0
                                                                    XRT(L) = -4.0D0*C1
                                                                    DO 109 I = 3, L, 2
                                                                      XRT(I - 1) = -4.0D0*C1*(SIN(FLOAT((I - 1))*DX))**2
                                                                      XRT(I) = XRT(I - 1)
109                                                                   CONTINUE
                                                                      CALL RFFTI(L, WX, ifac)
                                                                      GO TO 112
110                                                                   CALL SINQI(L, WX, ifac)
                                                                      GO TO 112
111                                                                   CALL COSQI(L, WX, ifac)
112                                                                   CONTINUE
                                                                      MRDEL = ((MP - 1)*(MP - 3)*(MP - 5))/3
                                                                      SCALY = M + MRDEL
                                                                      DY = PI/(2.*SCALY)
                                                                      GO TO(120, 115, 113, 114, 113), MP
113                                                                   DJ = 0.5D0
                                                                      SCALY = 2.0D0*SCALY
                                                                      GO TO 116
114                                                                   DJ = 1.0D0
                                                                      GO TO 116
115                                                                   DJ = 0.0D0
116                                                                   DO 117 J = 1, M
                                                                        YRT(J) = -4.0D0*C2*(SIN((FLOAT(J) - DJ)*DY))**2
117                                                                     CONTINUE
                                                                        SCALY = 2.0D0*SCALY
                                                                        GO TO(124, 118, 122, 119, 123), MP
118                                                                     CALL SINTI(M, WY, ifac)
                                                                        GO TO 124
119                                                                     CALL COSTI(M, WY, ifac)
                                                                        GO TO 124
120                                                                     YRT(1) = 0.0D0
                                                                        YRT(M) = -4.0D0*C2
                                                                        DO 121 J = 3, M, 2
                                                                          YRT(J - 1) = -4.0D0*C2*(SIN(FLOAT((J - 1))*DY))**2
                                                                          YRT(J) = YRT(J - 1)
121                                                                       CONTINUE
                                                                          CALL RFFTI(M, WY, ifac)
                                                                          GO TO 124
                                                                          !122 CALL SINQI (M,WY,ifac)
122                                                                       CALL SINQI(M, WY, ifac)
                                                                          GO TO 124
123                                                                       CALL COSQI(M, WY, ifac)
124                                                                       CONTINUE
                                                                          IFWRD = 1
                                                                          IS = 1
125                                                                       CONTINUE
!
! TRANSFORM X
!
                                                                          DO 141 J = 1, M
                                                                            DO 140 K = 1, N
                                                                              DO 126 I = 1, L
                                                                                T(I) = F(I, J, K)
126                                                                             CONTINUE
                                                                                GO TO(127, 130, 131, 134, 135), LP
127                                                                             GO TO(128, 129), IFWRD
128                                                                             CALL RFFTF(L, T, WX, ifac)
                                                                                GO TO 138
129                                                                             CALL RFFTB(L, T, WX, ifac)
                                                                                GO TO 138
130                                                                             CALL SINT(L, T, WX, ifac)
                                                                                GO TO 138
131                                                                             GO TO(132, 133), IFWRD
132                                                                             CALL SINQF(L, T, WX, ifac)
                                                                                GO TO 138
133                                                                             CALL SINQB(L, T, WX, ifac)
                                                                                GO TO 138
134                                                                             CALL COST(L, T, WX, ifac)
                                                                                GO TO 138
135                                                                             GO TO(136, 137), IFWRD
136                                                                             CALL COSQF(L, T, WX, ifac)
                                                                                GO TO 138
137                                                                             CALL COSQB(L, T, WX, ifac)
138                                                                             CONTINUE
                                                                                DO 139 I = 1, L
                                                                                  F(I, J, K) = T(I)
139                                                                               CONTINUE
140                                                                               CONTINUE
141                                                                               CONTINUE
                                                                                  GO TO(142, 164), IFWRD
!
! TRANSFORM Y
!
142                                                                               CONTINUE
                                                                                  DO 158 I = 1, L
                                                                                    DO 157 K = 1, N
                                                                                      DO 143 J = 1, M
                                                                                        T(J) = F(I, J, K)
143                                                                                     CONTINUE
                                                                                        GO TO(144, 147, 148, 151, 152), MP
144                                                                                     GO TO(145, 146), IFWRD
145                                                                                     CALL RFFTF(M, T, WY, ifac)
                                                                                        GO TO 155
146                                                                                     CALL RFFTB(M, T, WY, ifac)
                                                                                        GO TO 155
147                                                                                     CALL SINT(M, T, WY, ifac)
                                                                                        GO TO 155
148                                                                                     GO TO(149, 150), IFWRD
149                                                                                     CALL SINQF(M, T, WY, ifac)
                                                                                        GO TO 155
150                                                                                     CALL SINQB(M, T, WY, ifac)
                                                                                        GO TO 155
151                                                                                     CALL COST(M, T, WY, ifac)
                                                                                        GO TO 155
152                                                                                     GO TO(153, 154), IFWRD
153                                                                                     CALL COSQF(M, T, WY, ifac)
                                                                                        GO TO 155
154                                                                                     CALL COSQB(M, T, WY, ifac)
155                                                                                     CONTINUE
                                                                                        DO 156 J = 1, M
                                                                                          F(I, J, K) = T(J)
156                                                                                       CONTINUE
157                                                                                       CONTINUE
158                                                                                       CONTINUE
                                                                                          GO TO(159, 125), IFWRD
159                                                                                       CONTINUE
!
! SOLVE TRIDIAGONAL SYSTEMS IN Z
!
                                                                                          DO 163 I = 1, L
                                                                                            DO 162 J = 1, M
                                                                                              DO 160 K = 1, N
                                                                                                BB(K) = B(K) + XRT(I) + YRT(J)
                                                                                                T(K) = F(I, J, K)
160                                                                                             CONTINUE
                                                                                                CALL TRID(N, A, BB, C, T, D)
                                                                                                DO 161 K = 1, N
                                                                                                  F(I, J, K) = T(K)
161                                                                                               CONTINUE
162                                                                                               CONTINUE
163                                                                                               CONTINUE
                                                                                                  IFWRD = 2
                                                                                                  IS = -1
                                                                                                  GO TO 142
164                                                                                               CONTINUE
                                                                                                  DO I = 1, L
                                                                                                    DO J = 1, M
                                                                                                      DO K = 1, N
                                                                                               F(I, J, K) = F(I, J, K)/(SCALX*SCALY)
                                                                                                      END DO
                                                                                                    END DO
                                                                                                  END DO
                                                                                                  RETURN
                                                                                                  END SUBROUTINE

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

                                                                                                  SUBROUTINE trid(m, a, b, c, y, d)
                                                                                                    USE params, ONLY: DP, PI
                                                                                                    IMPLICIT NONE

                                                                                                    INTEGER, INTENT(IN) :: m
                                                                                                    REAL(DP), INTENT(IN) :: a(:)
                                                                                                    REAL(DP), INTENT(IN) :: b(:)
                                                                                                    REAL(DP), INTENT(IN) :: c(:)
                                                                                                    REAL(DP), INTENT(INOUT) :: y(:)
                                                                                                    REAL(DP), INTENT(OUT) :: d(:)

                                                                                                    INTEGER :: i, ip, mm1
                                                                                                    REAL(DP):: z
                                                                                                    mm1 = m - 1
                                                                                                    z = 1.0D0/b(1)
                                                                                                    d(1) = c(1)*z
                                                                                                    y(1) = y(1)*z
                                                                                                    DO i = 2, mm1
                                                                                                    z = 1.0D0/(b(i) - a(i)*d(i - 1))
                                                                                                      d(i) = c(i)*z
                                                                                                     y(i) = (y(i) - a(i)*y(i - 1))*z
                                                                                                    END DO
                                                                                                    z = b(m) - a(m)*d(mm1)
                                                                                                    IF (Z .NE. 0.0D0) THEN
                                                                                                      y(m) = (y(m) - a(m)*y(mm1))/z
                                                                                                    ELSE
                                                                                                      y(m) = 0.0D0
                                                                                                    END IF
                                                                                                    DO ip = 1, mm1
                                                                                                      i = m - ip
                                                                                                      y(i) = y(i) - d(i)*y(i + 1)
                                                                                                    END DO
                                                                                                    RETURN
                                                                                                  END SUBROUTINE trid

                                                                                                  END MODULE fishpack
