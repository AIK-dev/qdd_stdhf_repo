
!-----COMPLEX FFT PACKAGE FROM NETLIB -----------------------------------------
!     CONVERTED TO F90

! VERSION FOR SINGLE PRECISION

!------------------------------------------------------------------------------
! PREPARATION:
! ADJUST THE KDIM PARAMETER TO INTENDED PROBLEM SIZE;
! IT NEEDS TO CHANGED TO THE SAME VALUE IN EACH SUBROUTINE.

!------------------------------------------------------------------------------
! USAGE:

! -> INITIALIZATION OF ARRAYS 'WA(1:KDIM)' AND 'IFAC(1:KDIM)'
! BY
! CALL DCFTI1 (N,WSAVE,IFAC)
! WITH
! N = DIMENSION OF THE TRANSFORMATION,
! MAKE SURE THAT 2*N < KDIM !
! WSAVE = STORAGE FOR TRIGONOMETRIC COEFFICIENTS
! IFAC = STORAGE FOR POINTERS
! NOTE: THE ARRAYS 'WA' AND 'IFAC' MUST BE 'SAVE'D IN THE CALLING
! ROUTINE.

! -> FORWARD TRANSFORMATION BY
! CALL DCFTF1 (N,C,WRKFFT,WSAVE,IFAC)
! WITH
! N = DIMENSION OF THE TRANSFORMATION,
! MAKE SURE THAT 2*N < KDIM !
! C = COMPLEX ARRAY OF LENGTH 'N' TO BE TRANSFORMED
! WRKFFT = WORKSPACE, REAL ARRAY(1:2*N)
! WSAVE = STORAGE FOR TRIGONOMETRIC COEFFICIENTS
! IFAC = STORAGE FOR POINTERS

! -> BACKWARD TRANSFORMATION BY
! CALL DCFTB1 (N,C,WRKFFT,WSAVE,IFAC)
! WITH
! N = DIMENSION OF THE TRANSFORMATION,
! MAKE SURE THAT 2*N < KDIM !
! C = COMPLEX ARRAY OF LENGTH 'N' TO BE TRANSFORMED
! WRKFFT = WORKSPACE, REAL ARRAY(1:2*N)
! WSAVE = STORAGE FOR TRIGONOMETRIC COEFFICIENTS
! IFAC = STORAGE FOR POINTERS
! NOTE: RESULTING 'C' NEEDS YET TO BE NORMALIZED BY '1/N',
! THEN WE HAVE EXACT REPRODUCTION FROM SUBSEQUENT
! FORWARD AND BACKWARD TRANSFORMATIONS.

! -> STORAGE: AS USUAL IN FFT, CENTER POINT IS AT 'C(1)'; 'X' OR
! 'P' < 0 ARE FROM 'C(N)' DOWNWARD; BOUNDARY OF
! THE TRANSFORMATION INTERVALS IS AT 'C(N/2)'.

!------------------------------------------------------------------------------

SUBROUTINE dcfti1(n, wa, ifac)
  USE params, ONLY: DP, tPI
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: wa(kdim)
  INTEGER, INTENT(OUT) :: ifac(kdim)
  REAL(DP) :: arg, argh, argld, fi
  INTEGER :: i, i1, ib, ido, idot, ii, ip, ipm
  INTEGER :: j, k1, l1, l2, ld, nf, nl, nq, nr, ntry
  INTEGER :: ntryh(4)
  DATA ntryh(1), ntryh(2), ntryh(3), ntryh(4)/3, 4, 2, 5/

  nl = n
  nf = 0
  j = 1
  ntry = ntryh(j)

  DO WHILE (nl /= 1)
    nq = nl/ntry
    nr = nl - ntry*nq
    IF (nr /= 0) THEN
      j = j + 1
      IF (j <= 4) THEN
        ntry = ntryh(j)
      ELSE
        ntry = ntry + 2
      END IF
      CYCLE
    END IF

    nf = nf + 1
    ifac(nf + 2) = ntry
    nl = nq
    IF ((ntry == 2) .AND. (nf /= 1)) THEN
      DO i = 2, nf
        ib = nf - i + 2
        ifac(ib + 2) = ifac(ib + 1)
      END DO
      ifac(3) = 2
    END IF
  END DO

  ifac(1) = n
  ifac(2) = nf

  argh = tpi/REAL(n, DP)
  i = 2
  l1 = 1
  DO k1 = 1, nf
    ip = ifac(k1 + 2)
    ld = 0
    l2 = l1*ip
    ido = n/l2
    idot = ido + ido + 2
    ipm = ip - 1

    DO j = 1, ipm
      i1 = i
      wa(i - 1) = 1.d0
      wa(i) = 0.d0
      ld = ld + l1
      fi = 0.d0
      argld = (ld)*argh
      DO ii = 4, idot, 2
        i = i + 2
        fi = fi + 1.d0
        arg = fi*argld
        wa(i - 1) = COS(arg)
        wa(i) = SIN(arg)
      END DO
      IF (ip <= 5) CYCLE
      wa(i1 - 1) = wa(i - 1)
      wa(i1) = wa(i)
    END DO

    l1 = l2
  END DO

  RETURN
END SUBROUTINE dcfti1

SUBROUTINE dcftf1(n, c, ch, wa, ifac)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: c(kdim)
  REAL(DP), INTENT(IN OUT) :: ch(kdim)
  REAL(DP), INTENT(IN OUT) :: wa(kdim)
  INTEGER, INTENT(IN) :: ifac(kdim)

  INTEGER :: idl1, ido, idot, ip, iw, ix2, ix3, ix4
  INTEGER :: k1, l1, l2, na, nac, nf
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1
  DO k1 = 1, nf
    ip = ifac(k1 + 2)
    l2 = ip*l1
    ido = n/l2
    idot = ido + ido
    idl1 = idot*l1
    SELECT CASE (ip)
    CASE (4)
      ix2 = iw + idot
      ix3 = ix2 + idot
      IF (na /= 0) THEN
        CALL dpssf4(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
      ELSE
        CALL dpssf4(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
      END IF
      na = 1 - na
    CASE (2)
      IF (na /= 0) THEN
        CALL dpssf2(idot, l1, ch, c, wa(iw))
      ELSE
        CALL dpssf2(idot, l1, c, ch, wa(iw))
      END IF
      na = 1 - na
    CASE (3)
      ix2 = iw + idot
      IF (na /= 0) THEN
        CALL dpssf3(idot, l1, ch, c, wa(iw), wa(ix2))
      ELSE
        CALL dpssf3(idot, l1, c, ch, wa(iw), wa(ix2))
      END IF
      na = 1 - na
    CASE (5)
      ix2 = iw + idot
      ix3 = ix2 + idot
      ix4 = ix3 + idot
      IF (na /= 0) THEN
        CALL dpssf5(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      ELSE
        CALL dpssf5(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      END IF
      na = 1 - na
    CASE DEFAULT
      IF (na /= 0) THEN
        CALL dpssf(nac, idot, ip, l1, idl1, ch, ch, ch, c, c, wa(iw))
      ELSE
        CALL dpssf(nac, idot, ip, l1, idl1, c, c, c, ch, ch, wa(iw))
      END IF
      IF (nac /= 0) na = 1 - na
    END SELECT

    l1 = l2
    iw = iw + (ip - 1)*idot
  END DO

  IF (na /= 0) c(1:2*n) = ch(1:2*n)

  RETURN
END SUBROUTINE dcftf1

SUBROUTINE dcftb1(n, c, ch, wa, ifac)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: n
  REAL(DP), INTENT(OUT) :: c(kdim)
  REAL(DP), INTENT(IN OUT) :: ch(kdim)
  REAL(DP), INTENT(IN OUT) :: wa(kdim)
  INTEGER, INTENT(IN) :: ifac(kdim)

  INTEGER :: idl1, ido, idot, ip, iw, ix2, ix3, ix4
  INTEGER :: k1, l1, l2, na, nac, nf

  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1
  DO k1 = 1, nf
    ip = ifac(k1 + 2)
    l2 = ip*l1
    ido = n/l2
    idot = ido + ido
    idl1 = idot*l1
    SELECT CASE (ip)
    CASE (4)
      ix2 = iw + idot
      ix3 = ix2 + idot
      IF (na /= 0) THEN
        CALL dpssb4(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3))
      ELSE
        CALL dpssb4(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3))
      END IF
      na = 1 - na

    CASE (2)
      IF (na /= 0) THEN
        CALL dpssb2(idot, l1, ch, c, wa(iw))
      ELSE
        CALL dpssb2(idot, l1, c, ch, wa(iw))
      END IF
      na = 1 - na

    CASE (3)
      ix2 = iw + idot
      IF (na /= 0) THEN
        CALL dpssb3(idot, l1, ch, c, wa(iw), wa(ix2))
      ELSE
        CALL dpssb3(idot, l1, c, ch, wa(iw), wa(ix2))
      END IF
      na = 1 - na

    CASE (5)
      ix2 = iw + idot
      ix3 = ix2 + idot
      ix4 = ix3 + idot
      IF (na /= 0) THEN
        CALL dpssb5(idot, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      ELSE
        CALL dpssb5(idot, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4))
      END IF
      na = 1 - na

    CASE DEFAULT
      IF (na /= 0) THEN
        CALL dpssb(nac, idot, ip, l1, idl1, ch, ch, ch, c, c, wa(iw))
      ELSE
        CALL dpssb(nac, idot, ip, l1, idl1, c, c, c, ch, ch, wa(iw))
      END IF
      IF (nac /= 0) na = 1 - na

    END SELECT
    l1 = l2
    iw = iw + (ip - 1)*idot
  END DO

  IF (na /= 0) c(1:2*n) = ch(1:2*n)

  RETURN
END SUBROUTINE dcftb1

SUBROUTINE dpssb(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(OUT) :: nac
  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: ip
  INTEGER, INTENT(IN) :: l1
  INTEGER, INTENT(IN) :: idl1
  REAL(DP), INTENT(IN) :: cc(ido, ip, l1)
  REAL(DP), INTENT(OUT) :: c1(ido, l1, ip)
  REAL(DP), INTENT(OUT) :: c2(idl1, ip)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, ip)
  REAL(DP), INTENT(IN OUT) :: ch2(idl1, ip)
  REAL(DP), INTENT(IN) :: wa(kdim)

  REAL(DP) :: wai, war
  INTEGER :: i, idij, idj, idl, idlj, idot, idp, ik, inc, ipp2, ipph
  INTEGER :: j, jc, k, l, lc, nt

  idot = ido/2
  nt = ip*idl1
  ipp2 = ip + 2
  ipph = (ip + 1)/2
  idp = ip*ido

  IF (ido < l1) THEN

    DO j = 2, ipph
      jc = ipp2 - j
      DO i = 1, ido
        DO k = 1, l1
          ch(i, k, j) = cc(i, j, k) + cc(i, jc, k)
          ch(i, k, jc) = cc(i, j, k) - cc(i, jc, k)
        END DO
      END DO
    END DO

    DO i = 1, ido
      DO k = 1, l1
        ch(i, k, 1) = cc(i, 1, k)
      END DO
    END DO

  ELSE

    DO j = 2, ipph
      jc = ipp2 - j
      DO k = 1, l1
        DO i = 1, ido
          ch(i, k, j) = cc(i, j, k) + cc(i, jc, k)
          ch(i, k, jc) = cc(i, j, k) - cc(i, jc, k)
        END DO
      END DO
    END DO

    DO k = 1, l1
      DO i = 1, ido
        ch(i, k, 1) = cc(i, 1, k)
      END DO
    END DO

  END IF

  idl = 2 - ido
  inc = 0
  DO l = 2, ipph
    lc = ipp2 - l
    idl = idl + ido
    DO ik = 1, idl1
      c2(ik, l) = ch2(ik, 1) + wa(idl - 1)*ch2(ik, 2)
      c2(ik, lc) = wa(idl)*ch2(ik, ip)
    END DO
    idlj = idl
    inc = inc + ido
    DO j = 3, ipph
      jc = ipp2 - j
      idlj = idlj + inc
      IF (idlj > idp) idlj = idlj - idp
      war = wa(idlj - 1)
      wai = wa(idlj)
      DO ik = 1, idl1
        c2(ik, l) = c2(ik, l) + war*ch2(ik, j)
        c2(ik, lc) = c2(ik, lc) + wai*ch2(ik, jc)
      END DO
    END DO
  END DO

  DO j = 2, ipph
    DO ik = 1, idl1
      ch2(ik, 1) = ch2(ik, 1) + ch2(ik, j)
    END DO
  END DO

  DO j = 2, ipph
    jc = ipp2 - j
    DO ik = 2, idl1, 2
      ch2(ik - 1, j) = c2(ik - 1, j) - c2(ik, jc)
      ch2(ik - 1, jc) = c2(ik - 1, j) + c2(ik, jc)
      ch2(ik, j) = c2(ik, j) + c2(ik - 1, jc)
      ch2(ik, jc) = c2(ik, j) - c2(ik - 1, jc)
    END DO
  END DO

  nac = 1
  IF (ido == 2) RETURN
  nac = 0

  DO ik = 1, idl1
    c2(ik, 1) = ch2(ik, 1)
  END DO

  DO j = 2, ip
    DO k = 1, l1
      c1(1, k, j) = ch(1, k, j)
      c1(2, k, j) = ch(2, k, j)
    END DO
  END DO

  IF (idot > l1) THEN

    idj = 2 - ido
    DO j = 2, ip
      idj = idj + ido
      DO k = 1, l1
        idij = idj
        DO i = 4, ido, 2
          idij = idij + 2
          c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - wa(idij)*ch(i, k, j)
          c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij)*ch(i - 1, k, j)
        END DO
      END DO
    END DO

  ELSE

    idij = 0
    DO j = 2, ip
      idij = idij + 2
      DO i = 4, ido, 2
        idij = idij + 2
        DO k = 1, l1
          c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - wa(idij)*ch(i, k, j)
          c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij)*ch(i - 1, k, j)
        END DO
      END DO
    END DO

  END IF

  RETURN

END SUBROUTINE dpssb

SUBROUTINE dpssb2(ido, l1, cc, ch, wa1)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 2, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 2)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  INTEGER :: i, k
  REAL(DP) :: ti2, tr2

  IF (ido > 2) THEN
    DO k = 1, l1
      DO i = 2, ido, 2
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + cc(i - 1, 2, k)
        tr2 = cc(i - 1, 1, k) - cc(i - 1, 2, k)
        ch(i, k, 1) = cc(i, 1, k) + cc(i, 2, k)
        ti2 = cc(i, 1, k) - cc(i, 2, k)
        ch(i, k, 2) = wa1(i - 1)*ti2 + wa1(i)*tr2
        ch(i - 1, k, 2) = wa1(i - 1)*tr2 - wa1(i)*ti2
      END DO
    END DO
  ELSE
    DO k = 1, l1
      ch(1, k, 1) = cc(1, 1, k) + cc(1, 2, k)
      ch(1, k, 2) = cc(1, 1, k) - cc(1, 2, k)
      ch(2, k, 1) = cc(2, 1, k) + cc(2, 2, k)
      ch(2, k, 2) = cc(2, 1, k) - cc(2, 2, k)
    END DO
  END IF

  RETURN
END SUBROUTINE dpssb2

SUBROUTINE dpssb3(ido, l1, cc, ch, wa1, wa2)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 3, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 3)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP) :: ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui, taur, ti2, tr2
  INTEGER :: i, k
  DATA taur/-0.5D0/
  DATA taui/0.86602540378443864676372317075293618D0/

! ONE HALF SQRT(3) = .866025..... .

  IF (ido /= 2) THEN
    DO k = 1, l1
      DO i = 2, ido, 2
        tr2 = cc(i - 1, 2, k) + cc(i - 1, 3, k)
        cr2 = cc(i - 1, 1, k) + taur*tr2
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2
        ti2 = cc(i, 2, k) + cc(i, 3, k)
        ci2 = cc(i, 1, k) + taur*ti2
        ch(i, k, 1) = cc(i, 1, k) + ti2
        cr3 = taui*(cc(i - 1, 2, k) - cc(i - 1, 3, k))
        ci3 = taui*(cc(i, 2, k) - cc(i, 3, k))
        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3
        ch(i, k, 2) = wa1(i - 1)*di2 + wa1(i)*dr2
        ch(i - 1, k, 2) = wa1(i - 1)*dr2 - wa1(i)*di2
        ch(i, k, 3) = wa2(i - 1)*di3 + wa2(i)*dr3
        ch(i - 1, k, 3) = wa2(i - 1)*dr3 - wa2(i)*di3
      END DO
    END DO
  ELSE
    DO k = 1, l1
      tr2 = cc(1, 2, k) + cc(1, 3, k)
      cr2 = cc(1, 1, k) + taur*tr2
      ch(1, k, 1) = cc(1, 1, k) + tr2
      ti2 = cc(2, 2, k) + cc(2, 3, k)
      ci2 = cc(2, 1, k) + taur*ti2
      ch(2, k, 1) = cc(2, 1, k) + ti2
      cr3 = taui*(cc(1, 2, k) - cc(1, 3, k))
      ci3 = taui*(cc(2, 2, k) - cc(2, 3, k))
      ch(1, k, 2) = cr2 - ci3
      ch(1, k, 3) = cr2 + ci3
      ch(2, k, 2) = ci2 + cr3
      ch(2, k, 3) = ci2 - cr3
    END DO
  END IF

  RETURN
END SUBROUTINE dpssb3

SUBROUTINE dpssb4(ido, l1, cc, ch, wa1, wa2, wa3)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 4, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 4)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP), INTENT(IN) :: wa3(kdim)
  REAL(DP) :: ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4

  INTEGER :: i, k

  IF (ido /= 2) THEN
    DO k = 1, l1
      DO i = 2, ido, 2
        ti1 = cc(i, 1, k) - cc(i, 3, k)
        ti2 = cc(i, 1, k) + cc(i, 3, k)
        ti3 = cc(i, 2, k) + cc(i, 4, k)
        tr4 = cc(i, 4, k) - cc(i, 2, k)
        tr1 = cc(i - 1, 1, k) - cc(i - 1, 3, k)
        tr2 = cc(i - 1, 1, k) + cc(i - 1, 3, k)
        ti4 = cc(i - 1, 2, k) - cc(i - 1, 4, k)
        tr3 = cc(i - 1, 2, k) + cc(i - 1, 4, k)
        ch(i - 1, k, 1) = tr2 + tr3
        cr3 = tr2 - tr3
        ch(i, k, 1) = ti2 + ti3
        ci3 = ti2 - ti3
        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4
        ch(i - 1, k, 2) = wa1(i - 1)*cr2 - wa1(i)*ci2
        ch(i, k, 2) = wa1(i - 1)*ci2 + wa1(i)*cr2
        ch(i - 1, k, 3) = wa2(i - 1)*cr3 - wa2(i)*ci3
        ch(i, k, 3) = wa2(i - 1)*ci3 + wa2(i)*cr3
        ch(i - 1, k, 4) = wa3(i - 1)*cr4 - wa3(i)*ci4
        ch(i, k, 4) = wa3(i - 1)*ci4 + wa3(i)*cr4
      END DO
    END DO
  ELSE
    DO k = 1, l1
      ti1 = cc(2, 1, k) - cc(2, 3, k)
      ti2 = cc(2, 1, k) + cc(2, 3, k)
      tr4 = cc(2, 4, k) - cc(2, 2, k)
      ti3 = cc(2, 2, k) + cc(2, 4, k)
      tr1 = cc(1, 1, k) - cc(1, 3, k)
      tr2 = cc(1, 1, k) + cc(1, 3, k)
      ti4 = cc(1, 2, k) - cc(1, 4, k)
      tr3 = cc(1, 2, k) + cc(1, 4, k)
      ch(1, k, 1) = tr2 + tr3
      ch(1, k, 3) = tr2 - tr3
      ch(2, k, 1) = ti2 + ti3
      ch(2, k, 3) = ti2 - ti3
      ch(1, k, 2) = tr1 + tr4
      ch(1, k, 4) = tr1 - tr4
      ch(2, k, 2) = ti1 + ti4
      ch(2, k, 4) = ti1 - ti4
    END DO
  END IF

  RETURN
END SUBROUTINE dpssb4

SUBROUTINE dpssb5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 5, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 5)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP), INTENT(IN) :: wa3(kdim)
  REAL(DP), INTENT(IN) :: wa4(kdim)
  REAL(DP) :: ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, &
              di2, di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3, &
              ti4, ti5, tr11, tr12, tr2, tr3, tr4, tr5
  INTEGER ::i, k
  DATA tr11/0.30901699437494742410229341718281906D0/
  DATA ti11/0.95105651629515357211643933337938214D0/
  DATA tr12/-0.80901699437494742410229341718281906D0/
  DATA ti12/0.58778525229247312916870595463907277D0/

! SIN(PI/10) = .30901699.... .
! COS(PI/10) = .95105651.... .
! SIN(PI/5 ) = .58778525.... .
! COS(PI/5 ) = .80901699.... .

  IF (ido /= 2) THEN
    DO k = 1, l1
      DO i = 2, ido, 2
        ti5 = cc(i, 2, k) - cc(i, 5, k)
        ti2 = cc(i, 2, k) + cc(i, 5, k)
        ti4 = cc(i, 3, k) - cc(i, 4, k)
        ti3 = cc(i, 3, k) + cc(i, 4, k)
        tr5 = cc(i - 1, 2, k) - cc(i - 1, 5, k)
        tr2 = cc(i - 1, 2, k) + cc(i - 1, 5, k)
        tr4 = cc(i - 1, 3, k) - cc(i - 1, 4, k)
        tr3 = cc(i - 1, 3, k) + cc(i - 1, 4, k)
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2 + tr3
        ch(i, k, 1) = cc(i, 1, k) + ti2 + ti3
        cr2 = cc(i - 1, 1, k) + tr11*tr2 + tr12*tr3
        ci2 = cc(i, 1, k) + tr11*ti2 + tr12*ti3
        cr3 = cc(i - 1, 1, k) + tr12*tr2 + tr11*tr3
        ci3 = cc(i, 1, k) + tr12*ti2 + tr11*ti3
        cr5 = ti11*tr5 + ti12*tr4
        ci5 = ti11*ti5 + ti12*ti4
        cr4 = ti12*tr5 - ti11*tr4
        ci4 = ti12*ti5 - ti11*ti4
        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5
        ch(i - 1, k, 2) = wa1(i - 1)*dr2 - wa1(i)*di2
        ch(i, k, 2) = wa1(i - 1)*di2 + wa1(i)*dr2
        ch(i - 1, k, 3) = wa2(i - 1)*dr3 - wa2(i)*di3
        ch(i, k, 3) = wa2(i - 1)*di3 + wa2(i)*dr3
        ch(i - 1, k, 4) = wa3(i - 1)*dr4 - wa3(i)*di4
        ch(i, k, 4) = wa3(i - 1)*di4 + wa3(i)*dr4
        ch(i - 1, k, 5) = wa4(i - 1)*dr5 - wa4(i)*di5
        ch(i, k, 5) = wa4(i - 1)*di5 + wa4(i)*dr5
      END DO
    END DO
  ELSE
    DO k = 1, l1
      ti5 = cc(2, 2, k) - cc(2, 5, k)
      ti2 = cc(2, 2, k) + cc(2, 5, k)
      ti4 = cc(2, 3, k) - cc(2, 4, k)
      ti3 = cc(2, 3, k) + cc(2, 4, k)
      tr5 = cc(1, 2, k) - cc(1, 5, k)
      tr2 = cc(1, 2, k) + cc(1, 5, k)
      tr4 = cc(1, 3, k) - cc(1, 4, k)
      tr3 = cc(1, 3, k) + cc(1, 4, k)
      ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
      ch(2, k, 1) = cc(2, 1, k) + ti2 + ti3
      cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
      ci2 = cc(2, 1, k) + tr11*ti2 + tr12*ti3
      cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
      ci3 = cc(2, 1, k) + tr12*ti2 + tr11*ti3
      cr5 = ti11*tr5 + ti12*tr4
      ci5 = ti11*ti5 + ti12*ti4
      cr4 = ti12*tr5 - ti11*tr4
      ci4 = ti12*ti5 - ti11*ti4
      ch(1, k, 2) = cr2 - ci5
      ch(1, k, 5) = cr2 + ci5
      ch(2, k, 2) = ci2 + cr5
      ch(2, k, 3) = ci3 + cr4
      ch(1, k, 3) = cr3 - ci4
      ch(1, k, 4) = cr3 + ci4
      ch(2, k, 4) = ci3 - cr4
      ch(2, k, 5) = ci2 - cr5
    END DO
  END IF

  RETURN
END SUBROUTINE dpssb5

SUBROUTINE dpssf(nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(OUT) :: nac
  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: ip
  INTEGER, INTENT(IN) :: l1
  INTEGER, INTENT(IN) :: idl1
  REAL(DP), INTENT(IN) :: cc(ido, ip, l1)
  REAL(DP), INTENT(OUT) :: c1(ido, l1, ip)
  REAL(DP), INTENT(OUT) :: c2(idl1, ip)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, ip)
  REAL(DP), INTENT(IN OUT) :: ch2(idl1, ip)
  REAL(DP), INTENT(IN) :: wa(kdim)
  REAL(DP) :: wai, war

  INTEGER :: i, idij, idj, idl, idlj, idot, idp, ik, inc, ipp2, ipph
  INTEGER :: j, jc, k, l, lc, nt

  idot = ido/2
  nt = ip*idl1
  ipp2 = ip + 2
  ipph = (ip + 1)/2
  idp = ip*ido

  IF (ido < l1) THEN
    DO j = 2, ipph
      jc = ipp2 - j
      DO i = 1, ido
        DO k = 1, l1
          ch(i, k, j) = cc(i, j, k) + cc(i, jc, k)
          ch(i, k, jc) = cc(i, j, k) - cc(i, jc, k)
        END DO
      END DO
    END DO

    DO i = 1, ido
      DO k = 1, l1
        ch(i, k, 1) = cc(i, 1, k)
      END DO
    END DO
  ELSE
    DO j = 2, ipph
      jc = ipp2 - j
      DO k = 1, l1
        DO i = 1, ido
          ch(i, k, j) = cc(i, j, k) + cc(i, jc, k)
          ch(i, k, jc) = cc(i, j, k) - cc(i, jc, k)
        END DO
      END DO
    END DO

    DO k = 1, l1
      DO i = 1, ido
        ch(i, k, 1) = cc(i, 1, k)
      END DO
    END DO
  END IF

  idl = 2 - ido
  inc = 0
  DO l = 2, ipph
    lc = ipp2 - l
    idl = idl + ido
    DO ik = 1, idl1
      c2(ik, l) = ch2(ik, 1) + wa(idl - 1)*ch2(ik, 2)
      c2(ik, lc) = -wa(idl)*ch2(ik, ip)
    END DO
    idlj = idl
    inc = inc + ido
    DO j = 3, ipph
      jc = ipp2 - j
      idlj = idlj + inc
      IF (idlj > idp) idlj = idlj - idp
      war = wa(idlj - 1)
      wai = wa(idlj)
      DO ik = 1, idl1
        c2(ik, l) = c2(ik, l) + war*ch2(ik, j)
        c2(ik, lc) = c2(ik, lc) - wai*ch2(ik, jc)
      END DO
    END DO
  END DO

  DO j = 2, ipph
    DO ik = 1, idl1
      ch2(ik, 1) = ch2(ik, 1) + ch2(ik, j)
    END DO
  END DO

  DO j = 2, ipph
    jc = ipp2 - j
    DO ik = 2, idl1, 2
      ch2(ik - 1, j) = c2(ik - 1, j) - c2(ik, jc)
      ch2(ik - 1, jc) = c2(ik - 1, j) + c2(ik, jc)
      ch2(ik, j) = c2(ik, j) + c2(ik - 1, jc)
      ch2(ik, jc) = c2(ik, j) - c2(ik - 1, jc)
    END DO
  END DO

  nac = 1
  IF (ido == 2) RETURN
  nac = 0

  DO ik = 1, idl1
    c2(ik, 1) = ch2(ik, 1)
  END DO

  DO j = 2, ip
    DO k = 1, l1
      c1(1, k, j) = ch(1, k, j)
      c1(2, k, j) = ch(2, k, j)
    END DO
  END DO

  IF (idot > l1) THEN
    idj = 2 - ido
    DO j = 2, ip
      idj = idj + ido
      DO k = 1, l1
        idij = idj
        DO i = 4, ido, 2
          idij = idij + 2
          c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) + wa(idij)*ch(i, k, j)
          c1(i, k, j) = wa(idij - 1)*ch(i, k, j) - wa(idij)*ch(i - 1, k, j)
        END DO
      END DO
    END DO
  ELSE
    idij = 0
    DO j = 2, ip
      idij = idij + 2
      DO i = 4, ido, 2
        idij = idij + 2
        DO k = 1, l1
          c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) + wa(idij)*ch(i, k, j)
          c1(i, k, j) = wa(idij - 1)*ch(i, k, j) - wa(idij)*ch(i - 1, k, j)
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE dpssf

SUBROUTINE dpssf2(ido, l1, cc, ch, wa1)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 2, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 2)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP) :: ti2, tr2
  INTEGER :: i, k
  IF (ido > 2) THEN
    DO k = 1, l1
      DO i = 2, ido, 2
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + cc(i - 1, 2, k)
        tr2 = cc(i - 1, 1, k) - cc(i - 1, 2, k)
        ch(i, k, 1) = cc(i, 1, k) + cc(i, 2, k)
        ti2 = cc(i, 1, k) - cc(i, 2, k)
        ch(i, k, 2) = wa1(i - 1)*ti2 - wa1(i)*tr2
        ch(i - 1, k, 2) = wa1(i - 1)*tr2 + wa1(i)*ti2
      END DO
    END DO
  ELSE
    DO k = 1, l1
      ch(1, k, 1) = cc(1, 1, k) + cc(1, 2, k)
      ch(1, k, 2) = cc(1, 1, k) - cc(1, 2, k)
      ch(2, k, 1) = cc(2, 1, k) + cc(2, 2, k)
      ch(2, k, 2) = cc(2, 1, k) - cc(2, 2, k)
    END DO
  END IF

  RETURN
END SUBROUTINE dpssf2

SUBROUTINE dpssf3(ido, l1, cc, ch, wa1, wa2)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 3, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 3)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP) :: ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui, taur, ti2, tr2
  INTEGER :: i, k
  DATA taur/-0.5D0/
  DATA taui/-0.86602540378443864676372317075293618D0/

  IF (ido /= 2) THEN
    DO k = 1, l1
      DO i = 2, ido, 2
        tr2 = cc(i - 1, 2, k) + cc(i - 1, 3, k)
        cr2 = cc(i - 1, 1, k) + taur*tr2
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2
        ti2 = cc(i, 2, k) + cc(i, 3, k)
        ci2 = cc(i, 1, k) + taur*ti2
        ch(i, k, 1) = cc(i, 1, k) + ti2
        cr3 = taui*(cc(i - 1, 2, k) - cc(i - 1, 3, k))
        ci3 = taui*(cc(i, 2, k) - cc(i, 3, k))
        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3
        ch(i, k, 2) = wa1(i - 1)*di2 - wa1(i)*dr2
        ch(i - 1, k, 2) = wa1(i - 1)*dr2 + wa1(i)*di2
        ch(i, k, 3) = wa2(i - 1)*di3 - wa2(i)*dr3
        ch(i - 1, k, 3) = wa2(i - 1)*dr3 + wa2(i)*di3
      END DO
    END DO
  ELSE
    DO k = 1, l1
      tr2 = cc(1, 2, k) + cc(1, 3, k)
      cr2 = cc(1, 1, k) + taur*tr2
      ch(1, k, 1) = cc(1, 1, k) + tr2
      ti2 = cc(2, 2, k) + cc(2, 3, k)
      ci2 = cc(2, 1, k) + taur*ti2
      ch(2, k, 1) = cc(2, 1, k) + ti2
      cr3 = taui*(cc(1, 2, k) - cc(1, 3, k))
      ci3 = taui*(cc(2, 2, k) - cc(2, 3, k))
      ch(1, k, 2) = cr2 - ci3
      ch(1, k, 3) = cr2 + ci3
      ch(2, k, 2) = ci2 + cr3
      ch(2, k, 3) = ci2 - cr3
    END DO
  END IF

  RETURN
END SUBROUTINE dpssf3

SUBROUTINE dpssf4(ido, l1, cc, ch, wa1, wa2, wa3)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 4, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 4)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP), INTENT(IN) :: wa3(kdim)
  REAL(DP) :: ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4
  INTEGER :: i, k
  IF (ido /= 2) THEN
    DO k = 1, l1
      DO i = 2, ido, 2
        ti1 = cc(i, 1, k) - cc(i, 3, k)
        ti2 = cc(i, 1, k) + cc(i, 3, k)
        ti3 = cc(i, 2, k) + cc(i, 4, k)
        tr4 = cc(i, 2, k) - cc(i, 4, k)
        tr1 = cc(i - 1, 1, k) - cc(i - 1, 3, k)
        tr2 = cc(i - 1, 1, k) + cc(i - 1, 3, k)
        ti4 = cc(i - 1, 4, k) - cc(i - 1, 2, k)
        tr3 = cc(i - 1, 2, k) + cc(i - 1, 4, k)
        ch(i - 1, k, 1) = tr2 + tr3
        cr3 = tr2 - tr3
        ch(i, k, 1) = ti2 + ti3
        ci3 = ti2 - ti3
        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4
        ch(i - 1, k, 2) = wa1(i - 1)*cr2 + wa1(i)*ci2
        ch(i, k, 2) = wa1(i - 1)*ci2 - wa1(i)*cr2
        ch(i - 1, k, 3) = wa2(i - 1)*cr3 + wa2(i)*ci3
        ch(i, k, 3) = wa2(i - 1)*ci3 - wa2(i)*cr3
        ch(i - 1, k, 4) = wa3(i - 1)*cr4 + wa3(i)*ci4
        ch(i, k, 4) = wa3(i - 1)*ci4 - wa3(i)*cr4
      END DO
    END DO
  ELSE
    DO k = 1, l1
      ti1 = cc(2, 1, k) - cc(2, 3, k)
      ti2 = cc(2, 1, k) + cc(2, 3, k)
      tr4 = cc(2, 2, k) - cc(2, 4, k)
      ti3 = cc(2, 2, k) + cc(2, 4, k)
      tr1 = cc(1, 1, k) - cc(1, 3, k)
      tr2 = cc(1, 1, k) + cc(1, 3, k)
      ti4 = cc(1, 4, k) - cc(1, 2, k)
      tr3 = cc(1, 2, k) + cc(1, 4, k)
      ch(1, k, 1) = tr2 + tr3
      ch(1, k, 3) = tr2 - tr3
      ch(2, k, 1) = ti2 + ti3
      ch(2, k, 3) = ti2 - ti3
      ch(1, k, 2) = tr1 + tr4
      ch(1, k, 4) = tr1 - tr4
      ch(2, k, 2) = ti1 + ti4
      ch(2, k, 4) = ti1 - ti4
    END DO
  END IF

  RETURN
END SUBROUTINE dpssf4

SUBROUTINE dpssf5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 5, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 5)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP), INTENT(IN) :: wa3(kdim)
  REAL(DP), INTENT(IN) :: wa4(kdim)
  REAL(DP) :: ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, &
              di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3, ti4, &
              ti5, tr11, tr12, tr2, tr3, tr4, tr5
  INTEGER :: i, k
  DATA tr11/0.30901699437494742410229341718281906D0/
  DATA ti11/-0.95105651629515357211643933337938214D0/
  DATA tr12/-0.80901699437494742410229341718281906D0/
  DATA ti12/-0.58778525229247312916870595463907277D0/

  IF (ido /= 2) THEN
    DO k = 1, l1
      DO i = 2, ido, 2
        ti5 = cc(i, 2, k) - cc(i, 5, k)
        ti2 = cc(i, 2, k) + cc(i, 5, k)
        ti4 = cc(i, 3, k) - cc(i, 4, k)
        ti3 = cc(i, 3, k) + cc(i, 4, k)
        tr5 = cc(i - 1, 2, k) - cc(i - 1, 5, k)
        tr2 = cc(i - 1, 2, k) + cc(i - 1, 5, k)
        tr4 = cc(i - 1, 3, k) - cc(i - 1, 4, k)
        tr3 = cc(i - 1, 3, k) + cc(i - 1, 4, k)
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2 + tr3
        ch(i, k, 1) = cc(i, 1, k) + ti2 + ti3
        cr2 = cc(i - 1, 1, k) + tr11*tr2 + tr12*tr3
        ci2 = cc(i, 1, k) + tr11*ti2 + tr12*ti3
        cr3 = cc(i - 1, 1, k) + tr12*tr2 + tr11*tr3
        ci3 = cc(i, 1, k) + tr12*ti2 + tr11*ti3
        cr5 = ti11*tr5 + ti12*tr4
        ci5 = ti11*ti5 + ti12*ti4
        cr4 = ti12*tr5 - ti11*tr4
        ci4 = ti12*ti5 - ti11*ti4
        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5
        ch(i - 1, k, 2) = wa1(i - 1)*dr2 + wa1(i)*di2
        ch(i, k, 2) = wa1(i - 1)*di2 - wa1(i)*dr2
        ch(i - 1, k, 3) = wa2(i - 1)*dr3 + wa2(i)*di3
        ch(i, k, 3) = wa2(i - 1)*di3 - wa2(i)*dr3
        ch(i - 1, k, 4) = wa3(i - 1)*dr4 + wa3(i)*di4
        ch(i, k, 4) = wa3(i - 1)*di4 - wa3(i)*dr4
        ch(i - 1, k, 5) = wa4(i - 1)*dr5 + wa4(i)*di5
        ch(i, k, 5) = wa4(i - 1)*di5 - wa4(i)*dr5
      END DO
    END DO
  ELSE
    DO k = 1, l1
      ti5 = cc(2, 2, k) - cc(2, 5, k)
      ti2 = cc(2, 2, k) + cc(2, 5, k)
      ti4 = cc(2, 3, k) - cc(2, 4, k)
      ti3 = cc(2, 3, k) + cc(2, 4, k)
      tr5 = cc(1, 2, k) - cc(1, 5, k)
      tr2 = cc(1, 2, k) + cc(1, 5, k)
      tr4 = cc(1, 3, k) - cc(1, 4, k)
      tr3 = cc(1, 3, k) + cc(1, 4, k)
      ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
      ch(2, k, 1) = cc(2, 1, k) + ti2 + ti3
      cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
      ci2 = cc(2, 1, k) + tr11*ti2 + tr12*ti3
      cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
      ci3 = cc(2, 1, k) + tr12*ti2 + tr11*ti3
      cr5 = ti11*tr5 + ti12*tr4
      ci5 = ti11*ti5 + ti12*ti4
      cr4 = ti12*tr5 - ti11*tr4
      ci4 = ti12*ti5 - ti11*ti4
      ch(1, k, 2) = cr2 - ci5
      ch(1, k, 5) = cr2 + ci5
      ch(2, k, 2) = ci2 + cr5
      ch(2, k, 3) = ci3 + cr4
      ch(1, k, 3) = cr3 - ci4
      ch(1, k, 4) = cr3 + ci4
      ch(2, k, 4) = ci3 - cr4
      ch(2, k, 5) = ci2 - cr5
    END DO
  END IF

  RETURN
END SUBROUTINE dpssf5

SUBROUTINE dradb2(ido, l1, cc, ch, wa1)
  USE params, ONLY: DP
  IMPLICIT NONE
  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 2, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 2)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP) :: ti2, tr2
  INTEGER :: i, ic, idp2, k
  DO k = 1, l1
    ch(1, k, 1) = cc(1, 1, k) + cc(ido, 2, k)
    ch(1, k, 2) = cc(1, 1, k) - cc(ido, 2, k)
  END DO

  IF (ido < 2) THEN
    RETURN
  ELSE IF (ido /= 2) THEN
    idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        ch(i - 1, k, 1) = cc(i - 1, 1, k) + cc(ic - 1, 2, k)
        tr2 = cc(i - 1, 1, k) - cc(ic - 1, 2, k)
        ch(i, k, 1) = cc(i, 1, k) - cc(ic, 2, k)
        ti2 = cc(i, 1, k) + cc(ic, 2, k)
        ch(i - 1, k, 2) = wa1(i - 2)*tr2 - wa1(i - 1)*ti2
        ch(i, k, 2) = wa1(i - 2)*ti2 + wa1(i - 1)*tr2
      END DO
    END DO
    IF (MOD(ido, 2) == 1) RETURN
  END IF

  DO k = 1, l1
    ch(ido, k, 1) = cc(ido, 1, k) + cc(ido, 1, k)
    ch(ido, k, 2) = -(cc(1, 2, k) + cc(1, 2, k))
  END DO

  RETURN
END SUBROUTINE dradb2

SUBROUTINE dradb3(ido, l1, cc, ch, wa1, wa2)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 3, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 3)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP) :: ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui, taur, ti2, tr2
  INTEGER :: i, ic, idp2, k
  DATA taur/-0.5D0/
  DATA taui/0.86602540378443864676372317075293618D0/

  DO k = 1, l1
    tr2 = cc(ido, 2, k) + cc(ido, 2, k)
    cr2 = cc(1, 1, k) + taur*tr2
    ch(1, k, 1) = cc(1, 1, k) + tr2
    ci3 = taui*(cc(1, 3, k) + cc(1, 3, k))
    ch(1, k, 2) = cr2 - ci3
    ch(1, k, 3) = cr2 + ci3
  END DO

  IF (ido == 1) RETURN
  idp2 = ido + 2
  DO k = 1, l1
    DO i = 3, ido, 2
      ic = idp2 - i
      tr2 = cc(i - 1, 3, k) + cc(ic - 1, 2, k)
      cr2 = cc(i - 1, 1, k) + taur*tr2
      ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2
      ti2 = cc(i, 3, k) - cc(ic, 2, k)
      ci2 = cc(i, 1, k) + taur*ti2
      ch(i, k, 1) = cc(i, 1, k) + ti2
      cr3 = taui*(cc(i - 1, 3, k) - cc(ic - 1, 2, k))
      ci3 = taui*(cc(i, 3, k) + cc(ic, 2, k))
      dr2 = cr2 - ci3
      dr3 = cr2 + ci3
      di2 = ci2 + cr3
      di3 = ci2 - cr3
      ch(i - 1, k, 2) = wa1(i - 2)*dr2 - wa1(i - 1)*di2
      ch(i, k, 2) = wa1(i - 2)*di2 + wa1(i - 1)*dr2
      ch(i - 1, k, 3) = wa2(i - 2)*dr3 - wa2(i - 1)*di3
      ch(i, k, 3) = wa2(i - 2)*di3 + wa2(i - 1)*dr3
    END DO
  END DO

  RETURN
END SUBROUTINE dradb3

SUBROUTINE dradb4(ido, l1, cc, ch, wa1, wa2, wa3)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 4, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 4)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP), INTENT(IN) :: wa3(kdim)
  REAL(DP) :: ci2, ci3, ci4, cr2, cr3, cr4, sqrt2, ti1, ti2, ti3, ti4, &
              tr1, tr2, tr3, tr4
  INTEGER :: i, ic, idp2, k
  DATA sqrt2/1.41421356237309504880168872420970D0/

  DO k = 1, l1
    tr1 = cc(1, 1, k) - cc(ido, 4, k)
    tr2 = cc(1, 1, k) + cc(ido, 4, k)
    tr3 = cc(ido, 2, k) + cc(ido, 2, k)
    tr4 = cc(1, 3, k) + cc(1, 3, k)
    ch(1, k, 1) = tr2 + tr3
    ch(1, k, 2) = tr1 - tr4
    ch(1, k, 3) = tr2 - tr3
    ch(1, k, 4) = tr1 + tr4
  END DO

  IF (ido < 2) THEN
    RETURN
  ELSE IF (ido /= 2) THEN
    idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        ti1 = cc(i, 1, k) + cc(ic, 4, k)
        ti2 = cc(i, 1, k) - cc(ic, 4, k)
        ti3 = cc(i, 3, k) - cc(ic, 2, k)
        tr4 = cc(i, 3, k) + cc(ic, 2, k)
        tr1 = cc(i - 1, 1, k) - cc(ic - 1, 4, k)
        tr2 = cc(i - 1, 1, k) + cc(ic - 1, 4, k)
        ti4 = cc(i - 1, 3, k) - cc(ic - 1, 2, k)
        tr3 = cc(i - 1, 3, k) + cc(ic - 1, 2, k)
        ch(i - 1, k, 1) = tr2 + tr3
        cr3 = tr2 - tr3
        ch(i, k, 1) = ti2 + ti3
        ci3 = ti2 - ti3
        cr2 = tr1 - tr4
        cr4 = tr1 + tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4
        ch(i - 1, k, 2) = wa1(i - 2)*cr2 - wa1(i - 1)*ci2
        ch(i, k, 2) = wa1(i - 2)*ci2 + wa1(i - 1)*cr2
        ch(i - 1, k, 3) = wa2(i - 2)*cr3 - wa2(i - 1)*ci3
        ch(i, k, 3) = wa2(i - 2)*ci3 + wa2(i - 1)*cr3
        ch(i - 1, k, 4) = wa3(i - 2)*cr4 - wa3(i - 1)*ci4
        ch(i, k, 4) = wa3(i - 2)*ci4 + wa3(i - 1)*cr4
      END DO
    END DO
    IF (MOD(ido, 2) == 1) RETURN
  END IF

  CONTINUE
  DO k = 1, l1
    ti1 = cc(1, 2, k) + cc(1, 4, k)
    ti2 = cc(1, 4, k) - cc(1, 2, k)
    tr1 = cc(ido, 1, k) - cc(ido, 3, k)
    tr2 = cc(ido, 1, k) + cc(ido, 3, k)
    ch(ido, k, 1) = tr2 + tr2
    ch(ido, k, 2) = sqrt2*(tr1 - ti1)
    ch(ido, k, 3) = ti2 + ti2
    ch(ido, k, 4) = -sqrt2*(tr1 + ti1)
  END DO

  RETURN
END SUBROUTINE dradb4

SUBROUTINE dradb5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, 5, l1)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, 5)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP), INTENT(IN) :: wa3(kdim)
  REAL(DP), INTENT(IN) :: wa4(kdim)
  REAL(DP) :: ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, &
              di2, di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3, &
              ti4, ti5, tr11, tr12, tr2, tr3, tr4, tr5
  INTEGER :: i, ic, idp2, k
  DATA tr11/0.30901699437494742410229341718281906D0/
  DATA ti11/0.95105651629515357211643933337938214D0/
  DATA tr12/-0.80901699437494742410229341718281906D0/
  DATA ti12/0.58778525229247312916870595463907277D0/

  DO k = 1, l1
    ti5 = cc(1, 3, k) + cc(1, 3, k)
    ti4 = cc(1, 5, k) + cc(1, 5, k)
    tr2 = cc(ido, 2, k) + cc(ido, 2, k)
    tr3 = cc(ido, 4, k) + cc(ido, 4, k)
    ch(1, k, 1) = cc(1, 1, k) + tr2 + tr3
    cr2 = cc(1, 1, k) + tr11*tr2 + tr12*tr3
    cr3 = cc(1, 1, k) + tr12*tr2 + tr11*tr3
    ci5 = ti11*ti5 + ti12*ti4
    ci4 = ti12*ti5 - ti11*ti4
    ch(1, k, 2) = cr2 - ci5
    ch(1, k, 3) = cr3 - ci4
    ch(1, k, 4) = cr3 + ci4
    ch(1, k, 5) = cr2 + ci5
  END DO
  IF (ido == 1) RETURN

  idp2 = ido + 2
  DO k = 1, l1
    DO i = 3, ido, 2
      ic = idp2 - i
      ti5 = cc(i, 3, k) + cc(ic, 2, k)
      ti2 = cc(i, 3, k) - cc(ic, 2, k)
      ti4 = cc(i, 5, k) + cc(ic, 4, k)
      ti3 = cc(i, 5, k) - cc(ic, 4, k)
      tr5 = cc(i - 1, 3, k) - cc(ic - 1, 2, k)
      tr2 = cc(i - 1, 3, k) + cc(ic - 1, 2, k)
      tr4 = cc(i - 1, 5, k) - cc(ic - 1, 4, k)
      tr3 = cc(i - 1, 5, k) + cc(ic - 1, 4, k)
      ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2 + tr3
      ch(i, k, 1) = cc(i, 1, k) + ti2 + ti3
      cr2 = cc(i - 1, 1, k) + tr11*tr2 + tr12*tr3
      ci2 = cc(i, 1, k) + tr11*ti2 + tr12*ti3
      cr3 = cc(i - 1, 1, k) + tr12*tr2 + tr11*tr3
      ci3 = cc(i, 1, k) + tr12*ti2 + tr11*ti3
      cr5 = ti11*tr5 + ti12*tr4
      ci5 = ti11*ti5 + ti12*ti4
      cr4 = ti12*tr5 - ti11*tr4
      ci4 = ti12*ti5 - ti11*ti4
      dr3 = cr3 - ci4
      dr4 = cr3 + ci4
      di3 = ci3 + cr4
      di4 = ci3 - cr4
      dr5 = cr2 + ci5
      dr2 = cr2 - ci5
      di5 = ci2 - cr5
      di2 = ci2 + cr5
      ch(i - 1, k, 2) = wa1(i - 2)*dr2 - wa1(i - 1)*di2
      ch(i, k, 2) = wa1(i - 2)*di2 + wa1(i - 1)*dr2
      ch(i - 1, k, 3) = wa2(i - 2)*dr3 - wa2(i - 1)*di3
      ch(i, k, 3) = wa2(i - 2)*di3 + wa2(i - 1)*dr3
      ch(i - 1, k, 4) = wa3(i - 2)*dr4 - wa3(i - 1)*di4
      ch(i, k, 4) = wa3(i - 2)*di4 + wa3(i - 1)*dr4
      ch(i - 1, k, 5) = wa4(i - 2)*dr5 - wa4(i - 1)*di5
      ch(i, k, 5) = wa4(i - 2)*di5 + wa4(i - 1)*dr5
    END DO
  END DO

  RETURN
END SUBROUTINE dradb5

SUBROUTINE dradbg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
  USE params, ONLY: DP, tPI
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: ip
  INTEGER, INTENT(IN) :: l1
  INTEGER, INTENT(IN) :: idl1
  REAL(DP), INTENT(IN) :: cc(ido, ip, l1)
  REAL(DP), INTENT(IN OUT) :: c1(ido, l1, ip)
  REAL(DP), INTENT(OUT) :: c2(idl1, ip)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, ip)
  REAL(DP), INTENT(IN OUT) :: ch2(idl1, ip)
  REAL(DP), INTENT(IN) :: wa(kdim)
  REAL(DP) :: ai1, ai2, ar1, ar1h, ar2, ar2h, arg, dc2, dcp, ds2, dsp
  INTEGER :: i, ic, idij, idp2, ik, ipp2, ipph, is
  INTEGER :: j, j2, jc, k, l, lc, nbd

  arg = tpi/(ip)
  dcp = COS(arg)
  dsp = SIN(arg)
  idp2 = ido + 2
  nbd = (ido - 1)/2
  ipp2 = ip + 2
  ipph = (ip + 1)/2
  IF (ido < l1) THEN
    DO i = 1, ido
      DO k = 1, l1
        ch(i, k, 1) = cc(i, 1, k)
      END DO
    END DO
  ELSE
    DO k = 1, l1
      DO i = 1, ido
        ch(i, k, 1) = cc(i, 1, k)
      END DO
    END DO
  END IF

  DO j = 2, ipph
    jc = ipp2 - j
    j2 = j + j
    DO k = 1, l1
      ch(1, k, j) = cc(ido, j2 - 2, k) + cc(ido, j2 - 2, k)
      ch(1, k, jc) = cc(1, j2 - 1, k) + cc(1, j2 - 1, k)
    END DO
  END DO

  IF (ido /= 1) THEN
    IF (nbd < l1) THEN
      DO j = 2, ipph
        jc = ipp2 - j
        DO i = 3, ido, 2
          ic = idp2 - i
          DO k = 1, l1
            ch(i - 1, k, j) = cc(i - 1, 2*j - 1, k) + cc(ic - 1, 2*j - 2, k)
            ch(i - 1, k, jc) = cc(i - 1, 2*j - 1, k) - cc(ic - 1, 2*j - 2, k)
            ch(i, k, j) = cc(i, 2*j - 1, k) - cc(ic, 2*j - 2, k)
            ch(i, k, jc) = cc(i, 2*j - 1, k) + cc(ic, 2*j - 2, k)
          END DO
        END DO
      END DO
    ELSE
      DO j = 2, ipph
        jc = ipp2 - j
        DO k = 1, l1
          DO i = 3, ido, 2
            ic = idp2 - i
            ch(i - 1, k, j) = cc(i - 1, 2*j - 1, k) + cc(ic - 1, 2*j - 2, k)
            ch(i - 1, k, jc) = cc(i - 1, 2*j - 1, k) - cc(ic - 1, 2*j - 2, k)
            ch(i, k, j) = cc(i, 2*j - 1, k) - cc(ic, 2*j - 2, k)
            ch(i, k, jc) = cc(i, 2*j - 1, k) + cc(ic, 2*j - 2, k)
          END DO
        END DO
      END DO
    END IF
  END IF

  ar1 = 1D0
  ai1 = 0D0
  DO l = 2, ipph
    lc = ipp2 - l
    ar1h = dcp*ar1 - dsp*ai1
    ai1 = dcp*ai1 + dsp*ar1
    ar1 = ar1h
    DO ik = 1, idl1
      c2(ik, l) = ch2(ik, 1) + ar1*ch2(ik, 2)
      c2(ik, lc) = ai1*ch2(ik, ip)
    END DO
    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1
    DO j = 3, ipph
      jc = ipp2 - j
      ar2h = dc2*ar2 - ds2*ai2
      ai2 = dc2*ai2 + ds2*ar2
      ar2 = ar2h
      DO ik = 1, idl1
        c2(ik, l) = c2(ik, l) + ar2*ch2(ik, j)
        c2(ik, lc) = c2(ik, lc) + ai2*ch2(ik, jc)
      END DO
    END DO
  END DO

  DO j = 2, ipph
    DO ik = 1, idl1
      ch2(ik, 1) = ch2(ik, 1) + ch2(ik, j)
    END DO
  END DO

  DO j = 2, ipph
    jc = ipp2 - j
    DO k = 1, l1
      ch(1, k, j) = c1(1, k, j) - c1(1, k, jc)
      ch(1, k, jc) = c1(1, k, j) + c1(1, k, jc)
    END DO
  END DO

  IF (ido == 1) RETURN
  IF (nbd < l1) THEN
    DO j = 2, ipph
      jc = ipp2 - j
      DO i = 3, ido, 2
        DO k = 1, l1
          ch(i - 1, k, j) = c1(i - 1, k, j) - c1(i, k, jc)
          ch(i - 1, k, jc) = c1(i - 1, k, j) + c1(i, k, jc)
          ch(i, k, j) = c1(i, k, j) + c1(i - 1, k, jc)
          ch(i, k, jc) = c1(i, k, j) - c1(i - 1, k, jc)
        END DO
      END DO
    END DO
  ELSE
    DO j = 2, ipph
      jc = ipp2 - j
      DO k = 1, l1
        DO i = 3, ido, 2
          ch(i - 1, k, j) = c1(i - 1, k, j) - c1(i, k, jc)
          ch(i - 1, k, jc) = c1(i - 1, k, j) + c1(i, k, jc)
          ch(i, k, j) = c1(i, k, j) + c1(i - 1, k, jc)
          ch(i, k, jc) = c1(i, k, j) - c1(i - 1, k, jc)
        END DO
      END DO
    END DO
  END IF

  DO ik = 1, idl1
    c2(ik, 1) = ch2(ik, 1)
  END DO

  DO j = 2, ip
    DO k = 1, l1
      c1(1, k, j) = ch(1, k, j)
    END DO
  END DO

  IF (nbd > l1) THEN
    is = -ido
    DO j = 2, ip
      is = is + ido
      DO k = 1, l1
        idij = is
        DO i = 3, ido, 2
          idij = idij + 2
          c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - wa(idij)*ch(i, k, j)
          c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij)*ch(i - 1, k, j)
        END DO
      END DO
    END DO
  ELSE
    is = -ido
    DO j = 2, ip
      is = is + ido
      idij = is
      DO i = 3, ido, 2
        idij = idij + 2
        DO k = 1, l1
          c1(i - 1, k, j) = wa(idij - 1)*ch(i - 1, k, j) - wa(idij)*ch(i, k, j)
          c1(i, k, j) = wa(idij - 1)*ch(i, k, j) + wa(idij)*ch(i - 1, k, j)
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE dradbg

SUBROUTINE dradf2(ido, l1, cc, ch, wa1)
  USE params, ONLY: DP
  IMPLICIT NONE
  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, l1, 2)
  REAL(DP), INTENT(OUT) :: ch(ido, 2, l1)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP) :: ti2, tr2
  INTEGER :: i, ic, idp2, k

  DO k = 1, l1
    ch(1, 1, k) = cc(1, k, 1) + cc(1, k, 2)
    ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 2)
  END DO
  STOP
  IF (ido < 2) THEN
    RETURN
  ELSE IF (ido /= 2) THEN
    idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        tr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
        ti2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
        ch(i, 1, k) = cc(i, k, 1) + ti2
        ch(ic, 2, k) = ti2 - cc(i, k, 1)
        ch(i - 1, 1, k) = cc(i - 1, k, 1) + tr2
        ch(ic - 1, 2, k) = cc(i - 1, k, 1) - tr2
      END DO
    END DO

    IF (MOD(ido, 2) == 1) RETURN
  END IF

  DO k = 1, l1
    ch(1, 2, k) = -cc(ido, k, 2)
    ch(ido, 1, k) = cc(ido, k, 1)
  END DO

  RETURN
END SUBROUTINE dradf2

SUBROUTINE dradf3(ido, l1, cc, ch, wa1, wa2)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, l1, 3)
  REAL(DP), INTENT(OUT) :: ch(ido, 3, l1)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP) :: ci2, cr2, di2, di3, dr2, dr3, taui, taur, ti2, ti3, tr2, tr3
  INTEGER :: i, ic, idp2, k
  DATA taur/-0.5D0/
  DATA taui/0.86602540378443864676372317075293618D0/

  DO k = 1, l1
    cr2 = cc(1, k, 2) + cc(1, k, 3)
    ch(1, 1, k) = cc(1, k, 1) + cr2
    ch(1, 3, k) = taui*(cc(1, k, 3) - cc(1, k, 2))
    ch(ido, 2, k) = cc(1, k, 1) + taur*cr2
  END DO

  IF (ido == 1) RETURN
  idp2 = ido + 2
  DO k = 1, l1
    DO i = 3, ido, 2
      ic = idp2 - i
      dr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
      di2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
      dr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
      di3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
      cr2 = dr2 + dr3
      ci2 = di2 + di3
      ch(i - 1, 1, k) = cc(i - 1, k, 1) + cr2
      ch(i, 1, k) = cc(i, k, 1) + ci2
      tr2 = cc(i - 1, k, 1) + taur*cr2
      ti2 = cc(i, k, 1) + taur*ci2
      tr3 = taui*(di2 - di3)
      ti3 = taui*(dr3 - dr2)
      ch(i - 1, 3, k) = tr2 + tr3
      ch(ic - 1, 2, k) = tr2 - tr3
      ch(i, 3, k) = ti2 + ti3
      ch(ic, 2, k) = ti3 - ti2
    END DO
  END DO

  RETURN
END SUBROUTINE dradf3

SUBROUTINE dradf4(ido, l1, cc, ch, wa1, wa2, wa3)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, l1, 4)
  REAL(DP), INTENT(OUT) :: ch(ido, 4, l1)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP), INTENT(IN) :: wa3(kdim)
  REAL(DP) :: ci2, ci3, ci4, cr2, cr3, cr4, hsqt2, ti1, ti2, ti3, &
              ti4, tr1, tr2, tr3, tr4
  INTEGER :: i, ic, idp2, k
  DATA hsqt2/.70710678118654752440084436210485D0/

  DO k = 1, l1
    tr1 = cc(1, k, 2) + cc(1, k, 4)
    tr2 = cc(1, k, 1) + cc(1, k, 3)
    ch(1, 1, k) = tr1 + tr2
    ch(ido, 4, k) = tr2 - tr1
    ch(ido, 2, k) = cc(1, k, 1) - cc(1, k, 3)
    ch(1, 3, k) = cc(1, k, 4) - cc(1, k, 2)
  END DO

  IF (ido < 2) THEN
    RETURN
  ELSE IF (ido /= 2) THEN
    idp2 = ido + 2
    DO k = 1, l1
      DO i = 3, ido, 2
        ic = idp2 - i
        cr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
        ci2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
        cr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
        ci3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
        cr4 = wa3(i - 2)*cc(i - 1, k, 4) + wa3(i - 1)*cc(i, k, 4)
        ci4 = wa3(i - 2)*cc(i, k, 4) - wa3(i - 1)*cc(i - 1, k, 4)
        tr1 = cr2 + cr4
        tr4 = cr4 - cr2
        ti1 = ci2 + ci4
        ti4 = ci2 - ci4
        ti2 = cc(i, k, 1) + ci3
        ti3 = cc(i, k, 1) - ci3
        tr2 = cc(i - 1, k, 1) + cr3
        tr3 = cc(i - 1, k, 1) - cr3
        ch(i - 1, 1, k) = tr1 + tr2
        ch(ic - 1, 4, k) = tr2 - tr1
        ch(i, 1, k) = ti1 + ti2
        ch(ic, 4, k) = ti1 - ti2
        ch(i - 1, 3, k) = ti4 + tr3
        ch(ic - 1, 2, k) = tr3 - ti4
        ch(i, 3, k) = tr4 + ti3
        ch(ic, 2, k) = tr4 - ti3
      END DO
    END DO
    IF (MOD(ido, 2) == 1) RETURN
  END IF

  DO k = 1, l1
    ti1 = -hsqt2*(cc(ido, k, 2) + cc(ido, k, 4))
    tr1 = hsqt2*(cc(ido, k, 2) - cc(ido, k, 4))
    ch(ido, 1, k) = tr1 + cc(ido, k, 1)
    ch(ido, 3, k) = cc(ido, k, 1) - tr1
    ch(1, 2, k) = ti1 - cc(ido, k, 3)
    ch(1, 4, k) = ti1 + cc(ido, k, 3)
  END DO

  RETURN
END SUBROUTINE dradf4

SUBROUTINE dradf5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
  USE params, ONLY: DP
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: l1
  REAL(DP), INTENT(IN) :: cc(ido, l1, 5)
  REAL(DP), INTENT(OUT) :: ch(ido, 5, l1)
  REAL(DP), INTENT(IN) :: wa1(kdim)
  REAL(DP), INTENT(IN) :: wa2(kdim)
  REAL(DP), INTENT(IN) :: wa3(kdim)
  REAL(DP), INTENT(IN) :: wa4(kdim)
  REAL(DP) :: ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, &
              di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3, ti4, &
              ti5, tr11, tr12, tr2, tr3, tr4, tr5
  INTEGER :: i, ic, idp2, k

  DATA tr11/0.30901699437494742410229341718281906D0/
  DATA ti11/0.95105651629515357211643933337938214D0/
  DATA tr12/-0.80901699437494742410229341718281906D0/
  DATA ti12/0.58778525229247312916870595463907277D0/

  DO k = 1, l1
    cr2 = cc(1, k, 5) + cc(1, k, 2)
    ci5 = cc(1, k, 5) - cc(1, k, 2)
    cr3 = cc(1, k, 4) + cc(1, k, 3)
    ci4 = cc(1, k, 4) - cc(1, k, 3)
    ch(1, 1, k) = cc(1, k, 1) + cr2 + cr3
    ch(ido, 2, k) = cc(1, k, 1) + tr11*cr2 + tr12*cr3
    ch(1, 3, k) = ti11*ci5 + ti12*ci4
    ch(ido, 4, k) = cc(1, k, 1) + tr12*cr2 + tr11*cr3
    ch(1, 5, k) = ti12*ci5 - ti11*ci4
  END DO

  IF (ido == 1) RETURN
  idp2 = ido + 2
  DO k = 1, l1
    DO i = 3, ido, 2
      ic = idp2 - i
      dr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
      di2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
      dr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
      di3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
      dr4 = wa3(i - 2)*cc(i - 1, k, 4) + wa3(i - 1)*cc(i, k, 4)
      di4 = wa3(i - 2)*cc(i, k, 4) - wa3(i - 1)*cc(i - 1, k, 4)
      dr5 = wa4(i - 2)*cc(i - 1, k, 5) + wa4(i - 1)*cc(i, k, 5)
      di5 = wa4(i - 2)*cc(i, k, 5) - wa4(i - 1)*cc(i - 1, k, 5)
      cr2 = dr2 + dr5
      ci5 = dr5 - dr2
      cr5 = di2 - di5
      ci2 = di2 + di5
      cr3 = dr3 + dr4
      ci4 = dr4 - dr3
      cr4 = di3 - di4
      ci3 = di3 + di4
      ch(i - 1, 1, k) = cc(i - 1, k, 1) + cr2 + cr3
      ch(i, 1, k) = cc(i, k, 1) + ci2 + ci3
      tr2 = cc(i - 1, k, 1) + tr11*cr2 + tr12*cr3
      ti2 = cc(i, k, 1) + tr11*ci2 + tr12*ci3
      tr3 = cc(i - 1, k, 1) + tr12*cr2 + tr11*cr3
      ti3 = cc(i, k, 1) + tr12*ci2 + tr11*ci3
      tr5 = ti11*cr5 + ti12*cr4
      ti5 = ti11*ci5 + ti12*ci4
      tr4 = ti12*cr5 - ti11*cr4
      ti4 = ti12*ci5 - ti11*ci4
      ch(i - 1, 3, k) = tr2 + tr5
      ch(ic - 1, 2, k) = tr2 - tr5
      ch(i, 3, k) = ti2 + ti5
      ch(ic, 2, k) = ti5 - ti2
      ch(i - 1, 5, k) = tr3 + tr4
      ch(ic - 1, 4, k) = tr3 - tr4
      ch(i, 5, k) = ti3 + ti4
      ch(ic, 4, k) = ti4 - ti3
    END DO
  END DO

  RETURN
END SUBROUTINE dradf5

SUBROUTINE dradfg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
  USE params, ONLY: DP, tPI
  IMPLICIT NONE

  INTEGER, PARAMETER :: kdim = 512

  INTEGER, INTENT(IN) :: ido
  INTEGER, INTENT(IN) :: ip
  INTEGER, INTENT(IN) :: l1
  INTEGER, INTENT(IN) :: idl1
  REAL(DP), INTENT(OUT) :: cc(ido, ip, l1)
  REAL(DP), INTENT(IN OUT) :: c1(ido, l1, ip)
  REAL(DP), INTENT(IN OUT) :: c2(idl1, ip)
  REAL(DP), INTENT(OUT) :: ch(ido, l1, ip)
  REAL(DP), INTENT(OUT) :: ch2(idl1, ip)
  REAL(DP), INTENT(IN) :: wa(kdim)
  REAL(DP) :: ai1, ai2, ar1, ar1h, ar2, ar2h, arg, dc2, dcp, ds2, dsp
  INTEGER :: i, ic, idij, idp2, ik, ipp2, ipph, is
  INTEGER :: j, j2, jc, k, l, lc, nbd

  arg = tpi/(ip)
  dcp = COS(arg)
  dsp = SIN(arg)
  ipph = (ip + 1)/2
  ipp2 = ip + 2
  idp2 = ido + 2
  nbd = (ido - 1)/2
  IF (ido == 1) THEN
    DO ik = 1, idl1
      c2(ik, 1) = ch2(ik, 1)
    END DO
  ELSE
    DO ik = 1, idl1
      ch2(ik, 1) = c2(ik, 1)
    END DO
    DO j = 2, ip
      DO k = 1, l1
        ch(1, k, j) = c1(1, k, j)
      END DO
    END DO

    IF (nbd > l1) THEN
      is = -ido
      DO j = 2, ip
        is = is + ido
        DO k = 1, l1
          idij = is
          DO i = 3, ido, 2
            idij = idij + 2
            ch(i - 1, k, j) = wa(idij - 1)*c1(i - 1, k, j) + wa(idij)*c1(i, k, j)
            ch(i, k, j) = wa(idij - 1)*c1(i, k, j) - wa(idij)*c1(i - 1, k, j)
          END DO
        END DO
      END DO
    ELSE
      is = -ido
      DO j = 2, ip
        is = is + ido
        idij = is
        DO i = 3, ido, 2
          idij = idij + 2
          DO k = 1, l1
            ch(i - 1, k, j) = wa(idij - 1)*c1(i - 1, k, j) + wa(idij)*c1(i, k, j)
            ch(i, k, j) = wa(idij - 1)*c1(i, k, j) - wa(idij)*c1(i - 1, k, j)
          END DO
        END DO
      END DO
    END IF

    IF (nbd < l1) THEN
      DO j = 2, ipph
        jc = ipp2 - j
        DO i = 3, ido, 2
          DO k = 1, l1
            c1(i - 1, k, j) = ch(i - 1, k, j) + ch(i - 1, k, jc)
            c1(i - 1, k, jc) = ch(i, k, j) - ch(i, k, jc)
            c1(i, k, j) = ch(i, k, j) + ch(i, k, jc)
            c1(i, k, jc) = ch(i - 1, k, jc) - ch(i - 1, k, j)
          END DO
        END DO
      END DO
    ELSE
      DO j = 2, ipph
        jc = ipp2 - j
        DO k = 1, l1
          DO i = 3, ido, 2
            c1(i - 1, k, j) = ch(i - 1, k, j) + ch(i - 1, k, jc)
            c1(i - 1, k, jc) = ch(i, k, j) - ch(i, k, jc)
            c1(i, k, j) = ch(i, k, j) + ch(i, k, jc)
            c1(i, k, jc) = ch(i - 1, k, jc) - ch(i - 1, k, j)
          END DO
        END DO
      END DO
    END IF
  END IF

  DO j = 2, ipph
    jc = ipp2 - j
    DO k = 1, l1
      c1(1, k, j) = ch(1, k, j) + ch(1, k, jc)
      c1(1, k, jc) = ch(1, k, jc) - ch(1, k, j)
    END DO
  END DO

  ar1 = 1.d0
  ai1 = 0.d0
  DO l = 2, ipph
    lc = ipp2 - l
    ar1h = dcp*ar1 - dsp*ai1
    ai1 = dcp*ai1 + dsp*ar1
    ar1 = ar1h
    DO ik = 1, idl1
      ch2(ik, l) = c2(ik, 1) + ar1*c2(ik, 2)
      ch2(ik, lc) = ai1*c2(ik, ip)
    END DO
    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1
    DO j = 3, ipph
      jc = ipp2 - j
      ar2h = dc2*ar2 - ds2*ai2
      ai2 = dc2*ai2 + ds2*ar2
      ar2 = ar2h
      DO ik = 1, idl1
        ch2(ik, l) = ch2(ik, l) + ar2*c2(ik, j)
        ch2(ik, lc) = ch2(ik, lc) + ai2*c2(ik, jc)
      END DO
    END DO
  END DO

  DO j = 2, ipph
    DO ik = 1, idl1
      ch2(ik, 1) = ch2(ik, 1) + c2(ik, j)
    END DO
  END DO

  IF (ido < l1) THEN
    DO i = 1, ido
      DO k = 1, l1
        cc(i, 1, k) = ch(i, k, 1)
      END DO
    END DO
  ELSE
    DO k = 1, l1
      DO i = 1, ido
        cc(i, 1, k) = ch(i, k, 1)
      END DO
    END DO
  END IF

  DO j = 2, ipph
    jc = ipp2 - j
    j2 = j + j
    DO k = 1, l1
      cc(ido, j2 - 2, k) = ch(1, k, j)
      cc(1, j2 - 1, k) = ch(1, k, jc)
    END DO
  END DO

  IF (ido == 1) RETURN

  IF (nbd < l1) THEN
    DO j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      DO i = 3, ido, 2
        ic = idp2 - i
        DO k = 1, l1
          cc(i - 1, j2 - 1, k) = ch(i - 1, k, j) + ch(i - 1, k, jc)
          cc(ic - 1, j2 - 2, k) = ch(i - 1, k, j) - ch(i - 1, k, jc)
          cc(i, j2 - 1, k) = ch(i, k, j) + ch(i, k, jc)
          cc(ic, j2 - 2, k) = ch(i, k, jc) - ch(i, k, j)
        END DO
      END DO
    END DO
    RETURN
  ELSE
    DO j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      DO k = 1, l1
        DO i = 3, ido, 2
          ic = idp2 - i
          cc(i - 1, j2 - 1, k) = ch(i - 1, k, j) + ch(i - 1, k, jc)
          cc(ic - 1, j2 - 2, k) = ch(i - 1, k, j) - ch(i - 1, k, jc)
          cc(i, j2 - 1, k) = ch(i, k, j) + ch(i, k, jc)
          cc(ic, j2 - 2, k) = ch(i, k, jc) - ch(i, k, j)
        END DO
      END DO
    END DO
    RETURN
  END IF

END SUBROUTINE dradfg

INCLUDE "fftpack2.F90"
